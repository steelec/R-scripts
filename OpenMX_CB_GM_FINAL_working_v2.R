## -----------------------------------------------------------------------------
# 2016-04-13
# Christopher Steele
# Adapted from openMX code available from conference workshops (Sarah Medland)
# Matrix design MZ vs DZ vs non-twin siblings design, treating non-twin sibs as
# DZ twins
# - Analysis loops across columns of data as defined in allVars
# - Tests the ACE and CE models for each, tests for sig difference between the two
# - Controls for age and sex, and optional controlVar (set with booleans)
#
# Things to watch out for:
# 1) twin (MZ/DZ) data MUST occupy the first two columns of family-wise data
#    this is defined by the ordering
# 2) changing how M and F are coded changes the confidence intervals a bit (weird...)
# 3) 2 methods to use control var, TYPE2 is likely the most accurate, but slower :-(
# 4) used NPSOL optimiser because it allows me to calculate the CIs of the estimates
#    ends up with some "red" mxflags, but same results as default optimiser (no red)
# 5) only twins that have ages that are exactly the same are included in these analyses
## -----------------------------------------------------------------------------

## clean workspace
rm(list = ls())

require(OpenMx)
mxOption(NULL, "Default optimizer", "NPSOL") #changes the optimiser, allows me to calc the upper/lower bounds (CI) of genetic effect - the estimates seem to stay the same, using this one!
#mxOption(NULL, "Number of Threads", detectCores() - 1)

## set booleans for how we want to process this data
USE_CONTROL_VAR_LM =                   FALSE #regress out a control variable using lm (set variable as controlVar <- ???? )
USE_CONTROL_VAR_MX =                   FALSE #similar to above, but includes the variable in the openMX model (similar results)
                                            # does not apply to COMPARE_TWINS_DIRECTLY=TRUE, since I did not implement it
                                            # use LM version for similar results
USE_TWINS_ONLY  =                      FALSE #only uses the MZ and DZ twin data (produces the same results as with COMPARE_TWINS_DIRECTLY)
SKIP_ACE_COMPUTATION =                 FALSE  # go straight to calculating the summary stats
SKIP_EXTRA_PLOTTING =                  FALSE # this also does something :-)
COMPARE_TWINS_DIRECTLY_MX =            FALSE #uses a different approach within MX, concurrently outputs results to res_twins_only (no sig testing)
REMOVE_FAMILY_MEMBERS_CORR_PLOTS =     FALSE # creates correlations and plots when only a single individual from each family is selected (randomly)
USE_FAMILY_SUBSET =                    FALSE # remove some families from the analysis to check how stable it is (randomly)
NUM_FAMILES_REMOVE=                    30 #number of families to remove


## choose your input data
bx_dir='/home/chris/Documents/Dropbox/Projects/Working/CoBrA_Lab/TractREC/processing/bx/'
img_out_dir='/home/chris/Documents/Dropbox/Projects/Working/CoBrA_Lab/TractREC/processing/bx/images/'
#CB_lobs_S500='2015_12_HBM_HCP_allSibs_CB_LR_volume_FA_corrected_0p35_restr_unrestr_CORRECTED_SUBSET.csv'
#HcpData<- read.csv(paste(bx_dir,CB_lobs_S500,sep=""),header=T,stringsAsFactors = FALSE)
#all_structures_noFA_thresh_S500="2016_03_HCP_allSibs_all_labels_volume_restr_unrestr.csv"
#all_structures_SEJAL_S500='2016_03_HCP_allSibs_SEJAL_all_labels_volume_restr_unrestr.csv'
#HcpData<- read.csv(paste(bx_dir,all_structures_SEJAL_S500,sep=""),header=T,stringsAsFactors = FALSE)
#library(plyr)
#rename(HcpData,c("Sex","Gender"))
all_structures_S500 = "2016_03_HCP_allSibs_all_labels_volume_FA_corrected_0p35_restr_unrestr.csv"
all_volumes_all_metrics_S500_FA_fname = "/home/chris/Documents/Dropbox/Projects/Working/CoBrA_Lab/TractREC/processing/bx/S500/2016_03_HCP_FA_allLabels_FA_0p35_all.csv"
#HcpData<- read.csv(paste(bx_dir,all_structures_noFA_thresh_S500,sep=""),header=T,stringsAsFactors = FALSE)
HcpData<- read.csv(paste(bx_dir,all_structures_S500,sep=""),header=T,stringsAsFactors = FALSE)

## update Zygosity with version corrected for families that included a single MZ twin
HcpData$MZ_NotMZ_NotTwin <- HcpData$Zygosity #record what the individual was in a new var
HcpData$Zygosity <- HcpData$Recoded_Zygosity #XXX recoded to correct for families that only include a single MZ twin

## reorder the siblings such that the 1st and 2nd are always the twins
## XXX this is because the covariance matrix REQUIRES twins first
## TWINS (MZ or DZ) must always be in the 1st two columns of the family, so ordered 1 and 2 here
for (mother in unique(HcpData$Mother_ID)) {
  ## use "-" to reverse the zygosity data so that rank is decreasing
  ## add twin_stat*2 to push the twins higher even if their zygosity is .5 (DZ twins), so they will come out first
  sibs<-unique(HcpData$ID[HcpData$Mother_ID==mother])
  twin_stat<-HcpData$Twin_Stat[HcpData$Mother_ID==mother]=="Twin" #true or false that each of these is a twin
  ordered_zygosity=rank(-(HcpData$Zygosity[HcpData$ID %in% sibs] + twin_stat*rep(2,length(twin_stat))),ties.method="first") #rank by their zygosities, decreasing
  HcpData$sibling_number2[HcpData$Mother_ID==mother]<-ordered_zygosity
  
  #calculate the sibling count and Twins_Count for potential filtering later
  twin_sibs <- subset(HcpData,Mother_ID==mother & Twin_Stat == "Twin")$ID
  HcpData$Sibling_Count[HcpData$Mother_ID==mother] <- length(sibs)
  HcpData$Twins_Count[HcpData$Mother_ID==mother] <- length(twin_sibs)
}

#remove the individuals who are twins, but their twin sib is not in the dataset
HcpData <- subset(HcpData,Twins_Count != 1) #singleton twins removed
HcpData <- subset(HcpData,Twins_Count != 0) #families with no twins removed
#HcpData <- subset(HcpData,Sibling_Count > 1) #also need to have more than one sibling?

## create an ordered level for gender
HcpData$Gender_Txt <- HcpData$Gender
HcpData$Gender <- mxFactor(HcpData$Gender,levels=c('M','F'))
HcpData$Gender <- as.numeric(HcpData$Gender)

## zygosity now taken care of with a recoded variable to account for datasets where twins were split up
#tmp_zyg <- mxFactor(HcpData$Zygosity,labels=c(1,.5,.5),levels=c('MZ','NotMZ','NotTwin'),collapse=TRUE)
#HcpData$Zygosity <- as.numeric(levels(tmp_zyg)[tmp_zyg])


## select the columns from the dataframe that you are interested in running the analysis on

#allVars <- c('label_031_L_CA1_volume','label_032_L_subiculum_volume','label_131_R_CA1_volume','label_132_R_subiculum_volume')
#allVars <- c('label_002_L_III_volume', 'label_003_L_IV_volume', 'label_004_L_V_volume', 'label_005_L_VI_volume', 'label_006_L_Crus_I_volume', 'label_007_L_Crus_II_volume', 'label_008_L_VIIB_volume', 'label_009_L_VIIIA_volume', 'label_010_L_VIIIB_volume', 'label_011_L_IX_volume', 'label_012_L_X_volume', 'label_102_R_III_volume', 'label_103_R_IV_volume', 'label_104_R_V_volume', 'label_105_R_VI_volume', 'label_106_R_Crus_I_volume', 'label_107_R_Crus_II_volume', 'label_108_R_VIIB_volume', 'label_109_R_VIIIA_volume', 'label_110_R_VIIIB_volume', 'label_111_R_IX_volume', 'label_112_R_X_volume', 'label_001_Vermal_I_II_volume', 'label_014_Vermal_III_volume', 'label_015_Vermal_IV_volume', 'label_016_Vermal_V_volume', 'label_017_Vermal_VI_volume', 'label_018_Vermal_VIIA_volume', 'label_019_Vermal_VIIB_volume', 'label_020_Vermal_VIIIA_volume', 'label_021_Vermal_VIIIB_volume', 'label_022_Vermal_IX_volume', 'label_023_Vermal_X_volume','GM_total_volume')
#allVars <- c('label_002_L_III_volume', 'label_003_L_IV_volume', 'label_004_L_V_volume', 'label_005_L_VI_volume', 'label_006_L_Crus_I_volume', 'label_007_L_Crus_II_volume', 'label_008_L_VIIB_volume', 'label_009_L_VIIIA_volume', 'label_010_L_VIIIB_volume', 'label_011_L_IX_volume', 'label_012_L_X_volume', 'label_102_R_III_volume', 'label_103_R_IV_volume', 'label_104_R_V_volume', 'label_105_R_VI_volume', 'label_106_R_Crus_I_volume', 'label_107_R_Crus_II_volume', 'label_108_R_VIIB_volume', 'label_109_R_VIIIA_volume', 'label_110_R_VIIIB_volume', 'label_111_R_IX_volume', 'label_112_R_X_volume', 'label_001_Vermal_I_II_volume', 'label_014_Vermal_III_volume', 'label_015_Vermal_IV_volume', 'label_016_Vermal_V_volume', 'label_017_Vermal_VI_volume', 'label_018_Vermal_VIIA_volume', 'label_019_Vermal_VIIB_volume', 'label_020_Vermal_VIIIA_volume', 'label_021_Vermal_VIIIB_volume', 'label_022_Vermal_IX_volume', 'label_023_Vermal_X_volume')
allVars <- colnames(HcpData[grep("label_",colnames(HcpData))]) #select all column names that have "label_" in them
allVars <- allVars[allVars != "label_file"]
CBVars <-allVars[1:33] #use this to calculate the sum of all CB GM volume
allVars<-CBVars #also only care about them right now
HcpData$label_000_CB_GM_volume<-rowSums(HcpData[CBVars])
#allVars<-c("label_031_L_CA1_volume","label_000_CB_GM_volume","Height","Weight","MZ_NotMZ_NotTwin","Gender_Txt")
#allVars<-allVars[1:6]
## set your control variable
controlVar_txt<-"label_000_CB_GM_volume"
#controlVar_txt<-"label_080_ITV_volume"
#allVars <- colnames(HcpData)
#allVars<-c(allVars,"label_000_CB_GM_volume","Height","Weight","Handedness","MZ_NotMZ_NotTwin","Gender_Txt")
#allVars<-c("label_031_L_CA1_volume","Height","Weight","Handedness","MZ_NotMZ_NotTwin","Gender_Txt")
allVars <- unique(c(allVars, controlVar_txt)) #make sure that we don't list anything twice!

## create a subset of our data based on the boolean flags that we have set
# NOTE: twins are selected ONLY if the age is exactly the same
if (USE_TWINS_ONLY){
  twinned_ids=c()
  HcpData_twin_subset <- subset(HcpData, Twin_Stat == "Twin")
  ## we also need to clear out any individuals who are "Twin" but their twin is not in the dbase
  for (mother in unique(HcpData_twin_subset$Mother_ID)) {
    sibs<-unique(HcpData_twin_subset$ID[HcpData_twin_subset$Mother_ID==mother])
    #print(length(sibs))
    if (length(sibs)==2) {
      #if they have the same age, then they are twins!
      if ((HcpData_twin_subset$Age_in_Yrs[HcpData_twin_subset$ID==sibs[1]]==HcpData_twin_subset$Age_in_Yrs[HcpData_twin_subset$ID==sibs[2]])) {
        twinned_ids <- c(twinned_ids,sibs)
      }
    } else if (length(sibs)> 2) {
      print("WHOA! There are more than two twins here...!")
      for (sib in sibs) {
        sibs_sub=sibs[sibs!=sib] #subset of siblings
        for (other_sib in sibs_sub) {
          if (HcpData_twin_subset$Age_in_Yrs[HcpData_twin_subset$ID==sib] == HcpData_twin_subset$Age_in_Yrs[HcpData_twin_subset$ID==other_sib]) {
            twinned_ids <- c(twinned_ids,sib,other_sib)
            print(sib)
          } 
        }
      }
    }
  }
  print(length(twinned_ids))
  twinned_ids <- unique(twinned_ids) #double entries with the last one, though it should never get there
  print(length(twinned_ids))
  
  
  HcpData_twin_subset<- HcpData_twin_subset[HcpData_twin_subset$ID %in% twinned_ids, ]
  allVars <- unique(c('ID','Age_in_Yrs','Gender','Mother_ID','sibling_number2','Zygosity',allVars)) #don't duplicate
  HcpData_subset <- HcpData_twin_subset[allVars]
} else {
  # TODO: XXX : check in here and correct for individuals who don't belong!
  allVars <- unique(c('ID','Age_in_Yrs','Gender','Mother_ID','sibling_number2','Zygosity',allVars)) #don't duplicate
  HcpData_subset <- HcpData[allVars]
}

#remove Mother_ID and sibling_number2 from allVars since they will be ordering our matrix in wide format
allVars <- subset(allVars,!allVars %in% c("Mother_ID","sibling_number2",'ID','Age_in_Yrs','Gender','Zygosity','MZ_NotMZ_NotTwin'))

temp <- quote(controlVar_txt)
controlVar <- HcpData_subset[,eval(temp)]
rm(temp)
#HcpData$Menstrual_DaysSinceLast

HcpData_subset$Gender = HcpData_subset$Gender-1 #now 0 and 1

# XXX in progress to remove a subset of the mothers from the sample to assess stability of the estimate
unique_mothers=unique(HcpData$Mother_ID)
removed_mothers=sample(unique_mothers,NUM_FAMILES_REMOVE)
if (USE_FAMILY_SUBSET){
  HcpData_subset <- subset(HcpData_subset,!Mother_ID %in% removed_mothers)
}

#if we want to partial out the effect of some var prior to fitting, we can do that here
# just before converting to the family wide layout
if (USE_CONTROL_VAR_LM) {
  HcpData_subset_orig <- HcpData_subset # create a copy, since we will overwrite the current one
  for (theVarName in allVars) {
    if (!is.character(HcpData_subset[theVarName][,])) { #if the data is just characters, we skip it
      lm_fit <- lm(HcpData_subset[theVarName][,] ~ controlVar)
      HcpData_subset[theVarName] <- resid(lm_fit)
    }
  }
}

## reformat the data into families (each row = 1 family), ordered by sibling number
#  such that MZ is first (if present), then DZ (if present), then sibs (if present)
fam<-reshape(HcpData_subset,idvar="Mother_ID",timevar="sibling_number2",direction="wide")
names(fam)<-gsub("[.]","_",names(fam))
fam[is.na(fam)] <- 0 ## set values to zero where there is no data

## create a matrix to hold the results for each label
res <- matrix(data=NA, ncol=18, nrow=length(allVars))
colnames(res) <- c("estVA", "estVC", "estVE", "estVP", "estPropVA", "estPropVC", "estPropVE", "estVC_mCE", "estVE_mCE", "estVP_mCE", "estPropVC_mCE", "estPropVE_mCE", "LL_ACE", "LL_CE", "LRT_ACE_CE","ACE_vs_CE_p","A_CI_lower","A_CI_upper")
res_twins_only <- matrix(data=NA, ncol=8, nrow=length(allVars))
colnames(res_twins_only) <- c("estVA", "estVC", "estVE", "estVP", "estPropVA", "estPropVC", "estPropVE","LL_ACE")
allVars_txt<-gsub("_volume","",allVars)
#rownames(res) <- allVars
rownames(res) <- gsub("_"," ",gsub("^.*?_","",gsub("^.*?_","",allVars_txt)))
rownames(res_twins_only) <- gsub("_"," ",gsub("^.*?_","",gsub("^.*?_","",allVars_txt)))

## loop over each variable in allVars
if (!SKIP_ACE_COMPUTATION) {
  row_idx=0
  for(theVarName in allVars){
    row_idx=row_idx+1
    print(theVarName)
    if (!is.character(HcpData_subset[theVarName][,])) { #if the data is just characters, we skip it
      # Select Variables for Analysis
      nv <- 1 # number of variables
      #ntv <- nv*4 # tot num vars un the table #see below
      
      ## this is how to extract the values from the variable name that we have selected on the loop
      ## ends up setting the total number of vars depending on if we are only selecting twins or not
      if (!USE_TWINS_ONLY){
        selVars <- paste(theVarName,c(rep(1,nv),rep(2,nv),rep(3,nv),rep(4,nv)),sep="_")
      } else {
        selVars <- paste(theVarName,c(rep(1,nv),rep(2,nv)),sep="_")
      }
      
      temp <- quote(theVarName) #this and the following generate a list of names of the theVar ending in _1:4
      theVar <- HcpData_subset[,eval(temp)] # for use in indexing the fam dataframe
      
      ntv <- length(selVars) # number of total variables (b/c we have a max of 4 siblings)
      
      #calc mean and Path coeff
      svMe <- mean(theVar) # start value for the means
      svPa <- sqrt(var(theVar)/3) # start value for path coefficients (sqrt(variance/#ofpaths))
      #svPa <- 1.0 #XXX SET TO ONE TO TEST
      
      
      # ACE Model
      
      # Matrices declared to store a, c, and e Path Coefficients
      pathA <- mxMatrix( type="Lower", nrow=nv, ncol=nv,free=TRUE, values=svPa, label="a11", name="a" )
      pathC <- mxMatrix( type="Lower", nrow=nv, ncol=nv,free=TRUE, values=svPa, label="c11", name="c" )
      pathE <- mxMatrix( type="Lower", nrow=nv, ncol=nv,free=TRUE, values=svPa, label="e11", name="e" )
      
      # Matrices generated to hold A, C, and E computed Variance Components
      covA <- mxAlgebra( expression=a %*% t(a), name="A" )
      covC <- mxAlgebra( expression=c %*% t(c), name="C" )
      covE <- mxAlgebra( expression=e %*% t(e), name="E" )
      
      # Algebra to compute total variances
      covP <- mxAlgebra( expression=A+C+E, name="V" )
      
      # Algebra for expected Mean, beta sex and age and Variance/Covariance Matrices in MZ & DZ twins
      
      intercept <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE,values= svMe, label="mean", name="Mean" )
      
      
      # Matrix for moderating/interacting variable (sex)
      if (!USE_TWINS_ONLY){
        defSex <- mxMatrix( type="Full", nrow=1, ncol=4, free=FALSE, labels=c("data.Gender_1","data.Gender_2","data.Gender_3","data.Gender_4"), name="Sex")
      } else {
        defSex <- mxMatrix( type="Full", nrow=1, ncol=2, free=FALSE, labels=c("data.Gender_1","data.Gender_2"), name="Sex")
      }
      # Matrices declared to store linear Coefficients for covariate
      B_Sex  <- mxMatrix( type="Full", nrow=1, ncol=1, free=TRUE, values=0, label="betaSex", name="bSex" )
      meanSex <- mxAlgebra(  bSex%*%Sex, name="SexR")
      
      
      # Matrix for moderating/interacting variable (age)
      if (!USE_TWINS_ONLY){
        defAge    <- mxMatrix( type="Full", nrow=1, ncol=4, free=FALSE, labels=c("data.Age_in_Yrs_1","data.Age_in_Yrs_2","data.Age_in_Yrs_3","data.Age_in_Yrs_4"), name="Age")
      } else {
        defAge    <- mxMatrix( type="Full", nrow=1, ncol=2, free=FALSE, labels=c("data.Age_in_Yrs_1","data.Age_in_Yrs_2"), name="Age")
      }
      # Matrices declared to store linear Coefficients for covariate
      B_Age     <- mxMatrix( type="Full", nrow=1, ncol=1, free=TRUE, values= 0, label="betaAge", name="bAge" )
      meanAge   <- mxAlgebra(  bAge%*%Age, name="AgeR")
      
      if (USE_CONTROL_VAR_MX) {
        controlVars_txt <- paste('data',paste(controlVar_txt,c(rep(1,nv),rep(2,nv),rep(3,nv),rep(4,nv)),sep="_"),sep=".")
        # Matrix for moderating/interacting variable (control covariate)
        if (!USE_TWINS_ONLY){
          defControlCov    <- mxMatrix( type="Full", nrow=1, ncol=4, free=FALSE, labels=controlVars_txt, name="ControlCov")
        } else {
          controlVars_txt <- paste('data',paste(controlVar_txt,c(rep(1,nv),rep(2,nv)),sep="_"),sep=".")
          defControlCov    <- mxMatrix( type="Full", nrow=1, ncol=2, free=FALSE, labels=controlVars_txt[1:2], name="ControlCov")
        }
        # Matrices declared to store linear Coefficients for additional covariate
        B_ControlCov     <- mxMatrix( type="Full", nrow=1, ncol=1, free=TRUE, values= 0, label="betaControlCov", name="bControlCov" )
        meanControlCov   <- mxAlgebra(  bControlCov%*%ControlCov, name="ControlCovR")
      }
      
      if (USE_CONTROL_VAR_MX) {  
        expMean    <- mxAlgebra( Mean + SexR + AgeR + ControlCovR, name="expMean")
        defs      <- list( intercept, defSex, B_Sex, meanSex, defAge, B_Age, meanAge,
                           defControlCov, B_ControlCov, meanControlCov)
      } else {
        expMean    <- mxAlgebra( Mean + SexR + AgeR, name="expMean")
        defs      <- list( intercept, defSex, B_Sex, meanSex, defAge, B_Age, meanAge)
      }
      
      zygM <-mxMatrix	(type="Full", nrow=1, ncol=1, free=FALSE, labels="data.Zygosity_1", name = "Zygosity")
      
      # Algebra for Variance/Covariance Matrices in MZ & DZ twins, and max of 2 sibs in the family
      if (!USE_TWINS_ONLY){
        # Order in family must be twins followed by non-twins
        covM <- mxAlgebra(expression= rbind( cbind(V,               Zygosity%x%A+C,  0.5%x%A+C,  0.5%x%A+C), 
                                             cbind(Zygosity%x%A+C,  V,               0.5%x%A+C,  0.5%x%A+C), 
                                             cbind(0.5%x%A+C,       0.5%x%A+C,       V,          0.5%x%A+C), 
                                             cbind(0.5%x%A+C,       0.5%x%A+C,       0.5%x%A+C,  V)), name= "expCovM")
      } else {
        # much simpler, just have the MZ vs DZ, and the relatedness is calculated with the zygosity value
        covM <- mxAlgebra(expression= rbind( cbind(V,               Zygosity%x%A+C  ), 
                                             cbind(Zygosity%x%A+C,  V               )), name= "expCovM")
      }
      
      
      # Data objects for Multiple Groups
      
      dataM <- mxData( observed=fam, type="raw" )
      
      
      # Objective objects for Multiple Groups
      
      expM <- mxExpectationNormal( covariance="expCovM", means="expMean", dimnames=selVars )
      funML <- mxFitFunctionML()
      
      
      # Combine vars
      pars <- list( pathA, pathC, pathE, covA, covC, covE, covP )
      
      modelM <- mxModel( pars, defs, expMean, zygM, covM, dataM, expM, funML, name="Main" )
      
      fitML <- mxFitFunctionMultigroup("Main.fitfunction")
      
      
      #to get confidence intervals 
      avM <- mxAlgebra(A/V, name="StandA") #calculating standard proportion for genetic component
      ci <- mxCI(c('Main.A','a11','c11','e11', 'StandA')) #confidence interval of genetic component
      AceModel <- mxModel( "ACE", pars, modelM, fitML, avM, ci)
      
      
      
      # Run ADE model
      AceFit_covSexAge <- mxRun(AceModel, intervals=T)
      AceSumm_covSexAge <- summary(AceFit_covSexAge)
      AceSumm_covSexAge
      
      # extract CIs for the genetic component
      CI<-AceSumm_covSexAge$CI["ACE.StandA[1,1]",] #now resides in CI$lbound and CI$ubound
      
      #get output
      estVA <- mxEval(a*a, AceFit_covSexAge) #genetics variance 
      estVC <- mxEval(c*c, AceFit_covSexAge) # common environmental variance, c^2
      estVE <- mxEval(e*e, AceFit_covSexAge) # unique environmental variance, e^2
      estVP <- (estVA+estVC+estVE) # total variance
      estPropVA <- estVA/estVP # standardized additive genetic variance
      estPropVC <- estVC/estVP # standardized dominance variance
      estPropVE <- estVE/estVP
      
      estACE <- rbind(cbind(estVA,estVC,estVE),cbind(estPropVA,estPropVC,estPropVE)) # table of estimates
      LL_ACE <- mxEval(objective, AceFit_covSexAge) 
      # 
      # #CE model
      
      CeModel <- mxModel( AceFit_covSexAge, name="CE" )
      CeModel <- omxSetParameters( CeModel, labels="a11", free=FALSE, values=0 )
      CeFit <- mxRun(CeModel)  
      CeSumm_covSexAge <- summary(CeFit)
      CeSumm_covSexAge
      
      estVC_mCE <- mxEval(c*c, CeFit) # additive genetic variance, a^2
      estVE_mCE <- mxEval(e*e, CeFit) # unique environmental variance, e^2
      estVP_mCE <- (estVC_mCE+estVE_mCE) # total variance
      estPropVC_mCE <- estVC_mCE/estVP_mCE # standardized additive genetic variance
      estPropVE_mCE <- estVE_mCE/estVP_mCE # standardized unique environmental variance
      estCE_mCE <- rbind(cbind(estVC_mCE,estVE_mCE),cbind(estPropVC_mCE,estPropVE_mCE)) # table of estimates
      
      LL_CE <- mxEval(objective, CeFit) # likelihood of CE model
      LRT_ACE_CE <- LL_CE - LL_ACE
      
      comp = mxCompare(AceFit_covSexAge, CeFit)
      
      
      
      # create a pretty matrix of your results
      res[row_idx,1] <- estVA
      res[row_idx,2] <- estVC
      res[row_idx,3] <- estVE
      res[row_idx,4] <- estVP
      res[row_idx,5] <- estPropVA
      res[row_idx,6] <- estPropVC
      res[row_idx,7] <- estPropVE
      res[row_idx,8] <- estVC_mCE
      res[row_idx,9] <- estVE_mCE
      res[row_idx,10] <- estVP_mCE
      res[row_idx,11] <- estPropVC_mCE
      res[row_idx,12] <- estPropVE_mCE
      res[row_idx,13] <- LL_ACE
      res[row_idx,14] <- LL_CE
      res[row_idx,15] <- LRT_ACE_CE
      res[row_idx,16] <- comp$p[2] #this is the p-value of comparing the ACE to the CE
      res[row_idx,17] <- CI$lbound
      res[row_idx,18] <- CI$ubound
      
      
      ## now do the same thing with a twins only design, if selected
      if (COMPARE_TWINS_DIRECTLY_MX) {
        print("Comparing MZ to DZ directly with openMX")
        #using the path method
        #need to grab the data again, and re-convert to ensure that we have the MZ and DZ separate
        mzData <- subset(HcpData_subset, MZ_NotMZ_NotTwin=="MZ")
        dzData <- subset(HcpData_subset, MZ_NotMZ_NotTwin=="NotMZ")
        fam_MZ<-reshape(mzData,idvar="Mother_ID",timevar="sibling_number2",direction="wide")
        fam_DZ<-reshape(dzData,idvar="Mother_ID",timevar="sibling_number2",direction="wide")
        names(fam_MZ)<-gsub("[.]","_",names(fam_MZ))
        names(fam_DZ)<-gsub("[.]","_",names(fam_DZ))
        
        #selVars <- paste(theVarName,c(rep(1,nv),rep(2,nv)),sep="_")
        
        svMe <- mean(theVar) # start value for the means
        svPa <- sqrt(var(theVar)/3) # start value for path coefficients (sqrt(variance/#ofpaths))
        
        
        # ACE Model
        # Matrices declared to store a, d, and e Path Coefficients
        pathA     <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=svPa, 
                               label="a11", name="a" ) 
        pathC     <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=svPa, 
                               label="c11", name="c" )
        pathE     <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=svPa, 
                               label="e11", name="e" )
        
        # Matrices generated to hold A, D, and E computed Variance Components
        covA      <- mxAlgebra( a %*% t(a), name="A" )
        covC      <- mxAlgebra( c %*% t(c), name="C" ) 
        covE      <- mxAlgebra( e %*% t(e), name="E" )
        
        # Matrices for covariates and linear regression coefficients
        defAge    <- mxMatrix( type="Full", nrow=1, ncol=2, free=FALSE, 
                               labels=c("data.Age_in_Yrs_1","data.Age_in_Yrs_2"), name="Age" )
        pathBAge     <- mxMatrix( type="Full", nrow=1, ncol=1, free=TRUE, values=0, 
                                  label="betaAge", name="bAge" )
        
        defSex    <- mxMatrix( type="Full", nrow=1, ncol=2, free=FALSE, 
                               labels=c("data.Gender_1","data.Gender_2"), name="Sex" )
        pathBSex     <- mxMatrix( type="Full", nrow=1, ncol=1, free=TRUE, values=0, 
                                  label="betaSex", name="bSex" )
        
        # Algebra for expected Mean Matrices in MZ & DZ twins
        meanG     <- mxMatrix( type="Full", nrow=1, ncol=2, free=TRUE, values=svMe, 
                               labels="xbmi", name="mean" )
        expMean   <- mxAlgebra(  mean + cbind(bAge%x%Age) + cbind(bSex%x%Sex), name="expMean" )
        
        # Algebra for expected Variance/Covariance Matrices in MZ & DZ twins
        covP      <- mxAlgebra(  A+C+E, name="V" )
        
        expCovMZ <- mxAlgebra(  rbind( cbind(V, A+C), 
                                       cbind(A+C, V)), name="expCovMZ" )
        
        expCovDZ <- mxAlgebra(  rbind( cbind(V, 0.5%x%A+ C), 
                                       cbind(0.5%x%A+ C , V)), name="expCovDZ" )
        
        # Data objects for Multiple Groups
        dataMZ    <- mxData( observed=fam_MZ, type="raw" )
        dataDZ    <- mxData( observed=fam_DZ, type="raw" )
        
        if (length(selVars) > 2) {
          selVars_twins_only <- selVars[1:2]
        }
        # Objective objects for Multiple Groups
        objMZ     <- mxFIMLObjective( covariance="expCovMZ", means="expMean", dimnames=selVars_twins_only )
        objDZ     <- mxFIMLObjective( covariance="expCovDZ", means="expMean", dimnames=selVars_twins_only )
        
        # Combine Groups
        pars      <- list( pathA, pathC, pathE, covA, covC, covE, covP, pathBAge, pathBSex)
        modelMZ   <- mxModel( pars, defAge, defSex, meanG, expMean, expCovMZ, dataMZ, objMZ,
                              name="MZ" )
        modelDZ   <- mxModel( pars, defAge, defSex, meanG, expMean, expCovDZ, dataDZ, objDZ,
                              name="DZ" )
        minus2ll  <- mxAlgebra( MZ.objective + DZ.objective, name="m2LL" )
        obj       <- mxAlgebraObjective( "m2LL" )
        ACEModel  <- mxModel( "ACE", pars, modelMZ, modelDZ, minus2ll, obj )
        
        # ------------------------------------------------------------------------------
        
        # Run ACE model
        ACEFit    <- mxRun(ACEModel)
        ACESumm   <- summary(ACEFit)
        
        estVA <- mxEval(a*a, ACEFit) #genetics variance 
        estVC <- mxEval(c*c, ACEFit) # common environmental variance, c^2
        estVE <- mxEval(e*e, ACEFit) # unique environmental variance, e^2
        estVP <- (estVA+estVC+estVE) # total variance
        estPropVA <- estVA/estVP # standardized additive genetic variance
        estPropVC <- estVC/estVP # standardized dominance variance
        estPropVE <- estVE/estVP
        
        estACE <- rbind(cbind(estVA,estVC,estVE),cbind(estPropVA,estPropVC,estPropVE)) # table of estimates
        LL_ACE <- mxEval(objective, AceFit_covSexAge) 
        
        res_twins_only[row_idx,1] <- estVA
        res_twins_only[row_idx,2] <- estVC
        res_twins_only[row_idx,3] <- estVE
        res_twins_only[row_idx,4] <- estVP
        res_twins_only[row_idx,5] <- estPropVA
        res_twins_only[row_idx,6] <- estPropVC
        res_twins_only[row_idx,7] <- estPropVE
        res_twins_only[row_idx,8] <- LL_ACE
      }
    }
  }
  #plot(theVar,HcpData$GM_total_volume,title(main=theVarName))
  #abline(lm(HcpData$GM_total_volume~theVar), col="red") # regression line (y~x)
  error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
    if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
      stop("vectors must be same length")
    arrows(x,upper, x, lower, angle=90, code=3, length=length, ...)
  }
  
  par(mar=c(7.1,4.1,4.1,2.1)) #add more padding at the bottom of plots (5-1 (default) --> 7.1)
  #mycols <-c("#BD395F", "#BE2598", "#A23BBF", "#4C62CB", "#007DB8", "#00888F", "#008856", "#008100", "#5B7600", "#8B6600", "#AA5200")#, "#BD395F") #old colours
  mycols <-c("#42A0B1",
              "#CE5D34",
              "#6C87CC",
              "#6EC340",
              "#B572D8",
              "#C49A2F",
              "#D64FA3",
              "#48B57C",
              "#D14E67",
              "#678732",
              "#AE71A3")
#   svg(filename=paste(img_out_dir,"GeneticComponent_allLabels.svg",sep=""), 
#       width=25, 
#       height=10, 
#       pointsize=12)
#   dev.off()
  y=res[,5] #5 is genetic, then shared env, then unique env
  if (USE_CONTROL_VAR_LM | USE_CONTROL_VAR_MX){
    barx=barplot(y,beside=T,ylim=c(0,1),main=c("Genetic Component",gsub("_"," ",gsub("^.*?_","",gsub("^.*?_","",controlVar_txt)))),col=mycols,las=3) #,col=rainbow(11)
    error.bar(barx,y,res[,"A_CI_upper"],res[,"A_CI_lower"])
    text(barx,.98,ifelse(res[,"ACE_vs_CE_p"]*33<.05,"*",""),cex=2) #stars for sig, go at 1st posn when >.05
    
    barx=barplot(y[1:11],beside=T,ylim=c(0,1),main=c("Genetic Component",gsub("_"," ",gsub("^.*?_","",gsub("^.*?_","",controlVar_txt)))),col=mycols,las=3) #,col=rainbow(11)
    error.bar(barx,y[1:11],res[1:11,"A_CI_upper"],res[1:11,"A_CI_lower"])
    text(barx,.98,ifelse(res[1:11,"ACE_vs_CE_p"]*33<.05,"*",""),cex=2)
    barx=barplot(y[12:22],beside=T,ylim=c(0,1),main=c("Genetic Component",gsub("_"," ",gsub("^.*?_","",gsub("^.*?_","",controlVar_txt)))),col=mycols,las=3) #,col=rainbow(11)
    error.bar(barx,y[12:22],res[12:22,"A_CI_upper"],res[12:22,"A_CI_lower"])
    text(barx,.98,ifelse(res[12:22,"ACE_vs_CE_p"]*33<.05,"*",""),cex=2)
    barx=barplot(y[23:33],beside=T,ylim=c(0,1),main=c("Genetic Component",gsub("_"," ",gsub("^.*?_","",gsub("^.*?_","",controlVar_txt)))),col=mycols,las=3) #,col=rainbow(11)
    error.bar(barx,y[23:33],res[23:33,"A_CI_upper"],res[23:33,"A_CI_lower"])
    text(barx,.98,ifelse(res[23:33,"ACE_vs_CE_p"]*33<.05,"*",""),cex=2)
  } else {
    barx=barplot(y,beside=T,ylim=c(0,1),main="Genetic Component",col=mycols,las=3) #,col=rainbow(11)
    barx=barplot(y[1:11],beside=T,ylim=c(0,1),main="Genetic Component",col=mycols,las=3) #,col=rainbow(11)
    error.bar(barx,y[1:11],res[1:11,"A_CI_upper"],res[1:11,"A_CI_lower"])
    text(barx,.98,ifelse(res[1:11,"ACE_vs_CE_p"]*33<.05,"*",""),cex=2)
    barx=barplot(y[12:22],beside=T,ylim=c(0,1),main="Genetic Component",col=mycols,las=3) #,col=rainbow(11)
    error.bar(barx,y[12:22],res[12:22,"A_CI_upper"],res[12:22,"A_CI_lower"])
    text(barx,.98,ifelse(res[12:22,"ACE_vs_CE_p"]*33<.05,"*",""),cex=2)
    barx=barplot(y[23:33],beside=T,ylim=c(0,1),main="Genetic Component",col=mycols,las=3) #,col=rainbow(11)
    error.bar(barx,y[23:33],res[23:33,"A_CI_upper"],res[23:33,"A_CI_lower"])
    text(barx,.98,ifelse(res[23:33,"ACE_vs_CE_p"]*33<.05,"*",""),cex=2)
  }
  #dev.off()
  
  
  
  rm(theVar)
  rm(theVarName)
  
}

## summary statistics
## create a matrix to hold the results of summaries
res_demog <- matrix(data=NA, ncol=5, nrow=4)
colnames(res_demog) <- c("N","Females", "Males","Age", "Age Range")
rownames(res_demog) <- c("MZ","DZ","Non-twin Siblings","Total")
res_demog[1,1] <- dim(subset(HcpData_subset,MZ_NotMZ_NotTwin=="MZ"))[1]
res_demog[2,1] <- dim(subset(HcpData_subset,MZ_NotMZ_NotTwin=="NotMZ"))[1]
res_demog[3,1] <- dim(subset(HcpData_subset,MZ_NotMZ_NotTwin=="NotTwin"))[1]
res_demog[4,1] <- dim(HcpData_subset)[1]
res_demog[1,2] <- dim(subset(HcpData_subset,MZ_NotMZ_NotTwin=="MZ" & Gender_Txt=="F"))[1]
res_demog[2,2] <- dim(subset(HcpData_subset,MZ_NotMZ_NotTwin=="NotMZ" & Gender_Txt=="F"))[1]
res_demog[3,2] <- dim(subset(HcpData_subset,MZ_NotMZ_NotTwin=="NotTwin" & Gender_Txt=="F"))[1]
res_demog[4,2] <- dim(subset(HcpData_subset,Gender_Txt=="F"))[1]
res_demog[1,3] <- dim(subset(HcpData_subset,MZ_NotMZ_NotTwin=="MZ" & Gender_Txt=="M"))[1]
res_demog[2,3] <- dim(subset(HcpData_subset,MZ_NotMZ_NotTwin=="NotMZ" & Gender_Txt=="M"))[1]
res_demog[3,3] <- dim(subset(HcpData_subset,MZ_NotMZ_NotTwin=="NotTwin" & Gender_Txt=="M"))[1]
res_demog[4,3] <- dim(subset(HcpData_subset,Gender_Txt=="M"))[1]
res_demog[1,4] <- mean(subset(HcpData_subset,MZ_NotMZ_NotTwin=="MZ")$Age_in_Yrs)
res_demog[2,4] <- mean(subset(HcpData_subset,MZ_NotMZ_NotTwin=="NotMZ")$Age_in_Yrs)
res_demog[3,4] <- mean(subset(HcpData_subset,MZ_NotMZ_NotTwin=="NotTwin")$Age_in_Yrs)
res_demog[4,4] <- mean((HcpData_subset)$Age_in_Yrs)
res_demog[1,5] <- paste(range(subset(HcpData_subset,MZ_NotMZ_NotTwin=="MZ")$Age_in_Yrs),collapse="-")
res_demog[2,5] <- paste(range(subset(HcpData_subset,MZ_NotMZ_NotTwin=="NotMZ")$Age_in_Yrs),collapse="-")
res_demog[3,5] <- paste(range(subset(HcpData_subset,MZ_NotMZ_NotTwin=="NotTwin")$Age_in_Yrs),collapse="-")
res_demog[4,5] <- paste(range((HcpData_subset)$Age_in_Yrs),collapse="-")


## volume statistics for our results
res_vols<-matrix(data=NA,ncol=8,nrow=length(allVars))
allVars_txt<-gsub("_volume","",allVars)
colnames(res_vols) <- c("MZ mean","DZ mean","Non-twin Sibling mean","Sample mean","MZ std","DZ std","Non-twin Sibling std","Sample std")
rownames(res_vols) <- gsub("_"," ",gsub("^.*?_","",gsub("^.*?_","",allVars_txt)))

row_idx=0
for (theVarName in allVars) {
  row_idx=row_idx+1
  if (!is.character(HcpData_subset[theVarName][,])) {
  temp <- quote(theVarName) #this and the following generate a list of names of the theVar ending in _1:4
  theVar <- HcpData_subset[,eval(temp)]
  res_vols[row_idx,1] <- mean(theVar[HcpData_subset$MZ_NotMZ_NotTwin=="MZ"])
  res_vols[row_idx,2] <- mean(theVar[HcpData_subset$MZ_NotMZ_NotTwin=="NotMZ"])
  res_vols[row_idx,3] <- mean(theVar[HcpData_subset$MZ_NotMZ_NotTwin=="NotTwin"])
  res_vols[row_idx,4] <- mean(theVar)
  res_vols[row_idx,5] <- sd(theVar[HcpData_subset$MZ_NotMZ_NotTwin=="MZ"])
  res_vols[row_idx,6] <- sd(theVar[HcpData_subset$MZ_NotMZ_NotTwin=="NotMZ"])
  res_vols[row_idx,7] <- sd(theVar[HcpData_subset$MZ_NotMZ_NotTwin=="NotTwin"])
  res_vols[row_idx,8] <- sd(theVar)
  }
}


##--------------------------------------------------------------------------------------------------
if (!SKIP_EXTRA_PLOTTING) {
  # some plotting and intercorrelations
  
  ## sub-functions from http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization
  # Get lower triangle of the correlation matrix
  get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  
  # Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }
  
  # cluster the correlation matrix
  reorder_cormat <- function(cormat){
    # Use correlation between variables as distance
    dd <- as.dist((1-cormat)/2)
    hc <- hclust(dd)
    cormat <-cormat[hc$order, hc$order]
  }

  
  if (REMOVE_FAMILY_MEMBERS_CORR_PLOTS) {
    keep_ids=c()
    #remove all but one family member from the data
    #randomly?
    for (mother in unique(HcpData_subset$Mother_ID)) {
      sibs=HcpData_subset$ID[HcpData_subset$Mother_ID==mother]
      keep_ids <- c(keep_ids,sample(sibs,1))
    }
    corr_subset <- subset(HcpData_subset, ID %in% keep_ids)
  } else {
    corr_subset <- HcpData_subset
  }
  
  ## correlation matrix of our input numerical values
  is_numeric_row <- lapply(corr_subset,is.character)
  is_numeric_row <- !unlist(is_numeric_row,use.names = FALSE) #true false of the rows that are numeric
  corr_subset <- corr_subset[,is_numeric_row]
  
  corr_subset <- corr_subset[,colnames(corr_subset) %in% allVars] #only keep those that were in our input vars
  removeVars=c("Height","Weight","Handedness")
  corr_subset <- corr_subset[!(colnames(corr_subset) %in% removeVars)]

  colnames(corr_subset) <- gsub("_"," ",gsub("^.*?_","",gsub("^.*?_","",colnames(corr_subset))))
  colnames(corr_subset) <- gsub("volume","",colnames(corr_subset))
  
  library("Hmisc")
  library("ggplot2")
  library("reshape2")
  corr_res<-rcorr(as.matrix(corr_subset)) #1=r, #2=n, #3=p
  #ggplot(data = corr_subset) 
  #melted_cormat <- melt(corr_res$r)
  #ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + geom_tile()
  
  ## standard heatmap
  corr_res_p <- corr_res$r
  p_val <- 0.05/(dim(corr_res_p)[1]*dim(corr_res_p)[1]-dim(corr_res_p)[1])/2
  corr_res_p[!corr_res$P<p_val] <- NA
  cormat <- corr_res_p
  upper_tri <- get_upper_tri(cormat)
  diag(upper_tri) <- NA #remove the diagonal from the plot
  melted_cormat <- melt(upper_tri, na.rm = FALSE)
  
  # Heatmap
#  minn=min(melted_cormat$value)
#  maxx=max(melted_cormat$value)
  ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name="Pearson\nCorrelation") +
    theme_minimal()+ 
    theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                     size = 10, hjust = 1),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background = element_rect(fill = "lightgrey"))+
    coord_fixed()
#   ## Reorder the correlation matrix
  
#   cormat <- reorder_cormat(corr_res$r)
#   reor_upper_tri <- get_upper_tri(cormat)
#   # Melt the correlation matrix
#   melted_cormat <- melt(reor_upper_tri, na.rm = TRUE)
#   
#   # Create a ggheatmap
#   ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
#     geom_tile(color = "white")+
#     scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
#                          midpoint = 0, limit = c(-1,1), space = "Lab", 
#                          name="Pearson\nCorrelation") +
#     theme_minimal()+ # minimal theme
#     theme(axis.text.x = element_text(angle = 45, vjust = 1, 
#                                      size = 12, hjust = 1),
#           axis.title.x=element_blank(),
#           axis.title.y=element_blank())+
#     coord_fixed()
#   
#   # Print the heatmap
#   print(ggheatmap)
  
  ## some additional plots of MZ vs DZ
  # help here: http://www.sthda.com/english/wiki/ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page-r-software-and-data-visualization
  # library("cowplot")
  # ## correlation matrix of our input numerical values
  # MZ <- subset(fam,MZ_NotMZ_NotTwin_1=="MZ") #twins always in the first two sets of cols
  # DZ <- subset(fam,MZ_NotMZ_NotTwin_1=="NotMZ")
  # row_idx=0
  # for (theVarName in allVars) {
  #   row_idx=row_idx+1
  #   if (!is.character(HcpData_subset[theVarName][,])) {
  #     print(theVarName)
  #     #then we can plot it!
  #     xvar=paste(theVarName,"1",sep="_")
  #     yvar=paste(theVarName,"2",sep="_")
  #     p1 <- ggplot(MZ,aes(x=MZ[xvar][,],y=MZ[yvar][,]))+
  #       geom_point(colour="red")+
  #       geom_smooth(method=lm,colour="red")+
  #       ggtitle("MZ")+
  #       xlab("")+
  #       ylab("Volume (mm^3)")
  #     p2 <-ggplot(DZ,aes(x=DZ[xvar][,],y=DZ[yvar][,]))+
  #       geom_point(colour="blue")+
  #       geom_smooth(method=lm,colour="blue")+
  #       ggtitle("DZ")+
  #       xlab("")+
  #       ylab("")
  #     p3<-plot_grid(p1,p2,labels=theVarName,label_size = 14, hjust = -2.5, vjust = 1.5)
  #     print(p3)
  #     }
  # }
  # 
}