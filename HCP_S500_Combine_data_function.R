rm(list = ls()) #clean the workspace

library("Hmisc")
library("ggplot2")
library("reshape2")

#HCP S500 data files
restr='/home/chris/Documents/Dropbox/Projects/Working/CoBrA_Lab/TractREC/source/bx/HCP_S500_RESTRICTED_steelec_10_30_2015_11_14_9.csv'
unrestr='/home/chris/Documents/Dropbox/Projects/Working/CoBrA_Lab/TractREC/source/bx/HCP_S500_Bx_unrestricted_steelec_9_17_2015_9_46_45.csv'
unrestr_FS_vols='/home/chris/Documents/Dropbox/Projects/Working/CoBrA_Lab/TractREC/source/bx/HCP_S500_unrestricted_hcp_freesurfer.csv'

#bx_dir='/home/chris/Documents/Dropbox/Projects/Working/CoBrA_Lab/TractREC/processing/bx/'
#all_structures_S500_volume_restr_unrestr = "2016_03_HCP_allSibs_all_labels_volume_FA_corrected_0p35_restr_unrestr.csv"
#all_volumes_fname='2016_03_HCP_allLabels_FA_0p35_volume.csv'

bx_dir='/home/chris/Documents/Dropbox/Projects/Working/CoBrA_Lab/TractREC/processing/bx/S500/'
img_out_dir='/home/chris/Documents/Dropbox/Projects/Working/CoBrA_Lab/TractREC/processing/bx/S500/images/'



## chose the metric file(s?) to work with
all_volumes_all_metrics_S500_FA_fname='2016_03_HCP_FA_allLabels_FA_0p35_all.csv'
#all_volumes_all_metrics_S500_MD_fname='2016_03_HCP_MD_allLabels_FA_0p35_all.csv'
#all_volumes_all_metrics_S500_KM_fname='2016_03_HCP_KM_allLabels_FA_0p35_all.csv'
#all_volumes_all_metrics_S500_t1wdivt2w_fname='2016_03_HCP_t1wdivt2w_allLabels_FA_0p35_all.csv'


## --- SET SOME VARIABLES
REMOVE_FAMILY_MEMBERS_CORR_PLOTS =     TRUE # creates correlations and plots when only a single individual from each family is selected (randomly)
RANDOM_FAMILY_REMOVAL =                FALSE
USE_CONTROL_VAR_LM =                   FALSE #regress out age and sex
USE_BONF_CORRECTION =                  FALSE
SAVE_PLOTS =                           FALSE
USE_RIGHT_HANDERS_ONLY =               FALSE #handedness >50

#11 colours generated from http://tools.medialab.sciences-po.fr/iwanthue/
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

metric_file = all_volumes_all_metrics_S500_FA_fname #choose your metric file
img_fname_head="FA" # set the name of the metric (no spaces)
value_types = c('volume')#,'median')#,'median')#,'std') #valid metrics in this file to loop through
## --- 

df_restr<-read.csv(restr,header=T)
df_unrestr<-read.csv(unrestr,header=T)
df_unrestr_FS_vols<-read.csv(unrestr_FS_vols,header=T)
common_cols <- names(df_unrestr)[names(df_unrestr) %in% names(df_unrestr_FS_vols)] #a, when a is in b
#df_metric=read.csv(paste(bx_dir,all_volumes_fname,sep=""),header=T)
df_metric<-read.csv(paste(bx_dir,metric_file,sep=""),header=T)

## merge the restricted and unrestricted datasets, then merge with the extracted metric data
df_HCPbx<-merge(x=df_restr,y=df_unrestr,by="Subject")
df_HCPbx<-merge(x=df_HCPbx,y=df_unrestr_FS_vols,by=common_cols) #merge all common columns, since some data is repeated
data<-merge(x=df_metric,y=df_HCPbx,by.x="ID",by.y="Subject")

##  ---------------------------------------------------------------------------------------------------
## format the data so that we can work with it in OpenMX
## update Zygosity with version corrected for families that included a single MZ twin

#call this little function to do so
HCP.clean_data_df <- function(data_df, filter_twins=FALSE) {
  #recodes and cleans data for analysis
  
  data_df$MZ_NotMZ_NotTwin <- data_df$Zygosity #record what the individual was in a new var
  data_df$Recoded_Zygosity[data_df$Zygosity=="MZ"] <- 1
  data_df$Recoded_Zygosity[data_df$Zygosity=="NotMZ"] <- .5
  data_df$Recoded_Zygosity[data_df$Zygosity=="NotTwin"] <- .5
  data_df$Zygosity<-data_df$Recoded_Zygosity
  
  ## reorder the siblings such that the 1st and 2nd are always the twins
  ## XXX this is because the covariance matrix REQUIRES twins first
  ## TWINS (MZ or DZ) must always be in the 1st two columns of the family, so ordered 1 and 2 here
  for (mother in unique(data_df$Mother_ID)) {
    ## use "-" to reverse the zygosity data so that rank is decreasing
    ## add twin_stat*2 to push the twins higher even if their zygosity is .5 (DZ twins), so they will come out first
    sibs<-unique(data_df$ID[data_df$Mother_ID==mother])
    twin_stat<-data_df$Twin_Stat[data_df$Mother_ID==mother]=="Twin" #true or false that each of these is a twin
    ordered_zygosity=rank(-(data_df$Zygosity[data_df$ID %in% sibs] + twin_stat*rep(2,length(twin_stat))),ties.method="first") #rank by their zygosities, decreasing
    data_df$sibling_number2[data_df$Mother_ID==mother]<-ordered_zygosity
    
    #calculate the sibling count and Twins_Count for potential filtering later
    twin_sibs <- subset(data_df,Mother_ID==mother & Twin_Stat == "Twin")$ID
    data_df$Sibling_Count[data_df$Mother_ID==mother] <- length(sibs)
    data_df$Twins_Count[data_df$Mother_ID==mother] <- length(twin_sibs)
  }
  
  ## remove the individuals who are twins, but their twin sib is not in the dataset
  # only use this if you are interested in it!
  if (filter_twins) {
  data_df <- subset(data_df,Twins_Count != 1) #singleton twins removed
  data_df <- subset(data_df,Twins_Count != 0) #families with no twins removed
  data_df <- subset(data_df,Sibling_Count > 1) #also need to have more than one sibling?
  }
  #remove problematic individuals
  data_df <- subset(data_df,ID != 105014)
  
  
  ## make gender numerical
  data_df$Gender_Txt <- data_df$Gender
  data_df$Gender <- NA
  data_df$Gender[data_df$Gender_Txt=="M"] <- 1
  data_df$Gender[data_df$Gender_Txt=="F"] <- 2
  return (data_df)
}
HCP.get_unrelated_data <- function(data_df,individual_label="ID",family_grouping_label="Mother_ID",random=FALSE) {
  #remove all but one family member from the data
  #returns the first sibling by default
  all_fams <- unique(data_df[family_grouping_label])[,]
  keep_ids=c()
  for (fam in all_fams) {
    sibs<- data_df[ data_df[[family_grouping_label]] == fam , ][individual_label][,]
    if (random) {
      keep_ids <- c(keep_ids,sample(sibs,1))
    }
    else {
      keep_ids <- c(keep_ids,sibs[1]) # XXX JUST TAKING THE FIRST SIB SO THIS DOESN"T CHANGE  
    }
  }
  data_df <- subset(data_df,ID %in% keep_ids)
  return (data_df)
}
HCP.control_for_variables <- function(data_df,control_var_name_list,label_subset=NULL,label_skip=NULL) {
  ## regress out the effect of control_var_name_list on all variables or subset of variables in dataframe
  if (is.null(label_subset)) { 
    var_list <- names(data_df)
  }
  else {
    var_list <- label_subset
  }
  var_list_class <-sapply(data_df,class) #determine the class of the data
  for (theVarName in var_list) {
    if (is.numeric(data_df[,theVarName])) { #if the data is a number (int, float, etc), we use it
      if (!theVarName %in% label_skip) {
        print(theVarName)
        lm_fit <- lm(data_df[theVarName][,] ~ ., data_df[control_var_name_list])
        data_df[theVarName] <- resid(lm_fit)
      }
    }
  }
  return (data_df)
}
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
HCP.plotcorrmat <- function(cormat,title=NULL,value=NULL) {
  cormat.melt <- melt(cormat, na.rm=FALSE)
  p <- ggplot(data = cormat.melt , aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name="Pearson\nCorrelation") +
    theme_minimal()+ 
    theme(axis.text.x = element_text(angle = 90, vjust = .5, 
                                     size = 10, hjust = 1),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background = element_rect(fill = "lightgrey"))+
    coord_fixed()+
    ggtitle(paste(img_fname_head,value,sep=" "))
  return(p)
}
# cluster the correlation matrix
reorder_cormat <- function(cormat,method="ward.D2"){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd,method=method)
  cormat <-cormat[hc$order, hc$order]
  return(cormat)
}

data <- HCP.clean_data_df(data)
if (REMOVE_FAMILY_MEMBERS_CORR_PLOTS) {
  data <- HCP.get_unrelated_data(data,random = RANDOM_FAMILY_REMOVAL)
}
if (USE_RIGHT_HANDERS_ONLY) {
  data <- subset(data,Handedness>50) #here we can further select variables etc
}


#data$Gender <- as.numeric(data$Gender)
##  ---------------------------------------------------------------------------------------------------
# some plotting and intercorrelations
for (value in value_types){
  allVars <- colnames(data[grep("label_",colnames(data))]) #select all column names that have "label_" in them
  allVars <- allVars[allVars != "label_file"]
  allVars <- allVars[grep(value,allVars)]
  
  CBVars <-allVars[1:33] #use this to calculate the sum of all CB GM volume
  allVars<-CBVars #also only care about them right now
  
  figure_text_append_control <- "" #text to add to the figure output name
  

  if (USE_CONTROL_VAR_LM) {

#    allVars<-c(CBVars,'label_000_CB_GM_volume')
#    controlVar_txt<-'label_000_CB_GM_volume'
#    controlVars <- c("Age_in_Yrs","Gender",controlVar_txt)
    controlVars <- c("Age_in_Yrs","Gender")
    data_subset <- HCP.control_for_variables(data_subset,controlVars,label_subset=allVars)
    
  } else {
    data_subset <- data
  }
  
  corr_subset<-data_subset #don't need to do this anymore, but too lazy to change names
  
  ## correlation matrix of our input numerical values
  is_numeric_row <- lapply(corr_subset,is.character)
  is_numeric_row <- !unlist(is_numeric_row,use.names = FALSE) #true false of the rows that are numeric
  corr_subset <- corr_subset[,is_numeric_row]
  
  corr_subset <- corr_subset[,colnames(corr_subset) %in% allVars] #only keep those that were in our input vars
  removeVars=c("Height","Weight","Handedness")
  corr_subset <- corr_subset[!(colnames(corr_subset) %in% removeVars)]
  
  colnames(corr_subset) <- gsub("_"," ",gsub("^.*?_","",gsub("^.*?_","",colnames(corr_subset))))
  colnames(corr_subset) <- gsub(value,"",colnames(corr_subset))
  
  corr_res<-rcorr(as.matrix(corr_subset)) #1=r, #2=n, #3=p
  ## standard heatmap
  corr_res_p <- corr_res$r #we might correct these results by p-value
  
  if (USE_BONF_CORRECTION) {
    p_val <- 0.05/(dim(corr_res_p)[1]*dim(corr_res_p)[1]-dim(corr_res_p)[1])/2
    corr_res_p[!corr_res$P<p_val] <- NA
  }
  
  cormat <- corr_res_p
  upper_tri <- get_upper_tri(cormat)
  diag(upper_tri) <- NA #remove the diagonal from the plot
  melted_cormat <- melt(upper_tri, na.rm = FALSE)
  
  # Heatmap
  #  minn=min(melted_cormat$value)
  #  maxx=max(melted_cormat$value)
  p<-ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
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
    coord_fixed()+
    ggtitle(paste(img_fname_head,value,sep=" "))
  print(p)
  fname=paste(Sys.Date(),img_fname_head,"CB","corrmat",value,sep="_")
  print(fname)
  if (SAVE_PLOTS) { 
    ggsave(paste(img_out_dir,fname,".svg",sep=""),device=svg)
  }
  ## create and save violin plots for each hemi/vermis
  vars_LH <- allVars[1:11]
  vars_RH <- allVars[12:22]
  vars_vermis <- allVars[23:33]
  vars_grouped <- list(vars_LH,vars_vermis,vars_RH)
  
  if (!USE_CONTROL_VAR_LM) {
    loop_idx<-0
    for (var_group in vars_grouped) {
      loop_idx<- loop_idx+1
      vars_subset <- c("ID","Gender_Txt",var_group)
      plot_subset <- data_subset[,vars_subset]
      colnames(plot_subset) <- gsub(paste(" ",value,sep=""),"",gsub("_"," ",gsub("^.*?_","",gsub("^.*?_","",colnames(plot_subset)))))
      colnames(plot_subset)[colnames(plot_subset)=="Txt"]<-"Gender_Txt" #put this back, since we cut it out :-/
      plot_subset.m <- reshape2::melt(plot_subset, id.vars = c("ID","Gender_Txt")) #now one entry per row
      colnames(plot_subset.m) <- c('ID',"Gender","Lobule","Value")
      if (loop_idx == 2) { #we are vermal, different scale
        ylim=c(0,2500)
      } else {
        ylim=c(0,17500)
      }
      p <- ggplot(data=plot_subset.m,aes(x = Lobule, y = Value, fill=Gender))+ 
        geom_boxplot(notch=TRUE)+
        theme_minimal()+
        theme(axis.text.x=element_text(size=10,angle=90,vjust=.5))+
        xlab("")+
        ylab(value)+
        ylim(ylim)+
        #geom_jitter(position=position_jitter(width=.1, height=0))+
        #scale_y_log10(limits=c(1,20000))+
        ggtitle(img_fname_head)
      #ylim(c(0,5))
      #ggtitle(value)
      #  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10, hjust = 1))
      print(p)
      fname=paste(Sys.Date(),img_fname_head,"CB","boxplot",value,loop_idx,sep="_")
      print(fname)
      if (SAVE_PLOTS){ 
        ggsave(paste(img_out_dir,fname,".svg",sep=""),plot = p,device=svg)
      }
    }
    
  }
  
}

#descriptive stats for the dataset
data_subset$L_CB_Volume <- rowSums(data_subset[,vars_grouped[[1]]])
data_subset$Vermal_CB_Volume <- rowSums(data_subset[,vars_grouped[[2]]])
data_subset$R_CB_Volume <- rowSums(data_subset[,vars_grouped[[3]]])

descVars<-c(allVars,'L_CB_Volume','Vermal_CB_Volume','R_CB_Volume')

res_desc<-matrix(data=NA,ncol=6,nrow=length(descVars))
colnames(res_desc)<-c('mean','sd','median','min','max','n')
rownames(res_desc)<-descVars
res_desc[,1] <- colMeans(data_subset[,descVars])
res_desc[,2] <- apply(data_subset[,descVars], 2, sd)
res_desc[,3] <- apply(data_subset[,descVars], 2, median)
res_desc[,4:5] <- t(apply(data_subset[,descVars], 2, range))
res_desc[,6] <-apply(data_subset[,descVars], 2, length)


res_ttest_R_L<-matrix(data=NA,ncol = 8,nrow=length(vars_LH)+1)
colnames(res_ttest_R_L)<-c('LH mean','RH mean','mean of the diff','t','p','df','ci','description')
rownames(res_ttest_R_L)<-c(vars_LH,'L_CB_Volume')
idx <- 0
for (theVar in vars_LH) {
  idx <- idx + 1
  res_t <- t.test(data_subset[,vars_LH[idx]],data_subset[,vars_RH[idx]],paired=TRUE, alternative="two.sided")
  res_ttest_R_L[idx,1] <- mean(data_subset[,vars_LH[idx]])
  res_ttest_R_L[idx,2] <- mean(data_subset[,vars_RH[idx]])
  res_ttest_R_L[idx,3] <- res_t$estimate
  res_ttest_R_L[idx,4] <- res_t$statistic
  res_ttest_R_L[idx,5] <- res_t$p.value
  res_ttest_R_L[idx,6] <- res_t$parameter
  res_ttest_R_L[idx,7] <- paste(res_t$conf.int,collapse = " ")
  res_ttest_R_L[idx,8] <- paste(paste(vars_LH[idx],vars_LH[idx],sep=" vs. "),res_t$method,res_t$alternative,sep="; ")
  }
idx <- idx + 1
res_t <- t.test(data_subset[,'L_CB_Volume'],data_subset[,'R_CB_Volume'],paired=TRUE, alternative="two.sided")
res_ttest_R_L[idx,1] <- mean(data_subset[,'L_CB_Volume'])
res_ttest_R_L[idx,2] <- mean(data_subset[,'R_CB_Volume'])
res_ttest_R_L[idx,3] <- res_t$estimate
res_ttest_R_L[idx,4] <- res_t$statistic
res_ttest_R_L[idx,5] <- res_t$p.value
res_ttest_R_L[idx,6] <- res_t$parameter
res_ttest_R_L[idx,7] <- paste(res_t$conf.int,collapse = " ")
res_ttest_R_L[idx,8] <- paste(paste('L_CB_Volume','R_CB_Volume',sep=" vs. "),res_t$method,res_t$alternative,sep="; ")
rm(idx,res_t)

idx <- 0
for (theVar in vars_LH) {
  idx <- idx + 1
  var1 <- vars_LH[idx]
  var2 <- vars_RH[idx]
  
  }
# source('/home/chris/Documents/Dropbox/Projects/Working/CoBrA_Lab/TractREC/processing/bx/S500/func_CreateRadialPlot.R')
# #radial plot?
# # define the metrics
# metrics=c("FA","MD","KM","t1wdivt2w")
# value_types = c('volume')#,'volume','median')#,'std') #valid metrics in this file to loop through
# 
# df=read.csv(paste(bx_dir,all_volumes_all_metrics_S500_FA_fname,sep=""),header=T)
# colnames(df) <- paste("FA",colnames(df),sep="_")
# df_metrics=df
# 
# ## select the variables of interest
# allVars <- colnames(df[grep("label_",colnames(df))]) #select all column names that have "label_" in them
# allVars <- allVars[allVars != "FA_label_file"] 
# #allVars <- allVars[grep(value,allVars)]
# 
# 
# CBVars <-allVars[1:33] #use this to calculate the sum of all CB GM volume
# allVars<-CBVars #also only care about them right now
# 
# df=read.csv(paste(bx_dir,all_volumes_all_metrics_S500_MD_fname,sep=""),header=T)
# colnames(df) <- paste("MD",colnames(df),sep="_")
# df_metrics=merge(x=df_metrics,y=df,by.x="FA_ID",by.y="MD_ID")
# 
# df=read.csv(paste(bx_dir,all_volumes_all_metrics_S500_KM_fname,sep=""),header=T)
# colnames(df) <- paste("KM",colnames(df),sep="_")
# df_metrics=merge(x=df_metrics,y=df,by.x="FA_ID",by.y="KM_ID")
# 
# df=read.csv(paste(bx_dir,all_volumes_all_metrics_S500_t1wdivt2w_fname,sep=""),header=T)
# colnames(df) <- paste("t1wdivt2w",colnames(df),sep="_")
# df_metrics=merge(x=df_metrics,y=df,by.x="FA_ID",by.y="t1wdivt2w_ID")
# 
# #colnames(df_metrics[colnames(df_metrics)=="FA_ID"]) <- "ID"
# 
# ## merge the restricted and unrestricted datasets with the metric data
# df_all=merge(x=df_metrics,y=df_HCPbx,by.x="FA_ID",by.y="Subject")
# 
# # ## radar plots
# # for (theVarName in allVars) {
# #  theVar_metrics <- paste(metrics,theVarName,sep="_") 
# #  for (value in value_types) {
# #    subset_vars <- theVar_metrics[grep(value,theVar_metrics)]
# #    vars_means <- colMeans(df_all[subset_vars])
# #  #  CreateRadialPlot(data.frame(vars_means))
# #    vars_means.m=reshape2:::melt.matrix(vars_means)
# #    colnames(vars_means.m)[1]
# #    ggplot(vars_means.m,x=Var1,y=value)+
# #      geom_path(alpha=0.5,colour="black")+
# #      geom_point(alpha=0.2, colour="blue")
# #  }
# # }
# 


# ## -- correlation between cortical volumes and our CB volumes
FULL_MATRIX_REORDERED = FALSE #use the full matrix for bivariate correlation

FSVars<-names(data)[grep("FS",names(data))]
FSVars_GrayVol<-FSVars[grep("GrayVol",FSVars)]

selected_FSVars <-FSVars_GrayVol[grep("FS",FSVars_GrayVol)]
selected_CBVars <-CBVars[1:33] #11 per hemi/vermis

selected_IQVars <- names(data[grep("PMAT",names(data))])
selected_IQVars <- selected_IQVars[2:length(selected_IQVars)]

## alternative ordering of CBVars
#selected_CBVars <- c(mapply(c,CBVars[1:11],CBVars[12:22],c(CBVars[24:33],CBVars[23])),recursive=TRUE) #one after the other
#selected_CBVars <- c(mapply(c,CBVars[1:11],CBVars[12:22]),CBVars[23:33],recursive=TRUE) # L,R repeated, then vermal

## make some useful values (just volume for now)
data$label_000_CB_GM_volume <- rowSums(data[,selected_CBVars])
data$label_000_CB_L_volume <- rowSums(data[,selected_CBVars[1:10]])
data$label_000_CB_R_volume <- rowSums(data[,selected_CBVars[12:21]])
data$label_000_CB_Vermis_volume <- rowSums(data[,selected_CBVars[23:32]])
data$label_000_CB_Flocculus_volume <- rowSums(data[,selected_CBVars[c(11,22,33)]])

superset_CBVars <- c("label_000_CB_GM_volume",
                     "label_000_CB_L_volume",
                     "label_000_CB_R_volume",
                     "label_000_CB_Vermis_volume",
                     "label_000_CB_Flocculus_volume")


corr_subset <- HCP.control_for_variables(data,c("Age_in_Yrs","Gender"),label_subset=c(selected_CBVars,selected_FSVars,selected_IQVars,superset_CBVars),label_skip=NULL)


selected_data_CB <-corr_subset[,selected_CBVars] #select the CB volumes
selected_data_FS <-corr_subset[,selected_FSVars] #select the FS GM volumes

corr_subset <- corr_subset[,c(selected_CBVars,selected_FSVars)] # the full one, in case we asked for it



names(selected_data_FS) <- gsub("FS ","",
                                gsub("_"," ",
                                     gsub("GrayVol","",
                                          names(selected_data_FS))))
names(selected_data_CB) <- gsub("_"," ",
                                gsub("^.*?_","",
                                     gsub("^.*?_","",
                                          gsub("volume","",
                                               names(selected_data_CB)))))

names(corr_subset) <-gsub("\\d+ ","",
                          gsub("_"," ",
                               gsub("FS_","",
                                    gsub("label_","-----------",
                                         gsub("_volume","",
                                              gsub("GrayVol","",
                                                   names(corr_subset)))))))

## -- START PLOTTING OF CB vs FS CTX volumes -- ##
#corr_FS <- rcorr(as.matrix(corr_subset)) #1=r, #2=n, #3=p
corr_FS_v2 <- cor(selected_data_CB,selected_data_FS)
bonf_p <- .05/(length(selected_CBVars)*length(selected_FSVars))
#cormat <- corr_FS$r
cormat <-corr_FS_v2

if (FULL_MATRIX_REORDERED) {
  corr_FS <- rcorr(as.matrix(corr_subset)) #1=r, #2=n, #3=p
  cormat <- reorder_cormat(corr_FS$r)
  #cormat[cormat<0] <-NA
}

#p_vals <- corr_FS$P
#cormat[cormat>.2]<-NA
#melted_cormat <- melt(cormat, na.rm = FALSE)
p <- HCP.plotcorrmat(cormat,title=img_fname_head,value=value) #melts it for us, weeeee!
print(p)

## -- END PLOTTING OF CB vs CTX -- ##

## -- second order correlations of CB-CTX corr vs CB-CTX corr -- ##
# we take the cormat as input, and correlate it again
if (!FULL_MATRIX_REORDERED) {
  rotate <- function(x) t(apply(x, 2, rev))
  #2nd order connectivity matrix, correlate the vector of correlations...
  cormat_2 <- rcorr(t(as.matrix(cormat))) #rotate so that we are CB to CB, not CTX to CTX
  #melted_cormat_2 <- melt(cormat_2$r, na.rm = FALSE)
  p2 <- HCP.plotcorrmat(cormat_2$r,title=img_fname_head,value=value) #melts it for us, weeeee!
  print(p2)
  #plot(agnes(as.dist(1-cormat_2$r)))
  
}
## -- END PLOTTING OF FULL CB vs CB -- ##

## -- ward cluster the similarities of structural corr matrix --##
dist_mat <- as.dist(1-cormat_2$r)

CB_conn_corr <- hclust(dist_mat,method="ward.D2")
plot(CB_conn_corr)
CB_conn_corr.mat <-(reorder_cormat(cormat_2$r))
CB_conn_corr.mat.melt <- melt(CB_conn_corr.mat , na.rm = FALSE)

p3 <- HCP.plotcorrmat(CB_conn_corr.mat,title=img_fname_head,value=value) #melts it for us, weeeee!
print(p3)

# library(pvclust)
# cluster.bootstrap <- pvclust(cormat, 
#                              nboot=1000, 
#                              method.dist="correlation",
#                              parallel=TRUE,
#                              method.hclust = "ward.D2")
# plot(cluster.bootstrap)
# pvrect(cluster.bootstrap)

# #http://research.stowers-institute.org/efg/R/Visualization/cor-cluster/
# distmat <- 1-cormat
# dist <-as.dist(distmat)
# plot(hclust(dist), 
#      main="Dissimilarity = 1 - Correlation", xlab="")
# library(cluster) 
# plot(agnes(dist,method="ward"))
# 

#functional to calculate the partial correlation matrix
#just to convince myself that I also did this correctly with my regression and then correlation of residuals.
# XXX yes, they are the same to the 10-15 decimal point.
# returns r and p, in list (break out with res[[1]] and res[[2]])
# library(ppcor)
# HCP.calc_pcormat <- function(df,vars_to_correlate_x=NULL,vars_to_correlate_y=NULL,control_vars=c("Age_in_Yrs","Gender")) {
#   #calculate partial correlation and output a matrix
#   pcormat_r <- matrix(data=NA, ncol=length(vars_to_correlate_x), 
#                       nrow=length(vars_to_correlate_y))
#   pcormat_p <- matrix(data=NA, ncol=length(vars_to_correlate_x), 
#                       nrow=length(vars_to_correlate_y))
#   rownames(pcormat_r) <- (vars_to_correlate_y)
#   colnames(pcormat_r) <- (vars_to_correlate_x)
#   rownames(pcormat_p) <- (vars_to_correlate_y)
#   colnames(pcormat_p) <- (vars_to_correlate_x)
#   for (mycol in 1:length(vars_to_correlate_x)) {
#     for (myrow in 1:length(vars_to_correlate_y)) {
#       if (vars_to_correlate_x[myrow] == vars_to_correlate_y[mycol]) {
#         pcormat_r[myrow,mycol]<-1
#         pcormat_p[myrow,mycol]<-0  
#       } else {
#         myres <- pcor.test(df[,vars_to_correlate_y[myrow]],df[,vars_to_correlate_x[mycol]],df[,control_vars])
#         pcormat_r[myrow,mycol]<-myres$estimate
#         pcormat_p[myrow,mycol]<-myres$p.value
#       }
#     }
#   }
#   return(list(pcormat_r,pcormat_p))
# }
# 
# 
# b<-HCP.calc_pcormat(data,vars_to_correlate_y=selected_CBVars,vars_to_correlate_x=selected_CBVars, control_vars=c("Age_in_Yrs","Gender"))
# ps<-b[[2]]
# rs<-b[[1]]
# rs[ps>0.05/33^2]<-NA
# p4<-HCP.plotcorrmat(rs)
# print(p4)
# 
