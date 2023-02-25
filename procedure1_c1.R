###DSS 2.0, sDSS/zDSS/rDSS computations
###set up
#Install required libraries
packages.required <- c("plotly", "scales", "parallel", "foreach", "gridExtra", "grid", "graphics", "gplots",
         "ggplot2", "raster", "xtable","Rcpp","dplyr", "drc", "caTools", "gsubfn", "gtools", "data.table", "doSNOW","stringr",
         "MESS", "reshape", "reshape2", "matrixStats")
packages.new <- packages.required[!(packages.required %in% installed.packages()[,"Package"])]
if(length(packages.new)) install.packages(packages.new)


rm(list = ls())
setwd('D:/AML/DSS/')
source('DSS.R')

###1. calculate the DSS scores (DSS1, DSS2, DSS3) from the cell viability data
#load the example data: the cell viability of 10 drugs were tested in 2 FIMM patient samples
df_dose.responses = read.csv('./exampleData_procedure1.csv',header = T,sep = ",",check.names = F)

colnames(df_dose.responses) = c('drug', 'concentration','Patient.num' ,'response.raw')
df_dose.responses$ID <- paste0(df_dose.responses$drug,"_Patient ", df_dose.responses$Patient.num)
# compute cell inhibition
df_dose.responses$inhibition_percent <- 100 - df_dose.responses$response.raw*100  
# compute raw AUC
df.AUC.raw <- df_dose.responses %>%
  group_by(ID, Patient.num, drug) %>%
  summarise(AUC.raw = Log10.AUC(concentration, response.raw)) %>%
  as.data.frame()

# Compute BREEZE drug sensitivity metrics
xpr_tbl <- df_dose.responses
xpr_tbl$Concentration <- xpr_tbl$concentration
drug_wells_ <- df.AUC.raw
iter <- nrow(drug_wells_)
list.DSS_metrics <- list()
list.believe <- list()
list.IC50 <- list()
for(it in 1:iter){
  
  df.DSS.ctx <- CALC_IC50_EC50_DSS(it, drug_wells_, xpr_tbl, DSS_typ="2", readoutCTX = T, path)
  
  list.believe[[it]] <- data.frame(believe_DSS=df.DSS.ctx[[3]])
  list.DSS_metrics[[it]] <- data.frame(AUC=df.DSS.ctx[[1]]$AUC,DSS1=df.DSS.ctx[[5]],
                                       DSS2=df.DSS.ctx[[6]],DSS3=df.DSS.ctx[[7]])
  list.IC50[[it]] <- data.frame(IC50=df.DSS.ctx[[8]])
  
  rm( df.DSS.ctx)
  
  svMisc::progress(it/iter*100)
  
}

df.believe <- bind_rows(list.believe)
df.DSS_metrics <- bind_rows(list.DSS_metrics)
df.IC50_metrics <- bind_rows(list.IC50)

df_scores <- cbind(df.AUC.raw, df.DSS_metrics, df.IC50_metrics, df.believe)

###save DSS scores and other metrics
#write.table(df_scores, file = './Results_exampleData_procedure1.csv',sep = ',',append=F,row.names = F,col.names = T)


###2. load DSS data of 10 FIMM controls (FO5A plate)
controls_dss <- read.csv('./controls/FIMM_DSS_10_controls_FIMMID.txt',header = T, sep = '\t', row.names = 1,stringsAsFactors = F, check.names = F)
#compute mean, SD, median and MAD of control DSSs
controls_summary = as.data.frame(rbind(colMeans(as.matrix(controls_dss)),colSds(as.matrix(controls_dss)),colMedians(as.matrix(controls_dss)), colMads(as.matrix(controls_dss))))
rownames(controls_summary ) = c('mean', 'sd', 'median', 'mad')

### 3. compute patient sDSS, zDSS and rDSS
#here we select DSS2 as patient DSS
patients_dss <- as.data.frame(acast(df_scores,df_scores$Patient.num ~ df_scores$drug , value.var  = 'DSS2', ))

#3.1 compute patient sDSS
patients_sdss <- patients_dss - slice(controls_summary['mean', colnames(patients_dss)],rep(1:n(), each = nrow(patients_dss)))
#3.2 compute patient zDSS
patients_zdss <- (patients_dss - slice(controls_summary['mean', colnames(patients_dss)],rep(1:n(), each = nrow(patients_dss))))/(slice(controls_summary['sd', colnames(patients_dss)],rep(1:n(), each = nrow(patients_dss))) + 1)
#3.3 compute patient rDSS
patients_rdss <- (patients_dss - slice(controls_summary['median', colnames(patients_dss)],rep(1:n(), each = nrow(patients_dss))))/(slice(controls_summary['mad', colnames(patients_dss)],rep(1:n(), each = nrow(patients_dss))) + 1)


### save selective DSS scores
write.table(patients_dss, file = './Results_exampledata_DSS_procedure1.txt',sep = '\t',append=F,row.names = T,col.names = T)
write.table(patients_sdss, file = './Results_exampledata_sDSS_procedure1.txt',sep = '\t',append=F,row.names = T,col.names = T)
write.table(patients_zdss, file = './Results_exampledata_zDSS_procedure1.txt',sep = '\t',append=F,row.names = T,col.names = T)
write.table(patients_rdss, file = './Results_exampledata_rDSS_procedure1.txt',sep = '\t',append=F,row.names = T,col.names = T)
