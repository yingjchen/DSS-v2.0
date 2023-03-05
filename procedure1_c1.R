###DSS 2.0, sDSS/zDSS/rDSS computations
###set up
##Install required libraries
packages.required <- c("plotly", "scales", "parallel", "foreach", "gridExtra", "grid", "graphics", "gplots",
                       "ggplot2", "raster", "xtable","Rcpp","dplyr", "drc", "caTools", "gsubfn", "gtools", "data.table", "doSNOW","stringr",
                       "MESS", "reshape", "reshape2", "matrixStats")
packages.new <- packages.required[!(packages.required %in% installed.packages()[,"Package"])]
if(length(packages.new)) install.packages(packages.new)

##load the packages
lapply(packages.required, library, character.only = !0)


rm(list = ls())
setwd('D:/AML/DSS/')

#url_to_dsscalc <- 'https://github.com/xxx' 
#source(url_to_dsscalc)
source('./exampledata/github repo/DSS.R')

###1. calculate the DSS scores (DSS1, DSS2, DSS3) from the cell viability data
##load the example data: the cell viability of 10 drugs were tested in 2 FIMM patient samples
path_to_exampledata <- './exampledata/github repo/exampledata_procedure1_c1.csv'

df_dose.responses <- read.csv(path_to_exampledata, header = T,sep = ",",check.names = F)

colnames(df_dose.responses) <- c('drug', 'Concentration','Patient.num' ,'response.raw')
df_dose.responses$ID <- paste0(df_dose.responses$drug,"_Patient ", df_dose.responses$Patient.num)
# compute cell inhibition
df_dose.responses$inhibition_percent <- 100 - df_dose.responses$response.raw*100  

df_dose.responses.grouped <- df_dose.responses %>% 
  group_by(ID, Patient.num, drug) %>% 
  summarise(conc = max(Concentration)) %>%
  as.data.frame()  

##compute BREEZE drug sensitivity metrics
iter <- nrow(df_dose.responses.grouped)
list.DSS_metrics <- list()
list.believe <- list()
list.IC50 <- list()

for(it in 1:iter){
  df.DSS.ctx <- CALC_IC50_EC50_DSS(it, df_dose.responses.grouped, df_dose.responses, DSS_typ="2", readoutCTX = T)
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
df_scores <- cbind( df_dose.responses.grouped, df.DSS_metrics, df.IC50_metrics, df.believe)

##here we select DSS2 as patient DSS
patients.dss <- as.data.frame(acast(df_scores,df_scores$Patient.num ~ df_scores$drug , value.var  = 'DSS2'))

##save DSS scores and other metrics
#write.table(df_scores, file = './Results_exampleData_procedure1.csv',sep = ',',append=F,row.names = F,col.names = T)


###2. load DSS data of 10 FIMM controls (FO5A plate)
path_to_controldss <- './exampledata/github repo/File_1_Drugname_response_DSS_10Healthy.txt'
controls.dss <- read.csv(path_to_controldss, header = T, sep = '\t', row.names = 1,stringsAsFactors = F, check.names = F)

##compute mean, SD, median and MAD of control DSSs
controls.summary <- as.data.frame(rbind(colMeans(as.matrix(controls.dss)),colSds(as.matrix(controls.dss)),colMedians(as.matrix(controls.dss)), colMads(as.matrix(controls.dss))))
rownames(controls.summary ) <- c('mean', 'sd', 'median', 'mad')

path_to_controldss <- './exampledata/github repo/File_1_Drugname_response_DSS_10Healthy.txt'
controls.dss <- read.csv(path_to_controldss, header = T, sep = '\t', row.names = 1,stringsAsFactors = F, check.names = F)

controls.summary <- as.data.frame(rbind(colMeans(as.matrix(controls.dss)),colSds(as.matrix(controls.dss)),colMedians(as.matrix(controls.dss)), colMads(as.matrix(controls.dss))))
rownames(controls.summary ) <- c('mean', 'sd', 'median', 'mad')


###3.search the drugs of interest from the drug library (527 drugs) 
##3.1 map the drugs with drug synonyms
path_to_druglibrary <- './exampledata/github repo/File_2_Drugname_library_527D.txt'
df_drug.library <- read.csv(path_to_druglibrary, header = T, sep = '\t', row.names = NULL,stringsAsFactors = F, check.names = F)

colnames_drug <- toupper(colnames(patients.dss))
patients.dss.map <- patients.dss

for (i in df_drug.library$drug){
  drugsets = as.character(df_drug.library[df_drug.library$drug == i, "synonyms"])
  drugsynonyms = unique(unlist(strsplit(drugsets,";")))
  patients.dss.map['drugname', colnames_drug %in% drugsynonyms] = as.character(i)
}

##filter out unmatched drugs
patients.dss <- patients.dss[, !(is.na(patients.dss.map['drugname',]))]
colnames(patients.dss) <- patients.dss.map['drugname',]

##the drugs can also may be mapped with InChIKey from the drug library


###4. compute patient sDSS, zDSS and rDSS
##4.1 compute patient sDSS
patients.sdss <- patients.dss - slice(controls.summary['mean', colnames(patients.dss)],rep(1:n(), each = nrow(patients.dss)))

##4.2 compute patient zDSS
patients.zdss <- (patients.dss - slice(controls.summary['mean', colnames(patients.dss)],rep(1:n(), each = nrow(patients.dss))))/(slice(controls.summary['sd', colnames(patients.dss)],rep(1:n(), each = nrow(patients.dss))) + 1)

##4.3 compute patient rDSS
patients.rdss <- (patients.dss - slice(controls.summary['median', colnames(patients.dss)],rep(1:n(), each = nrow(patients.dss))))/(slice(controls.summary['mad', colnames(patients.dss)],rep(1:n(), each = nrow(patients.dss))) + 1)


###save selective DSS scores
write.table(patients.dss, file = './Results_exampledata_DSS_procedure1.txt',sep = '\t',append=F,row.names = T,col.names = T)
write.table(patients.sdss, file = './Results_exampledata_sDSS_procedure1.txt',sep = '\t',append=F,row.names = T,col.names = T)
write.table(patients.zdss, file = './Results_exampledata_zDSS_procedure1.txt',sep = '\t',append=F,row.names = T,col.names = T)
write.table(patients.rdss, file = './Results_exampledata_rDSS_procedure1.txt',sep = '\t',append=F,row.names = T,col.names = T)
