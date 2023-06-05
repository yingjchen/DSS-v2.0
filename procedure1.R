###DSS 2.0, sDSS/zDSS/rDSS computations
###set up
##Install required libraries
packages.required <- packages.required <- c("matrixStats","dplyr","reshape","reshape2", "scales", "drc", "caTools", "ggplot2", "data.table", 
                                            "stringr","MESS", "BiocManager","svMisc", "egg", "pheatmap")
packages.bio <- c("sva", "pcaMethods")
packages.new <- packages.required[!(packages.required %in% installed.packages()[,"Package"])]
if(length(packages.new)) install.packages(packages.new)
if (!requireNamespace(packages.bio, quietly = TRUE))
    BiocManager::install(packages.bio)


##load the packages
lapply(packages.required, library, character.only = T)
lapply(packages.bio, library, character.only = T)

##replace the working directory replace '/path/to/working/directory/' with the desired path
path_to_working_directory <- '/path/to/working/directory/'
setwd(dir = path_to_working_directory)

download.file(url = 'https://github.com/yingjchen/DSS-v2.0/archive/refs/heads/main.zip', destfile = 'DSS-v2.0-main.zip')
unzip('DSS-v2.0-main.zip') 
setwd(dir = file.path(path_to_working_directory, 'DSS-v2.0-main'))

source('./DSS.R')
source('./HelperFunctions.R')


###1. calculate the DSS scores (DSS1, DSS2, DSS3) from the cell viability data
##load the example data: the cell viability of ~500 drugs were tested in 3 FIMM patient samples
path_to_exampledata <- './exampleData_procedure1.csv'
df_dose.responses <- read.csv(path_to_exampledata, header = T,sep = ',', check.names = F)

#set the viability argument to ‘FALSE’ if input data is cell growth inhibition data
df_dose.responses.list <- DOSE_RESPONSE_PROCESS(df_dose.responses, viability = T)

##compute BREEZE drug sensitivity metrics, i.e. DSS1, DSS2, DSS3, Breeze AUC and relative IC50
#set the graph argument to 'TRUE' to generate the fitted dose-response curves (inhibition vs dose) under the directory ~/IC50/; the graph argument initialized with a default value: FALSE
df.metrics <- CALC_METRICS(df_dose.responses.list[[1]], df_dose.responses.list[[2]], graph = FALSE)

##here we select DSS2 as patient DSS
patients.dss <- as.data.frame(acast(df.metrics,df.metrics$Patient.num ~ df.metrics$drug , value.var  = 'DSS2'))
#plot a heatmap showing DSS2 of the drugs with highest standard deviations across the samples;
#The argument proportion can be adjusted to specify the proportion of most variable drugs (e.g. 10% drugs selected with the code provided below);
#the graph argument initialized with a default value: 1
HEATMAP_SD(patients.dss, proportion = 0.1)

##save DSS scores and other metrics
#write.table(df.metrics, file = './Results_exampleData_procedure1.csv',sep = ',',append=F,row.names = F,col.names = T)


###2. load DSS data of 10 FIMM controls (FO5A plate)
path_to_controldss <- './controls/File_1_Drugname_response_DSS_10Healthy.txt'
controls.dss <- read.csv(path_to_controldss, header = T, sep = '\t', row.names = 1,stringsAsFactors = F, check.names = F)

##compute mean, SD, median and MAD of control DSSs
controls.summary <- as.data.frame(rbind(colMeans(as.matrix(controls.dss)),colSds(as.matrix(controls.dss)),colMedians(as.matrix(controls.dss)), colMads(as.matrix(controls.dss))))
rownames(controls.summary ) <- c('mean', 'sd', 'median', 'mad')



###3.search the drugs of interest from the drug library (527 drugs) 
##3.1 map the drugs with drug synonyms
path_to_druglibrary <- './controls/File_2_Drugname_library_527D.txt'
df_drug.library <- read.csv(path_to_druglibrary, header = T, sep = '\t', row.names = NULL,stringsAsFactors = F, check.names = F)
patients.dss <- DRUG_FILTER_SYNONYMS(patients.dss, df_drug.library)

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


###5.1. Plotting the drug response distributions, combined plot (optional)
sample_id <- 'AML_013_01'
sample_dss <- SAMPLE_DSS_CONCAT(patients.dss, patients.sdss, patients.zdss, patients.rdss, sample_id = sample_id)
CHEMO_TAREGTED_PLOT(sample_dss, metric = 'DSS')
###5.2. Plotting the drug response distributions, density plot(optional)
SELECTIVE_SCORE_PLOT(sample_dss)
