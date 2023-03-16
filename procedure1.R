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
df_dose.responses <- read.csv(path_to_exampledata, header = T,sep = ",",check.names = F)

df_dose.responses.list <- DOSE_RESPONSE_PROCESS(df_dose.responses, viability = T)

##compute BREEZE drug sensitivity metrics, i.e. DSS1, DSS2, DSS3, Breeze AUC and relative IC50
df.metrics <- CALC_METRICS(df_dose.responses.list[[1]], df_dose.responses.list[[2]])

##here we select DSS2 as patient DSS
patients.dss <- as.data.frame(acast(df.metrics,df.metrics$Patient.num ~ df.metrics$drug , value.var  = 'DSS2'))

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


###5. Plotting the drug response distributions (optional)
sample_id <- 'AML_013_01'
sample_dss <- SAMPLE_DSS_CONCAT(sample_id = sample_id)
p1 <- ggplot(sample_dss) +  
  geom_histogram(aes(x = DSS, y = after_stat(density), fill= drugclass), binwidth = 1, color = "black", alpha = .6, position="identity", show.legend = T) +
  labs(title = sample_id,  x="DSS", y = "Density") + theme_classic()
p2 <- ggplot(sample_dss, aes(y = drugclass, x = DSS, fill= drugclass)) +
  geom_boxplot(aes(x = drugclass, y = DSS, fill = drugclass), outlier.shape = NA, alpha = .6,  colour = "black", show.legend = F) +
  geom_jitter(aes(color = drugclass), size = 3, alpha = .8,  show.legend = F)+
  theme_classic() + coord_flip()
p3 <- ggarrange(p1, p2, nrow = 2)
ggsave("./example_DSS_distribution.pdf", p3, height = 10,width = 10)
