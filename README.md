# DSS-v2.0
DSS-v2.0 is a new pipeline to better quantify selective drug responses between cancer patients and healthy controls.


# Instructions
R version 3.5.1 or newer is required.

# Under Linux/Unix
```
# download the example data and R scripts
git clone https://github.com/yingjchen/DSS-v2.0.git

# change the directory
cd ./DSS-v2.0

# start the R program
R
```
# Start by loading libraries and functions
```r
lapply(c("matrixStats","dplyr","reshape","reshape2", "scales", "drc", "caTools","ggplot2", "data.table", "stringr","MESS", "BiocManager","svMisc", "egg", "pheatmap", "sva", "pcaMethods"), library, character.only = T)

source('./DSS.R')
source('./HelperFunctions.R')
```

# Selective DSS calculation example
## Load the example data and compute DSS with *CALC_METRICS* function updated from BREEZE (<https://github.com/potdarswapnil/Breeze>)

```r
# load the ex vivo dose-response profiles (cell viability at five drug concentrations)
df_dose.responses <- read.csv('./exampleData_procedure1.csv', header = T,sep = ",",check.names = F)

head(df_dose.responses)
```

    ##    DRUG_NAME CONCENTRATION SCREEN_NAME    CELL_VIABILITY
    ## 1 Nelarabine         10000  AML_013_01            0.3012
    ## 2 Nelarabine          1000  AML_013_01            0.5565
    ## 3 Nelarabine           100  AML_013_01            0.7578
    ## 4 Nelarabine            10  AML_013_01            0.9133
    ## 5 Nelarabine             1  AML_013_01            0.8610
    ## 6 Decitabine         10000  AML_013_01            0.4829


```r
# calculate the percentage of growth inhibition 
df_dose.responses.list <- DOSE_RESPONSE_PROCESS(df_dose.responses, viability = T)

# calculate DSS metrics (DSS1, DSS2, DSS3), AUC and relative IC50
df.metrics <- CALC_METRICS(df_dose.responses.list[[1]], df_dose.responses.list[[2]])
```

## Import the control sample DSS profiles
```r
controls.dss <- read.csv('./controls/File_1_Drugname_response_DSS_10Healthy.txt', header = T, sep = '\t', row.names = 1,stringsAsFactors = F, check.names = F)
```
## Calculate the selective drug response scores(sDSS, zDSS, and rDSS)
```r
# compute the descriptive statistics of DSSs for each drug over 10 controls
controls.summary <- as.data.frame(rbind(colMeans(as.matrix(controls.dss)),colSds(as.matrix(controls.dss)),colMedians(as.matrix(controls.dss)), colMads(as.matrix(controls.dss))))

# define names of statistics
rownames(controls.summary ) <- c('mean', 'sd', 'median', 'mad')

# let's set DSS2 as the drug response metrics of patient samples (as an example)
patients.dss <- as.data.frame(acast(df.metrics,df.metrics$Patient.num ~ df.metrics$drug , value.var  = 'DSS2'))
patients.dss[, 1:3]
```
    ##                    1-methyl-D-tryptophan 4-hydroxytamoxifen 8-amino-adenosine
    ## AML_003_01                     0                9.2              25.4
    ## AML_004_01                     0                2.3              22.6
    ## AML_013_01                     0                4.3              23.9
```r
# normalize and scale patient-specific responses to drugs with control DSS profiles
patients.sdss <- patients.dss - slice(controls.summary['mean', colnames(patients.dss)],rep(1:n(), each = nrow(patients.dss)))

patients.zdss <- (patients.dss - slice(controls.summary['mean', colnames(patients.dss)],rep(1:n(), each = nrow(patients.dss))))/(slice(controls.summary['sd', colnames(patients.dss)],rep(1:n(), each = nrow(patients.dss))) + 1)

patients.rdss <- (patients.dss - slice(controls.summary['median', colnames(patients.dss)],rep(1:n(), each = nrow(patients.dss))))/(slice(controls.summary['mad', colnames(patients.dss)],rep(1:n(), each = nrow(patients.dss))) + 1)
```
