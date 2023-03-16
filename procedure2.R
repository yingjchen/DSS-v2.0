###batch effect correction and heatmap visualization

##define the required packages 
packages.required <- packages.required <- c("matrixStats","dplyr", "ggplot2", "data.table", "BiocManager", "pheatmap")
packages.bio <- c("sva", "pcaMethods")


##########Install required libraries, if needed##########
packages.new <- packages.required[!(packages.required %in% installed.packages()[,"Package"])]
if(length(packages.new)) install.packages(packages.new)
if (!requireNamespace(packages.bio, quietly = TRUE))
    BiocManager::install(packages.bio)
##########Install required libraries, if needed##########


##load the packages
lapply(packages.required, library, character.only = T)
lapply(packages.bio, library, character.only = T)


##########Download the github repo and set the working directory, if needed##########
##the working directory needs to be changed
path_to_working_directory <- '/path/to/working/directory/'
setwd(dir = path_to_working_directory)

download.file(url = 'https://github.com/yingjchen/DSS-v2.0/archive/refs/heads/main.zip', destfile = 'DSS-v2.0-main.zip')
unzip('DSS-v2.0-main.zip') 
setwd(dir = file.path(path_to_working_directory, 'DSS-v2.0-main'))
##########Download the github repo and set the working directory, if needed##########


path_to_exampledata <- './exampleData_procedure2.csv'
df.dss <- read.csv(path_to_exampledata, header = T,sep = ",",  row.names = 1, check.names = F)

##make PPCA plot of DSS
df.dss.1 <- df.dss[, 1 : (ncol(df.dss) - 3)]
res_ppca <- pca(data.matrix(df.dss.1), method = 'ppca', nPcs = 2, seed = 1)
score_ppca <-  as.data.frame(scores(res_ppca))
score_ppca$group <- paste(df.dss$cohort,  df.dss$status,  sep = ' ')
ggplot(score_ppca, aes(x = PC1, y = PC2, color = group)) +
  geom_point() + labs(title = "DSS",  x = "PC1", y = "PC2") +
  theme_classic()
ggsave("./example_DSS_ppca.pdf", height = 10, width = 10)

##perform ComBat correction  
df.dss.combat <- ComBat(dat = t(df.dss.1), batch = as.factor(df.dss$cohort), mod = NULL, par.prior = F, prior.plots = F)

##make PPCA plot of ComBat DSS
res_ppca_combat <- pca(data.matrix(t(df.dss.combat)), method = 'ppca', nPcs = 2, seed = 1)
score_ppca_combat <-  as.data.frame(scores(res_ppca_combat))
score_ppca_combat$group <- paste(df.dss$cohort,  df.dss$status,  sep = ' ')
ggplot(score_ppca_combat, aes(x = PC1, y = PC2, color = group)) +
  geom_point() + labs(title = "ComBat DSS",  x = "PC1", y = "PC2") +
  theme_classic()
ggsave("./example_ComBatDSS_ppca.pdf", height = 10, width = 10)

###Heatmap visualization
df.dss.2 <- as.data.frame(df.dss[df.dss$status != 'controls', ])
df.dss.2$datasource <- paste(df.dss.2$cohort, df.dss.2$plate,  sep = ' ')
r_ <- data.frame(datasource = as.factor(df.dss.2$datasource))
rownames(r_) <- rownames(df.dss.2)

##make heatmap of DSS
p1 <- pheatmap(df.dss.2[, 1:(ncol(df.dss.2) - 4)], annotation_row = r_, show_colnames = F, show_rownames = F,cluster_rows = T, clustering_distance_cols = "minkowski", color= colorRampPalette(c("lightgrey","steelblue"))(100 ))
ggsave("./example_DSS_heatmap.pdf", p1, height = 10,width = 10)

##make heatmap of ComBat DSS
df.dss.3 <- as.data.frame(t(df.dss.combat)[df.dss$status != 'controls', ])
p2 <- pheatmap(df.dss.3, annotation_row = r_, show_colnames = F, show_rownames = F,cluster_rows = T, clustering_distance_cols = "minkowski",color = colorRampPalette(c("lightgrey","steelblue"))(100 ))
ggsave("./example_ComBat_DSS_heatmap.pdf", p2, height = 10,width = 10)

##calculate ComBat rDSS
controls.dss <- t(df.dss.combat)[df.dss$status == 'controls', ]
controls.summary <- as.data.frame(rbind(colMeans(as.matrix(controls.dss)),colSds(as.matrix(controls.dss)),colMedians(as.matrix(controls.dss)), colMads(as.matrix(controls.dss))))
rownames(controls.summary ) <- c('mean', 'sd', 'median', 'mad')
patients.rdss <- (df.dss.3 - slice(controls.summary['median', colnames(df.dss.3)],rep(1:n(), each = nrow(df.dss.3 ))))/(slice(controls.summary['mad', colnames(df.dss.3 )],rep(1:n(), each = nrow(df.dss.3 ))) + 1)

##make heatmap of ComBat rDSS
p3 <- pheatmap(t(patients.rdss), annotation_col = r_, show_colnames = F, show_rownames = F,cluster_cols = T, clustering_distance_cols = "minkowski",color = colorRampPalette(c("lightgrey","steelblue"))(100 ))
ggsave("./example_ComBat_rDSS_heatmap.pdf", p3, height = 10,width = 10)
