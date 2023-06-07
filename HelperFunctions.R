
DOSE_RESPONSE_PROCESS <- function(dose_responses, viability = T){
  colnames(dose_responses) <- c('drug', 'Concentration','Patient.num' ,'response.raw')
  dose_responses$ID <- paste0(dose_responses$drug, "_", dose_responses$Patient.num)
  if (viability){
    dose_responses$inhibition_percent <- 100 - dose_responses$response.raw*100
    
  }else{
    dose_responses$inhibition_percent <- dose_responses$response.raw*100  
  }
  dose_responses_grouped <- dose_responses %>% 
    group_by(ID, Patient.num, drug) %>% 
    summarise(conc = max(Concentration)) %>%
    as.data.frame()
  
  return(list(dose_responses, dose_responses_grouped))
}



CALC_METRICS <- function(dose_responses, dose_responses_grouped, graph = F){
  start.time <- Sys.time()
  iter <- nrow(dose_responses_grouped)
  list.AUC_DSS <- list()
  list.believe <- list()
  list.IC50 <- list()
  unlink("./IC50", recursive = T) ;dir.create("./IC50")
  
  for(it in 1:iter){
    df.DSS.ctx <- CALC_IC50_EC50_DSS(it, dose_responses_grouped, dose_responses, DSS_typ = '2', graph = graph)
    list.believe[[it]] <- data.frame(believe_DSS = df.DSS.ctx[[3]])
    list.AUC_DSS[[it]] <- data.frame(AUC = df.DSS.ctx[[1]]$AUC,DSS1 = df.DSS.ctx[[5]],
                                     DSS2 = df.DSS.ctx[[6]],DSS3 = df.DSS.ctx[[7]])
    list.IC50[[it]] <- data.frame(IC50 = df.DSS.ctx[[8]])
    rm( df.DSS.ctx)
    svMisc::progress(it/iter*100)
  }
  df.believe <- bind_rows(list.believe)
  list.AUC_DSS <- bind_rows(list.AUC_DSS)
  df.IC50_metrics <- bind_rows(list.IC50)
  df.metric <- cbind( dose_responses_grouped, list.AUC_DSS, df.IC50_metrics, df.believe)
  message("Finished DSS computations in ", round(Sys.time() - start.time, 2), units(Sys.time() - start.time))
  return(df.metric)
  
}



DRUG_FILTER_SYNONYMS <- function(inputdrug, druglibrary){
  start.time <- Sys.time()
  drug_names <- toupper(colnames(inputdrug))
  inputdrug.map <- inputdrug
  
  for (i in druglibrary$drug){
    drugsets = as.character(druglibrary[druglibrary$drug == i, "synonyms"])
    drugsynonyms = unique(unlist(strsplit(drugsets,";")))
    inputdrug.map['drugname', drug_names %in% drugsynonyms] = as.character(i)
  }
  
  inputdrug <- inputdrug[, !(is.na(inputdrug.map['drugname',]))]
  colnames(inputdrug) <- inputdrug.map['drugname',]
  message("Finished drug mapping in ", round(Sys.time() - start.time, 2), units(Sys.time() - start.time),  paste(', ' ,as.character(ncol(patients.dss)), 'drugs in our compound library') )
  return(inputdrug)
}



SAMPLE_DSS_CONCAT <- function(df.dss, df.sdss, df.zdss, df.rdss, sample_id){
  sample_dss<- data.frame(DSS = as.numeric(df.dss[sample_id, ]), sDSS = as.numeric(df.sdss[sample_id, ]), zDSS = as.numeric(df.zdss[sample_id,]), rDSS = as.numeric(df.rdss[sample_id,]), row.names = colnames(df.dss) )

  sample_dss$drugclass <- unlist(lapply(rownames(sample_dss), function(d_){ 
    as.character(df_drug.library[df_drug.library$drug == as.character(d_), 'drugclass'])
  }))
  return(sample_dss)
}



HEATMAP_SD <- function(df, proportion = 1, filename = ""){
  if (proportion <= 0 | proportion > 1 | !is.numeric(proportion)){
    stop("The value of the argument proportion should be larger than 0 and no larger than 1")
  }
  
  if (nrow(df) > 1 & ceiling(ncol(df)*proportion) > 1){
    # we need at least two rows/columns for clustering
    drug_variable <- order(-colSds(as.matrix(df)))[1 : ceiling(ncol(df)*proportion)]
    p0 <- pheatmap(df[, drug_variable],  show_colnames = T, show_rownames = T,  clustering_distance_cols = "minkowski")
    if (filename == "") {ggsave(paste0("./step8_example_Breeze_DSS_", ceiling(ncol(df)*proportion), "drugs_heatmap.pdf"), p0, height = 10,width = 10)}else
    {ggsave(filename, p0, height = 10,width = 10)}
    message("Finished heatmap for DSS of ", ceiling(ncol(df)*proportion), " most variable drugs across ", nrow(df)," samples")
    
  }else if (nrow(df) != 1 & nrow(df) != 1 & ceiling(ncol(df)*proportion) == 1){
    p0 <- pheatmap(df,  show_colnames = T, show_rownames = T,  clustering_distance_cols = "minkowski")
    if (filename == "") {ggsave(paste0("./step8_example_Breeze_DSS_", ncol(df), "drugs_heatmap.pdf"), p0, height = 10,width = 10)}else
    {ggsave(filename, p0, height = 10,width = 10)}
    message("Finished heatmap for DSS of ", ncol(df), " drugs across ", nrow(df)," samples")
    
  }else {
    #bar plot of all the compounds or samples, if there is one compound or one sample
    if (nrow(df) == 1) df <- as.data.frame(t(df))
    df_order <-  data.frame(DSS = df[ order(as.matrix(df), decreasing=F), ], ID = row.names(df)[order(as.matrix(df), decreasing=F)])
    p0 <- ggplot(df_order) + geom_bar(aes(y = DSS, x = ID, fill = 'red'), stat = "identity", show.legend = F)  +
      labs(title = "",  x = "", y = "DSS") + scale_x_discrete(limits = df_order$ID) + theme_classic()
    
    if (filename == "") {ggsave("./step8_example_Breeze_DSS_barplot.pdf", p0, height = 10, width = 10)}else
    {ggsave(filename, p0, height = 10, width = 10)}
    message("Finished bar plot for DSS of ", nrow(df), " drugs across ", ncol(df)," samples")
  }
}


CHEMO_TAREGTED_PLOT <- function(df, metric, filename = ""){
  if (metric == 'DSS'){
    p1 <- ggplot(df) +  
      geom_histogram(aes(x = DSS, y = after_stat(density), fill= drugclass), binwidth = 1, color = "black", alpha = .6, position="identity", show.legend = T) +
      labs(title = sample_id,  x="DSS", y = "Density") + theme_classic()
    p2 <- ggplot(sample_dss, aes(x = drugclass, y = DSS, fill= drugclass)) +
      geom_boxplot(aes(x = drugclass, y = DSS, fill = drugclass), outlier.shape = NA, alpha = .6,  colour = "black", show.legend = F) +
      geom_jitter(aes(color = drugclass), size = 3, alpha = .8,  show.legend = F)+
      theme_classic() + coord_flip()
    p3 <- ggarrange(p1, p2, nrow = 2)
    
  }else if (metric == 'sDSS'){
    p1 <- ggplot(df) +  
      geom_histogram(aes(x = sDSS, y = after_stat(density), fill= drugclass), binwidth = 1, color = "black", alpha = .6, position="identity", show.legend = T) +
      labs(title = sample_id,  x="sDSS", y = "Density") + theme_classic()
    p2 <- ggplot(sample_dss, aes(x = drugclass, y = sDSS, fill= drugclass)) +
      geom_boxplot(aes(x = drugclass, y = sDSS, fill = drugclass), outlier.shape = NA, alpha = .6,  colour = "black", show.legend = F) +
      geom_jitter(aes(color = drugclass), size = 3, alpha = .8,  show.legend = F)+
      theme_classic() + coord_flip()
    p3 <- ggarrange(p1, p2, nrow = 2)

  }else if (metric == 'zDSS'){
    p1 <- ggplot(df) +  
      geom_histogram(aes(x = zDSS, y = after_stat(density), fill= drugclass), binwidth = 1, color = "black", alpha = .6, position="identity", show.legend = T) +
      labs(title = sample_id,  x="zDSS", y = "Density") + theme_classic()
    p2 <- ggplot(sample_dss, aes(x = drugclass, y = zDSS, fill= drugclass)) +
      geom_boxplot(aes(x = drugclass, y = zDSS, fill = drugclass), outlier.shape = NA, alpha = .6,  colour = "black", show.legend = F) +
      geom_jitter(aes(color = drugclass), size = 3, alpha = .8,  show.legend = F)+
      theme_classic() + coord_flip()
    p3 <- ggarrange(p1, p2, nrow = 2)
    
  }else if (metric == 'rDSS'){
    p1 <- ggplot(df) +  
      geom_histogram(aes(x = rDSS, y = after_stat(density), fill= drugclass), binwidth = 1, color = "black", alpha = .6, position="identity", show.legend = T) +
      labs(title = sample_id,  x="rDSS", y = "Density") + theme_classic()
    p2 <- ggplot(sample_dss, aes(x = drugclass, y = rDSS, fill= drugclass)) +
      geom_boxplot(aes(x = drugclass, y = rDSS, fill = drugclass), outlier.shape = NA, alpha = .6,  colour = "black", show.legend = F) +
      geom_jitter(aes(color = drugclass), size = 3, alpha = .8,  show.legend = F)+
      theme_classic() + coord_flip()
    p3 <- ggarrange(p1, p2, nrow = 2)
    
  }else{
    stop("The argument metric should be one of 'DSS', 'sDSS', 'zDSS', or 'rDSS'")
  }
  if (filename == "") {ggsave(paste0("./step13_chemo_taregted_example_", metric, "_distribution.pdf"), p3, height = 10, width = 10)}else
  {ggsave(filename, p3, height = 10, width = 10)}
  message("Finished data distribution plots of ", metric, " in chemo and targeted drugs")
}

SELECTIVE_SCORE_PLOT <- function(df, filename = ""){
  df_new <- data.frame(score = c(df$DSS, df$sDSS, df$zDSS, df$rDSS ), class = factor(rep(c('DSS', 'sDSS','zDSS', 'rDSS' ), c(nrow(df), nrow(df), nrow(df), nrow(df))), levels = c('DSS', 'sDSS','zDSS', 'rDSS' )))
  p1 <- ggplot(df_new) + 
    geom_density(aes(x = score, fill = class), alpha = 0.3, lwd = 1, bw = 2,  show.legend = T) +
    labs(title = sample_id,  x="", y = "Density") + theme_classic()
  if (filename == "") {ggsave("./step13_density_example_scores_distribution.pdf", p1, height = 10, width = 10)}else
  {ggsave(filename, p1, height = 10, width = 10)}
  message("Finished data distribution plots of DSS, sDSS, zDSS and rDSS in one sample")
}
  
PCA_FUNC <- function(df){
  if (sum(is.na(df)) > 0){
    res_pca <- pca(data.matrix(df), method = 'ppca', nPcs = 2, seed = 1,scale = NULL, center = TRUE)
    score_pca <-  as.data.frame(scores(res_pca))
    message(round(sum(is.na(df))/ (0.01 *nrow(df) * ncol(df)), 2), "% missing values exist, using PPCA analysis")
  } else{
    res_pca <- prcomp(data.matrix(df), scale = FALSE, center = TRUE)
    score_pca <-  as.data.frame(res_pca$x[, c('PC1', 'PC2')])
    message("No missing values exist, using PCA analysis")
  }
  colnames(score_pca) <- c('PC1', 'PC2')
  return(score_pca)
}

