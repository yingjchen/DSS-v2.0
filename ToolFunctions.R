

DOSE_RESPONSE_PROCESS <- function(dose_responses, viability = T){
  colnames(dose_responses) <- c('drug', 'Concentration','Patient.num' ,'response.raw')
  dose_responses$ID <- paste0(dose_responses$drug,"_Patient ", dose_responses$Patient.num)
  if (viability){
    dose_responses$inhibition_percent <- 100 - dose_responses$response.raw*100
    
  }else{
    dose_responses$inhibition_percent <- dose_responses$response.raw*100  
  }
  
  return(dose_responses)
}




CALC_METRICS <- function(dose_responses, dose_responses_grouped){
  iter <- nrow(dose_responses_grouped)
  list.AUC_DSS <- list()
  list.believe <- list()
  list.IC50 <- list()
  for(it in 1:iter){
    df.DSS.ctx <- CALC_IC50_EC50_DSS(it, dose_responses_grouped, dose_responses, DSS_typ = '2')
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
  return(df.metric)
}



DRUG_FILTER_SYNONYMS <- function(inputdrug, druglibrary){
  drug_names <- toupper(colnames(inputdrug))
  inputdrug.map <- inputdrug
  
  for (i in druglibrary$drug){
    drugsets = as.character(druglibrary[druglibrary$drug == i, "synonyms"])
    drugsynonyms = unique(unlist(strsplit(drugsets,";")))
    inputdrug.map['drugname', drug_names %in% drugsynonyms] = as.character(i)
  }
  
  inputdrug <- inputdrug[, !(is.na(inputdrug.map['drugname',]))]
  colnames(inputdrug) <- inputdrug.map['drugname',]
  return(inputdrug)
}




SAMPLE_DSS_CONCAT <- function(sample_id){
  sample_dss<- data.frame(DSS = as.numeric(patients.dss[sample_id, ]), sDSS = as.numeric(patients.sdss[sample_id, ]), zDSS = as.numeric(patients.zdss[sample_id,]), row.names = colnames(patients.dss) )
  
  
  sample_dss$drugclass <- unlist(lapply(rownames(sample_dss), function(d_){ 
    as.character(df_drug.library[df_drug.library$drug == as.character(d_), 'drugclass'])
  }))
  return(sample_dss)
}
