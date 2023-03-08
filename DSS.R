
# Libraries 
lapply(c("matrixStats","dplyr","reshape","reshape2","plotly", "scales", "parallel", "foreach", "gridExtra", "grid", "graphics", "gplots",
          "xtable","Rcpp"), library, character.only = !0)
lapply(c("drc", "caTools", "ggplot2", "gsubfn", "gtools", "data.table", "stringr",
         "MESS"), library, character.only = !0)
# DSS1/DSS2/DSS3/AUC/DSS/EC50 computations  

#Raw AUC functions
Log10.AUC <- function(conc, resp){
  resp = resp[order(conc)]
  conc = conc[order(conc)]
  conc = log10(conc)
  
  A = 0
  for(j in 2:length(conc)){
    a <- (resp[j]) *(conc[j]-conc[j-1])/(max(conc)-min(conc))
    A = A + a
    rm(a,j)
  }
  A
}

#DSS/EC50 functions adapted from BREEZE (https://github.com/potdarswapnil/Breeze)
dss<-function(ic50,slope,max,min.conc.tested,max.conc.tested,y=10,DSS.type=2,concn_scale=1e-9){
  #rdata should be in in format containing IC50, SLOPE, MAX,MIN.Concentration,MAX.Concentration
  
  a=as.numeric(unname(max))
  
  b=as.numeric(unname(slope))
  d=0 # min response
  ic50 = as.numeric(unname(ic50))
  min.conc.tested = as.numeric(unname(min.conc.tested))
  max.conc.tested = as.numeric(unname(max.conc.tested))
  Min.Conc<- log10(min.conc.tested*concn_scale) #
  Max.Conc<- max.conc.tested
  x2<-log10(Max.Conc*concn_scale)
  
  
  if(is.na(ic50)||is.na(b)||is.na(a)||is.na(Min.Conc)||is.na(Max.Conc)){
    dss<-NA
  }
  else if(isTRUE(ic50>=Max.Conc)){
    dss<-0
  }
  else if(isTRUE(b==0)){
    dss<-0
  }
  else{
    if(a>100){ a<-100  }
    if(isTRUE(b<0)){ b<--b  }
    c<-log10(ic50*concn_scale)
    if(a>y){
      if(y!=0){
        x1<-(c - ((log(a-y)-log(y-d))/(b*log(10))))
        if(isTRUE(x1 < Min.Conc)){x1<-Min.Conc}
        else if(isTRUE(x1 > x2)){x1<-x2}
      }
      else {x1<-Min.Conc}
      
      # This is a logistic function used in Dotmatics.com
      # y = d+(a-d)/(1+10^(b*(c-x)))
      #inverse function
      # x = c - ((log(a-y)-log(d-y))/(b*log(10)))
      
      int_y=(((((a-d)*log(1+10^(b*(c-x2))))/(b*log(10)))+a*x2)-((((a-d)*log(1+10^(b*(c-x1))))/(b*log(10)))+a*x1)) - (y*(x2-x1))
      
      total_area<-(x2-Min.Conc)*(100-y)
      
      if(DSS.type==1){
        norm_area<-((int_y/total_area)*100)#DSS1
      }
      if(DSS.type==2){
        #	if(a>100){a<-100}
        norm_area<-((int_y/total_area)*100)/log10(a)#DSS2 #AUC1
        if(isTRUE(norm_area > 50)){ norm_area <- 0}
      }
      if(DSS.type==3){
        #	if(a>100){a<-100}
        norm_area<-((int_y/total_area)*100)*(log10(100)/log10(a))*((x2-x1)/(x2-Min.Conc)) #DSS3 #AUC5
      }
      if(isTRUE(norm_area < 0|norm_area > 100)){
        dss<-0
      }else{
        dss<-round(norm_area,digits=4)}
    } else {dss<-0}
  }
  return (dss)
}

CALC_IC50_EC50_DSS <- compiler::cmpfun(function(i, drug_wells_, xpr_tbl, DSS_typ, readoutCTX = F, path = ""){
  
  tryCatch({
    #gc(T);
    TEC50 = ifelse(readoutCTX, "TC50", "EC50"); drug_wells = drug_wells_[i,];
    
    #find indices of wells with drugs
    idx_filt <- xpr_tbl$ID %in% drug_wells$ID #& xpr_tbl$ProductName %in% drug_wells$ProductName
    #extract inhib. and viab. for wells with drugs in current plate
    inhibition = inhibition2 <- xpr_tbl$inhibition_percent[idx_filt]; viability2 = 100 - inhibition2; # with 2 used for ploting of real values.
    # if there are identical values in inhibition, add a bit noise
    if(all(inhibition <= 0)) inhibition <- rep(0, length(inhibition))
    if(any(duplicated(inhibition))) inhibition <- seq(from = 0, length.out = length(inhibition), by = 0.01) + inhibition;
    
    viability = 100-inhibition; believe_ = T;
    
    # extract concentrations, unique drug names and product ids for wells with drugs in current plate
    dose <- as.numeric(xpr_tbl$Concentration[idx_filt])
    drug_name <- unique(as.character(xpr_tbl$ID)[idx_filt])
    product_id <- unique(as.character(xpr_tbl$ID)[idx_filt])
    
    #combine the data and sort by dose.
    mat_tbl <- data.frame(inhibition,dose,logconc = log10(dose),viability, inhibition2, viability2)
    mat_tbl <- mat_tbl[order(mat_tbl[,"dose"]),]
    
    #print(paste0(product_id,",   ",drug_name));print(mat_tbl);
    
    
    if(DSS_typ == "AUC"){
      
      mat_tbl$indexx = 1:nrow(mat_tbl)
      model = approx(x = mat_tbl$indexx, y = mat_tbl$inhibition2, xout = seq(1,nrow(mat_tbl),length.out = 100), method="linear")
      loess_fit <- loess(y ~ x, model)
      model$y = predict(loess_fit)
      
      AUC <- round(sum(diff(model$x) * (head(model$y,-1)+tail(model$y,-1)))/2 / 5, 2)
      
      perInh <- t(matrix(mat_tbl[,"inhibition"],dimnames=
                           list(paste0(rep("D", length(mat_tbl[,"inhibition"])), 1:length(mat_tbl[,"inhibition"])))))
      IC50_dataframe <- data.frame(ID=product_id,DRUG_NAME=drug_name,ANALYSIS_NAME="IC50", IC50="",SLOPE="",MAX=max(model$y),MIN=min(model$y), 
                                   Min.Conc.tested=min(mat_tbl$dose),Max.Conc.tested=max(mat_tbl$dose), IC50_std_error="", perInh, GRAPH="", 
                                   DSS = as.numeric(AUC), sDSS = "", SE_of_estimate = "")
      EC50_dataframe <- data.frame(ID=product_id,DRUG_NAME=drug_name,ANALYSIS_NAME="EC50", EC50="",SLOPE="",MAX=100-max(model$y),MIN=100-min(model$y), 
                                   Min.Conc.tested=min(mat_tbl$dose),Max.Conc.tested=max(mat_tbl$dose),TEC50_std_error="",perInh, GRAPH="", 
                                   DSS = as.numeric(AUC), sDSS = "", SE_of_estimate = "")
      
      # icpl <- ggplot2::ggplot(mat_tbl, aes(indexx, inhibition2)) + scale_x_continuous(breaks=1:nrow(mat_tbl),labels=mat_tbl$dose) +
      #   geom_point(color = "blue", size = 2.8) + geom_line(data = data.frame(x = model$x, y = model$y), aes(x, y), color="blue", size = 0.8) +
      #   ggtitle(paste0(drug_name," (AUC:",AUC,")\n"))+
      #   theme(legend.title = element_text(size = 9)) + theme_bw() + labs(y = "response", x = "conc(nM)")  +
      #   theme(plot.background = element_rect(fill = "transparent",colour = NA),
      #         panel.background =element_rect(fill = "transparent",colour = NA), plot.title = element_text(hjust = 0.5))
      
      
      # graphics.off()
      # filename_ = file.path(path,"IC50", paste0(product_id,"_IC50_curve_drug.png"))
      # png(filename = filename_,width=190,height=190, bg = "transparent")
      # print(icpl)
      # dev.off()
      
      # ecpl <- ggplot2::ggplot(mat_tbl, aes(indexx, viability2)) + scale_x_continuous(breaks=1:nrow(mat_tbl),labels=mat_tbl$dose) +
      #   geom_point(color = "blue", size = 2.8) + geom_line(data = data.frame(x = model$x, y = 100-model$y), aes(x, y), color="blue", size = 0.8) +
      #   ggtitle(paste0(drug_name," (AUC:",AUC,")\n"))+
      #   theme(legend.title = element_text(size = 9)) + theme_bw() + labs(y = "response", x = "conc(nM)")  +  ylim(-25, 125) +
      #   theme(plot.background = element_rect(fill = "transparent",colour = NA),
      #         panel.background =element_rect(fill = "transparent",colour = NA), plot.title = element_text(hjust = 0.5))
      #
      
      # graphics.off()
      # png(filename = file.path(path, "Curve_fits", TEC50, paste0(product_id,"_", TEC50,"_curve_drug.png")),width=190,height=190, bg = "transparent")
      # print(ecpl)
      # dev.off()
      
      #TEC50base64 <- gsub("\r?\n|\r", " ", base64::img(filename_))
      
      #return list with 3 nodes - 1 row for IC50 table and 1 row for EC50 table and EC50 image in base64
      list(IC50_dataframe, EC50_dataframe, #TEC50base64,
           believe_ = T,
           AUC)
      
    }else if(nrow(mat_tbl) <= 3 || (length(unique(mat_tbl$dose)) <= 3)){
      
      print("Less than 3 rows... skipping...")
      NULL
    } else {
      
      #############################
      #############    IC50
      
      estimate_param <- tryCatch({drm(inhibition ~ logconc, data = mat_tbl, fct = LL.4(fixed = c(NA, NA, NA,NA),names = c("SLOPE","MIN","MAX","IC50")),logDose=10,control = drmc(errorm = F))},
                                 warning=function(w){drm(inhibition ~ logconc, data = mat_tbl, fct = L.4(fixed = c(NA, NA, NA,NA), names = c("SLOPE","MIN","MAX","IC50")),logDose=10)},
                                 error=function(e){drm(inhibition ~ logconc, data = mat_tbl, fct = L.4(fixed = c(NA, NA, NA,NA), names = c("SLOPE","MIN","MAX","IC50")),logDose=10)})
      # (extract and name coefficients)
      coef_estim <- coef(estimate_param); names(coef_estim) <- c("SLOPE","MIN","MAX","IC50")
      # see http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4696819/
      coef_estim["SLOPE"] <- coef_estim["SLOPE"]*-1
      
      # if curve decreases or IC50 is higher than max (i.e. IC50 is "outlier"), set IC50 to max conc.
      coef_estim["IC50"] <- ifelse(coef_estim["MAX"]<=coef_estim["MIN"] | coef_estim["IC50"]>max(mat_tbl$dose,na.rm=T), max(mat_tbl$dose,na.rm=T),coef_estim["IC50"])
      # if IC50 is less than 0 set it to min. conc. and if even min. conc. < 0, then set IC50 to mean of all conc.
      coef_estim["IC50"] <- ifelse(coef_estim["IC50"]<0,min(mat_tbl$dose,na.rm=T),coef_estim["IC50"])
      coef_estim["IC50"] <- ifelse(coef_estim["IC50"]<0,mean(mat_tbl$dose,na.rm=T),coef_estim["IC50"])
      # similar to previous step but now compare log10(IC50) with log(min. conc.).
      coef_estim["IC50"] <- log10(coef_estim["IC50"])
      coef_estim["IC50"] <- ifelse(coef_estim["IC50"]<min(mat_tbl$logconc),max(mat_tbl$logconc),coef_estim["IC50"])
      # if all inhib. < 0 set IC50 to max. log. conc !!!!! not obvious why!
      coef_estim["IC50"] <- ifelse(all(mat_tbl$inhibition<0),max(mat_tbl$logconc,na.rm=T),coef_estim["IC50"])
      #(Trying to fix curves that need outlier kickout)
      coef_estim["MIN"] <- 0; coef_estim["MAX"] <- max(mat_tbl$inhibition,na.rm=T)
      #(Fix off minimums) Find lowest inhibition value. If it is not in (0:100), fix it whether to 0 or 99.
      min_lower <- ifelse(min(mat_tbl$inhibition,na.rm=T) > 0,min(mat_tbl$inhibition,na.rm=T),0)
      min_lower <- ifelse(min_lower >= 100,99,min_lower)
      #similar to previous step but for MAX
      coef_estim["MAX"] <- ifelse(coef_estim["MAX"]>100,100,coef_estim["MAX"])
      coef_estim["MAX"] <- ifelse(coef_estim["MAX"]<0,100,coef_estim["MAX"])
      #max_lower and max_upper - lower and upper bounds for 'nl2sol' algorithm in nonlinear least-squares
      max_lower <- ifelse(max(mat_tbl$inhibition,na.rm=T)>100,coef_estim["MAX"],max(mat_tbl$inhibition,na.rm=T))
      max_lower <- ifelse(max_lower < 0,coef_estim["MAX"],max(mat_tbl$inhibition,na.rm=T))
      max_lower <- ifelse(max_lower < 0,0,max_lower)
      max_lower <- ifelse(max_lower > 100,100,max_lower)
      #(Fix upper maximum for negative slopes)
      run_avg <- caTools::runmean(mat_tbl$inhibition, 10)
      max_upper <- ifelse(any(run_avg[-nrow(mat_tbl)]>run_avg[nrow(mat_tbl)]),max(mat_tbl$inhibition[run_avg>run_avg[nrow(mat_tbl)]]),coef_estim["MAX"])
      max_upper <- ifelse(any(mat_tbl$inhibition > max_upper),mean(mat_tbl$inhibition[mat_tbl$inhibition > max_upper])+5,max_upper)
      max_upper <- ifelse(max_upper < 0,coef_estim["MAX"],max_upper)
      max_upper <- ifelse(max_upper > 100,100,max_upper) #coef_estim["MAX"]
      max_upper <- ifelse(max_lower > max_upper,coef_estim["MAX"],max_upper)
      # left it as it was, just rewritten a bit (ALEKS). not clear how values 25, 60 and 5 are chosen.
      mean_inh_last = mean(tail(mat_tbl$inhibition,2),na.rm=T)
      if(mean_inh_last < 60) {
        if(mean_inh_last > 25) coef_estim["IC50"] <- mean(mat_tbl$logconc,na.rm=T)
        else if(mean_inh_last < 25) coef_estim["IC50"] <- max(mat_tbl$logconc,na.rm=T)}
      if(mean(mat_tbl$inhibition[1:3],na.rm=T)<5) coef_estim["IC50"] <- max(mat_tbl$logconc,na.rm=T)
      #add a bit of positive noise to MAX if it is the same as MIN.
      if(unname(coef_estim["MIN"]) == unname(coef_estim["MAX"])) coef_estim["MAX"] <- coef_estim["MAX"] + 0.001
      
      #adaptive nonlinear Least-Squares algorithm NL2SOL to estimate parameters.
      nls_result_ic50_old <- function(){
        tryCatch({
          nls(inhibition ~ MIN + (MAX - MIN)/ (1 + (10^(SLOPE * (IC50 - logconc)))), data=mat_tbl, algorithm="port", start=list(SLOPE=1,MIN=coef_estim["MIN"][[1]],MAX=coef_estim["MAX"][[1]],IC50=coef_estim["IC50"][[1]]), lower=list(SLOPE=0.1,MIN=0,MAX=max_lower, IC50=min(mat_tbl$logconc)), upper=list(SLOPE=2.5,MIN=0,MAX=max_upper, IC50=max(mat_tbl$logconc)),control=list(warnOnly=T,minFactor = 1/2048))
        }, error = function(e) {
          
          # allows higher residual sum-of-squares
          minpack.lm::nlsLM(inhibition ~ MIN + (MAX - MIN)/ (1 + (10^(SLOPE * (IC50 - logconc)))), data=mat_tbl,
                            start=list(SLOPE=1, MIN=coef_estim["MIN"][[1]],MAX=coef_estim["MAX"][[1]],IC50=coef_estim["IC50"][[1]]),
                            lower=c(SLOPE=0.1, MIN=0,MAX=max_lower, IC50=min(mat_tbl$logconc)),
                            upper=c(SLOPE=2.5, MIN=0,MAX=max_upper, IC50=max(mat_tbl$logconc)))
        })
      }
      
      # IC50 first
      nls_result_ic50 <- nls_result_ic50_old();
      
      # IC50 second
      nls_result_ic50_2 <- tryCatch({
        # allows higher residual sum-of-squares
        nls(inhibition ~ MIN + (MAX - MIN)/ (1 + (10^(SLOPE * (IC50 - logconc)))), data=mat_tbl, algorithm="port",  start=list(SLOPE=1,MIN=coef_estim["MIN"][[1]],MAX=coef_estim["MAX"][[1]], IC50=median(mat_tbl$logconc)),lower=list(SLOPE=0.1,MIN=0,MAX=max_lower, IC50=min(mat_tbl$logconc)),upper=list(SLOPE=2.5,MIN=0,MAX=max_upper, IC50=max(mat_tbl$logconc)),control=list(warnOnly=T,minFactor = 1/2048))
      },warning = function(w) {
        nls_result_ic50_old()
      },error = function(e) {
        nls_result_ic50_old()
      })
      
      #element (4, 4) is zero, so the inverse cannot be computed
      nls_result_ic50 = tryCatch({summary(nls_result_ic50); nls_result_ic50},error=function(e){nls_result_ic50_2})
      
      #Calculate the standard error scores
      sumIC50 = list(summary(nls_result_ic50), summary(nls_result_ic50_2))
      
      ic50std_resid <- round(sqrt(sum((sumIC50[[1]]$residuals)^2)/(length(sumIC50[[1]]$residuals)-1)),1);
      ic50std_resid2 <- round(sqrt(sum((sumIC50[[2]]$residuals)^2)/(length(sumIC50[[2]]$residuals)-1)),1);
      
      # continue with the best
      switch_ = which.min(c(ic50std_resid, ic50std_resid2))
      nls_result_ic50 = list(nls_result_ic50, nls_result_ic50_2)[[switch_]]
      
      
      #if SLOPE <= 0.2, decrease IC50, change lower bound for SLOPE to 0.1 and repeat.
      if(coef(nls_result_ic50)["SLOPE"] <= 0.2)
      {
        if(mean_inh_last > 60)
          coef_estim["IC50"] <- min(mat_tbl$logconc,na.rm=T)
        nls_result_ic50 <- nls(inhibition ~ MIN + (MAX - MIN)/ (1 + (10^(SLOPE * (IC50 - logconc)))), data=mat_tbl, algorithm="port",start=list(SLOPE=1, MIN=unname(coef_estim["MIN"]),MAX=unname(coef_estim["MAX"]),IC50=unname(coef_estim["IC50"])),lower=list(SLOPE=0.1,MIN=min_lower,MAX=max_lower,IC50=min(mat_tbl$logconc)),upper=list(SLOPE=2.5, MIN=0,MAX=max_upper, IC50=max(mat_tbl$logconc)),control=list(warnOnly=T,minFactor = 1/2048))
      }
      
      #Calculate the standard error scores
      sumIC50 = summary(nls_result_ic50);
      ic50std_Error <- sumIC50$coefficients["IC50","Std. Error"]; #tec50std_Error <- sumTEC50$coefficients["TEC50","Std. Error"]
      ic50std_resid <- round(sqrt(sum((sumIC50$residuals)^2)/(length(sumIC50$residuals)-1)),1);
      max_signal <- max(mat_tbl$dose,na.rm=T); min_signal <- min(mat_tbl$dose,na.rm=T)
      
      #############################
      #############   Final modification & STD error
      
      #prepare final data and convert IC50 back from log scale (inverse)
      coef_ic50 <- coef(nls_result_ic50)[c("IC50", "SLOPE","MAX","MIN")]; coef_ic50["IC50"] <- 10^coef_ic50["IC50"]
      #(Fix ic50 for curves in wrong direction)
      coef_ic50["IC50"] <- ifelse(coef_ic50["SLOPE"]<0,max_signal,coef_ic50["IC50"])
      #(Fix based on MAX)
      coef_ic50["IC50"] <- ifelse(coef_ic50["MAX"]<0,max_signal,coef_ic50["IC50"])
      coef_ic50["IC50"] <- ifelse(coef_ic50["MAX"]<10,max_signal,coef_ic50["IC50"])
      coef_ic50["MAX"] <- ifelse(coef_ic50["MAX"]<0,0,coef_ic50["MAX"])
      #(Fix over sensitive drugs)
      coef_ic50["IC50"] <- ifelse(all(c(max(mat_tbl$inhibition,na.rm=T),min(mat_tbl$inhibition,na.rm=T))>50),min_signal,coef_ic50["IC50"])
      
      # for ploting
      x <- seq(min(mat_tbl$logconc),max(mat_tbl$logconc), length=100)
      yic <- predict(nls_result_ic50, data.frame(logconc=x))
      auc <- MESS::auc(x,yic)
      
      ##average replicates
      mat_tblCp <- mat_tbl[, c("inhibition", "dose")]
      cols_ <- colnames(mat_tblCp)[!grepl("inhibition", colnames(mat_tblCp))] # columns which should be equal to average PI
      X <- as.data.table(mat_tblCp)
      mat_tblCp <- as.data.frame(X[,list(inhibition = mean(inhibition)),cols_], stringAsFactors = !1)
      
      
      perInh <- t(matrix(mat_tblCp[,"inhibition"],dimnames=
                           list(paste0(rep("D", length(mat_tblCp[,"inhibition"])), 1:length(mat_tblCp[,"inhibition"])))))
      
      coef_tec50 = coef_ic50;
      coef_tec50["IC50"] <- ifelse(coef_tec50["MAX"] > 25, coef_tec50["IC50"], max(mat_tbl$dose,na.rm=T))
      if(readoutCTX){
        names(coef_tec50) <- c("TC50","SLOPE","MAX","MIN"); ytec <- yic; perViaTox <- perInh;
      } else{
        names(coef_tec50) <- c("EC50","SLOPE","MAX","MIN");
        coef_tec50["SLOPE"] = -1 * coef_tec50["SLOPE"]; # min - 0, max - 77 in ec50 it is max - 100, min - 23
        tmp = coef_tec50["MAX"]; coef_tec50["MAX"] = 100 - coef_tec50["MIN"]; coef_tec50["MIN"] = 100 - tmp; ytec <- 100 - yic;
        perViaTox <- 100 - perInh;
      }
      
      
      #############################
      #############    DSS
      
      dss_score <- round(as.numeric(dss(coef_ic50["IC50"],coef_ic50["SLOPE"],coef_ic50["MAX"],min_signal,max_signal, DSS.type=as.integer(DSS_typ))),1);
      dss_score1 <- round(as.numeric(dss(coef_ic50["IC50"],coef_ic50["SLOPE"],coef_ic50["MAX"],min_signal,max_signal, DSS.type=1)),1);
      dss_score2 <- round(as.numeric(dss(coef_ic50["IC50"],coef_ic50["SLOPE"],coef_ic50["MAX"],min_signal,max_signal, DSS.type=2)),1);
      dss_score3 <- round(as.numeric(dss(coef_ic50["IC50"],coef_ic50["SLOPE"],coef_ic50["MAX"],min_signal,max_signal, DSS.type=3)),1);
      coef_ic50 <- c(coef_ic50,Min.Conc.tested=min_signal,Max.Conc.tested=max_signal,IC50_std_error=ic50std_Error)
      coef_tec50 <- c(coef_tec50,Min.Conc.tested=min_signal,Max.Conc.tested=max_signal,TEC50_std_error=ic50std_Error)
     
      ####
      # Absolute IC50
      xIC50ABS <- seq(min(mat_tbl$logconc),max(mat_tbl$logconc)*15, length=5000)
      yicIC50ABS <- predict(nls_result_ic50, data.frame(logconc=xIC50ABS))
      if(all(yicIC50ABS < 50)) coef_ic50ABS= Inf else coef_ic50ABS = 10**xIC50ABS[which.min(abs(yicIC50ABS - 50))]
      ####
      
      #dataframe for IC50
      IC50_dataframe <- data.frame(ID=product_id,DRUG_NAME=drug_name,ANALYSIS_NAME="IC50", t(as.matrix(coef_ic50)), perInh,
                                   GRAPH=NA, DSS = as.numeric(dss_score), sDSS = "", SE_of_estimate = as.numeric(ic50std_resid),AUC=auc)
      #dataframe for EC50
      TEC50_dataframe <- data.frame(ID=product_id,DRUG_NAME=drug_name,ANALYSIS_NAME=TEC50,t(as.matrix(coef_tec50)), perViaTox,
                                    GRAPH=NA, DSS = as.numeric(dss_score), sDSS = "", SE_of_estimate = as.numeric(ic50std_resid),AUC=auc)
      
      #round by 2 dex. all the numeric colums
      numeric_cols <- sapply(IC50_dataframe, is.numeric); IC50_dataframe[,numeric_cols] <- round(IC50_dataframe[,numeric_cols],1)
      numeric_cols <- sapply(TEC50_dataframe, is.numeric); TEC50_dataframe[,numeric_cols] <- round(TEC50_dataframe[,numeric_cols],1)
      
      # plot IC50
      #mat_tbl$inhibition = xpr_tbl$inhibition_percent[idx_filt]; # if we have all values < 0, they will be replaced
      #mat_tbl$viability = 100 - mat_tbl$inhibition;  # we are replacing them back here.
      # icpl <- ggplot2::ggplot(mat_tbl, aes(logconc, inhibition2)) + scale_x_continuous(breaks=mat_tbl$logconc,labels=mat_tbl$dose) +
      #   geom_point(color = "blue", size = 2.8) + geom_line(data = data.frame(x = x, y = yic), aes(x, yic), color="blue", size = 0.8) +
      #   geom_vline(xintercept = log10(coef_ic50["IC50"]), colour="grey", size = 0.8) + ggtitle(paste0(strtrim(drug_name, 15)," (dss:",dss_score,")\n"))+
      #   theme_bw() + labs(y = "% inhibition", x = "conc(nM)")  +  ylim(-25, 125) +
      #   geom_text(mapping=aes(x2,y2,label = text2), data=data.frame(x2=log10(coef_ic50["IC50"])*0.95, y2=115, text2="IC50"), color="grey", parse=T) +
      #   theme(plot.background = element_rect(fill = "transparent",colour = NA),
      #         panel.background =element_rect(fill = "transparent",colour = NA), plot.title = element_text(hjust = 0.5, size = 12.5))
      
      # graphics.off()
      # filename_ = file.path(path, "IC50", paste0(product_id,"_IC50_curve_drug.png"))
      # png(filename = filename_,width=190,height=190, bg = "transparent")
      # print(icpl)
      # dev.off()
      
      if(readoutCTX) aes_ <- aes(logconc, inhibition2) else aes_ <- aes(logconc, viability2)
      # ecpl <- ggplot2::ggplot(mat_tbl, aes_) + scale_x_continuous(breaks=mat_tbl$logconc,labels=mat_tbl$dose) +
      #   geom_point(color = "blue", size = 2.8) + geom_line(data = data.frame(x = x, y = ytec), aes(x, ytec), color="blue", size = 0.8) +
      #   geom_vline(xintercept = log10(coef_tec50[TEC50]), colour="grey", size = 0.8) + ggtitle(paste0(strtrim(drug_name, 15)," (dss:",dss_score,")\n"))+
      #   theme_bw() + labs(y = ifelse(readoutCTX, "% toxicity", "% viability"), x = "conc(nM)") +  ylim(-25, 125) +
      #   geom_text(mapping=aes(x2,y2,label = text2), data=data.frame(x2=log10(coef_tec50[TEC50])*0.95, y2=115, text2=TEC50), color="grey", parse=T) +
      #   theme(plot.background = element_rect(fill = "transparent",colour = NA),
      #         panel.background =element_rect(fill = "transparent",colour = NA),
      #         plot.title = element_text(hjust = 0.5, size = 12.5))
      
      #browser()
      # graphics.off()
      # png(filename = file.path(path, "Curve_fits", paste0(product_id,"_", TEC50,"_curve_drug.png")),width=190,height=190, bg = "transparent")
      # print(ecpl)
      # dev.off()
      
      # TEC50base64 <- gsub("\r?\n|\r", " ", base64::img(filename_))
      
      # browser()
      # check believe
      if(IC50_dataframe$DSS > 1){
        if(IC50_dataframe$SE_of_estimate > 40) believe_ = F; # if SE > 40
        #if(IC50_dataframe[[paste0("D", length(mat_tbl$inhibition))]] - IC50_dataframe[[paste0("D", length(mat_tbl$inhibition)-1)]] < -25) believe_ = F; # if last is less effective than penultimate
        if(sumIC50$residuals[[1]] > 20 && sumIC50$residuals[[2]] > 20) believe_ = F; # if residuals for first 2 points more than 10
        if(max(sumIC50$residuals) > 30) believe_ = F; # if any of individual points deviates more than 30
      }
      
      resid_ = as.numeric(sumIC50$residuals); resp_ = as.numeric(perInh); cond = 0;
      if(sum(abs(resid_)>15)>1) {believe_ = !1; cond = 2};if(tail(resp_,2)[[1]]-tail(resp_,1) > 7 && tail(resp_,1)-tail(resp_,3)[[1]] > 7) {believe_ = !1; cond = 3}
      if(any(abs(resid_)>10) && (ic50std_resid>10)) {believe_ = !1; cond = 5};if(ic50std_resid>100) {believe_ = !1; cond = 6}
      if(sum(abs(resid_)>10)>2){believe_ = !1; cond = 8};if(sum(abs(resid_)>10)>1 && any(abs(resid_)>15)) {believe_ = !1; cond = 9};
      if(sum(resp_ < 2) >= (length(resp_)-1)) {believe_ = !0; cond = 1}
      #print(paste0(product_id,"   ",believe_,"   ",cond));
      Readout=ifelse(readoutCTX, "Toxicity", "Viability");
      cp_table=mat_tbl
      cp_table$Readout=Readout
      #cp_table$screen=screen
      cp_table$product_id=product_id
      cp_table$drug_name=drug_name
      cp_table$DSS=dss_score
      cp_table$IC50=IC50_dataframe$IC50
      cp_table$EC50=unlist(TEC50_dataframe[TEC50])
      dose_response_fit=data.frame(#screen=screen,
        product_id=product_id,drug_name=drug_name,x=x,yic=yic,ytec=ytec,Readout=Readout)
      #return list with 3 nodes - 1 row for IC50 table and 1 row for EC50 table and EC50 image in base64
      list(IC50_dataframe, TEC50_dataframe, #cp_table,#TEC50base64,
           believe_,
           dss_score,
           dss_score1,
           dss_score2,
           dss_score3,
           coef_ic50["IC50"],
           coef_tec50["TC50"])
    }
  }, error = function(e) {
    print(paste0("error in ", product_id, ", in ", xpr_tbl$screen_id[idx_filt][[1]]));
    print(e);
  })
})

