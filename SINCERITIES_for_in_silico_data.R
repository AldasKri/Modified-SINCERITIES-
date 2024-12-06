####################
# PACKAGES required:
# kSamples
# glmnet
# ppcor
# pracma
# R.matlab
####################
library(kSamples)
library(glmnet)
library(ppcor)
library(R.matlab)
library(transport)

## *** Data loading ***
SIGN <- 1
  for (x in 1:7){
    SIGN <- 1
    AUROC_mean <- numeric(0)
    AUROC_sd <- numeric(0)
    AUPR_mean <- numeric(0)
    AUPR_sd <- numeric(0)
      mat <- readMat('In silico single cell data/20_nets_20genes_8UNEVENtime_sigma01B_no_initial_points2.mat')
      time <- as.vector(mat$time.points)
      numGENES <- as.vector(mat$n)
      AUROC <- vector()
      AUPR <- vector()
      for (numEXAMPLES in 1:dim(mat$networks)[3]) {
        data_time_series <- mat$data.tot.array[[numEXAMPLES]][[1]]
        singleCELLdata <- list()
        for (i in 1:mat$num.time.points) {
          singleCELLdata[[i]] <- data_time_series[,i,]
        }
        genes <- vector(length=numGENES)
        for (i in 1:numGENES) {
          genes[i] <- sprintf('Gene %d',i)
        }
        totDATA <- matrix(nrow = 0, ncol = dim(data_time_series)[3])
        for (i in 1:mat$num.time.points) {
          totDATA <- rbind(totDATA,data_time_series[,i,])
        }
        
        DATA <- list(time=time, numGENES=numGENES, singleCELLdata=singleCELLdata, genes=genes, totDATA=totDATA)
        
        ## *** SINCERITIES ***
        
        library(kSamples)
        library(glmnet)
        library(ppcor)
        SINCERITITES <- dget("SINCERITIES functions/SINCERITIES.R")
        
        
        result <- SINCERITITES(DATA,distance=x,method = 6,noDIAG = 1,SIGN = 1)
        adj_matrix <- result$adj_matrix
        
        ## *** Performance Evaluation ***
        
        #Gold standard GRN
        a <- mat$networks[,,numEXAMPLES]
        a[row(a)==col(a)] <- 0
        if(SIGN==0){
          a[which(a!=0)] <- 1
        }
        
        #Final ranked list, AUROC and AUPR
        adj_matrix <- adj_matrix/max(adj_matrix)
        library(pracma)
        auc_from_ranks_TC_sign <- dget("SINCERITIES functions/auc_from_ranks_TC_sign.R")
        AUCresult <- auc_from_ranks_TC_sign(adj_matrix,a,1000)
        
        AUROC[numEXAMPLES] <- AUCresult$AUROC
        AUPR[numEXAMPLES] <- AUCresult$AUPR
        final_ranked_predictions <- dget("SINCERITIES functions/final_ranked_predictions.R")
        table <- final_ranked_predictions(adj_matrix,DATA$genes,SIGN=0,fileNAME=sprintf('results4insilicoNETWORK%d',numEXAMPLES),saveFile = TRUE)
      }
      AUC <- cbind(AUROC,AUPR)
      m <- apply(AUC,2,mean)
      s <- apply(AUC,2,std)
      
      AUROC_mean <- cbind(AUROC_mean,m[[1]])
      AUROC_sd <- cbind(AUROC_sd,s[[1]])
      AUPR_mean <- cbind(AUPR_mean,m[[2]])
      AUPR_sd <- cbind(AUPR_sd,s[[2]])
  
    distance_names <- c('KS','CvM','AD','KL','JS','Hell','Wass')
    print(distance_names[x])
    print(paste("AUROC mean:",round(mean(AUROC_mean),3)," AUROC sd:",round(mean(AUROC_sd),3)))
    print(paste("AUPR mean:",round(mean(AUPR_mean),3)," AUPR sd",round(mean(AUPR_sd),3)))
  }
  