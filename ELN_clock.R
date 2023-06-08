source("~/code/Dmathy/Damthy.R")
# Build a S4 class for training set
setClass("train",slots = list(m="matrix",age="numeric",s="data.frame",ctf="matrix",GSE="character",age_his="list"))
# Build a S4 class for testing set
setClass("test",slots = list(m="matrix",age="numeric",GSE="character",age_his="list"))


library(glmnet)
#1# Extracting needed data from pp S4 object
# pp is an object from Dmathy, seed is a random seed, dmct is a data matrix with rows as age-CpGs and columns as cell types
setGeneric("extract",function(obj,...){
  standardGeneric("extract")
})
setMethod("extract","pp",function(obj,seed,...){
  s <- obj@s
  ctf <- obj@ctf.o[[1]]
  m <- obj@m # A matrix with CpG names as rows and sample names as columns.
  age <- obj@s$age # The age of each sample
  m_idx <- seq(1,ncol(m),1) # A vector of row indices
  #Generate a sequence based on age range with a step size of 5 years
  age_seq <- seq(min(age),max(age),length.out = round((max(age)-min(age))/5)+1)
  ##Randomly sampling 70% samples in m in each age step of age_seq
  age_group <- list()
  for (i in 1:(length(age_seq)-1)){
    age_group <- append(age_group,list(c(age_seq[i],age_seq[i+1])))
  }
  #Looping the age groups and randomly sampling 70% samples in each age group
  m_s <- c()
  for (i in 1:length(age_group)){
    cat("age:",age_group[[i]],"\n")
    booler <- age>=age_group[[i]][1] & age<age_group[[i]][2]
    n <- round(sum(booler)*0.76) #Changing how many you want to take as training set
    cat("sample:",n," ","total:",sum(booler),"\n","\n")
    set.seed(seed)
    m_se <- sample(m_idx[which(booler)],n)
    m_s <- append(m_s,m_se)
  }
  m_t <- setdiff(m_idx,m_s)
  # training set and its age
  m_tr <- m[,m_s]
  age_tr <- age[m_s]
  s_tr <- s[m_s,]
  ctf_tr <- ctf[m_s,]
  # Divide the histogram by a 5-year interval
  hist_tr <- hist(age_tr,xlab = "age",main = "histogram of training set", breaks = age_seq)
  train_obj <- new("train", m = m_tr, age = age_tr, s = s_tr,ctf = ctf_tr,GSE = obj@name, age_his = list(hist_tr))
  # testing set and its age
  m_te <- m[,m_t]
  age_te <- age[m_t]
  # break the testing set into 5 age groups
  hist_te <- hist(age_te,xlab = "age",main = "histogram of testing set", breaks = age_seq)
  test_obj <- new("test", m = m_te, age = age_te, GSE = obj@name, age_his = list(hist_te))
  return(list(train_obj,test_obj))
})


#2# Using CellDMC to find CTS age-DMCs in training set
setGeneric("CTSAC",function(obj,...){
  standardGeneric("CTSAC")
})
setMethod("CTSAC","train",function(obj,idx=NULL,covs=NULL,n=100,...){
  library(EpiDISH)
  if(is.null(idx)){
    ctf <- obj@ctf
    m <- obj@m
    age <- obj@age
  }else{
    ctf <- obj@ctf[idx,]
    m <- obj@m[,idx]
    age <- obj@age[idx]
    if(!is.null(covs)){
      covs <- covs[idx,]
    }
  }
  celldmc.o <- CellDMC(m, age ,ctf, 
                       cov.mod = covs, 
                       adjPMethod = "fdr",
                       adjPThresh = 0.05,
                       sort = FALSE,
                       mc.cores = n)
  dmct<- celldmc.o$dmct
  return(dmct)
})

setMethod("CTSAC","pp",function(obj,idx=NULL,covs=NULL,n=100,...){
  library(EpiDISH)
  if(is.null(idx)){
    ctf <- obj@ctf.o[[1]]
    m <- obj@m
    age <- obj@s$age
  }else{
    ctf <- obj@ctf.o[[1]][idx,]
    m <- obj@m[,idx]
    age <- obj@s$age[idx]
    if(!is.null(covs)){
      covs <- covs[idx,]
    }
  }
  celldmc.o <- CellDMC(m, age ,ctf, 
                       cov.mod = covs, 
                       adjPMethod = "fdr",
                       adjPThresh = 0.05,
                       sort = FALSE,
                       mc.cores = n)
  dmct<- celldmc.o$dmct
  return(dmct)
})

#3# Using the training set to build a clock by elastic net regression
setGeneric("clock",function(obj,...){
  standardGeneric("clock")
})
setMethod("clock","train",function(obj,covs,seqs,nP=10,...){
  print("Start building liver clock...")
  library(EpiDISH)
  m_tr <- obj@m
  age_tr <- obj@age
  ctf <- obj@ctf
  idx <- seq(1,ncol(m_tr),1)
  n.idx <- idx
  bags <- list();
  for(p in 1:nP){
    bags[[p]] <- sample(idx,ncol(m_tr)/nP,replace=FALSE)
    idx <- setdiff(idx,bags[[p]])
  }
  
  #for each 9 partitions run CellDMC, record hepatocyte & chol. age-DMCs
  hepDMCT.lm <- list()
  cholDMCT.lm <- list()
  test.li <- list()
  train.li <- list()
  dmcts <- list()
  
  for(p in 1:nP){
    test.idx <- bags[[p]]
    train.idx <- setdiff(n.idx,test.idx)
    test.li[[p]] <- test.idx
    train.li[[p]] <- train.idx
    dmcts[[p]] <- CTSAC(obj = obj,idx = train.idx,covs = covs)[,c(2,4)]
    hepDMCT.lm[[p]] <- rownames(dmcts[[p]])[which(dmcts[[p]][,2] != 0)]
    cholDMCT.lm[[p]] <- rownames(dmcts[[p]])[which(dmcts[[p]][,1] != 0)]
    cat("Got: ",length(hepDMCT.lm[[p]])," hep age-DMCs","\n","Got: ",length(cholDMCT.lm[[p]])," chol age-DMCs","\n")
    
  }
  predH.m <- matrix(NA,nrow=length(seqs),ncol=ncol(m_tr))
  predC.m <- matrix(NA,nrow=length(seqs),ncol=ncol(m_tr))
  # cross validation
  for(p in 1:nP){
    hep.idx <- match(hepDMCT.lm[[p]],rownames(m_tr))
    chol.idx <- match(cholDMCT.lm[[p]],rownames(m_tr))
    
    model_chol <- glmnet(t(m_tr[chol.idx,train.li[[p]]]),y=age_tr[train.li[[p]]],alpha=1,standardize=TRUE,lambda=seqs)
    model_hep <- glmnet(t(m_tr[hep.idx,train.li[[p]]]),y=age_tr[train.li[[p]]],alpha=1,standardize=TRUE,lambda=seqs)
    
    match(test.li[[p]],n.idx) -> maptoN.idx
    predH.m[,maptoN.idx] <- t(predict(model_hep,newx=t(m_tr[hep.idx,test.li[[p]]])))
    match(test.li[[p]],n.idx) -> maptoN.idx
    predC.m[,maptoN.idx] <- t(predict(model_chol,newx=t(m_tr[chol.idx,test.li[[p]]])))
  }
  
  rmseH.v <- vector();
  rmseC.v <- vector();
  for(li in 1:nrow(predH.m)){
    rmseH.v[li] <- sqrt(mean((predH.m[li,] - age_tr)^2));
  }
  la_hep <- rev(seqs)[which.min(rmseH.v)]
  for(li in 1:nrow(predC.m)){
    rmseC.v[li] <- sqrt(mean((predC.m[li,] - age_tr)^2))
  } 
  la_chol <- rev(seqs)[which.min(rmseC.v)]
  if(la_chol == min(seqs) | la_hep == min(seqs)){
    stop("The training has failed this time, please run it again!")
  }
  age_hep <- predC.m[which.min(rmseH.v),]
  pdf("OptimizationCurve-Hep.pdf",width=5,height=3);
  par(mar=c(4,4,2,1))
  plot(rev(seqs),rmseH.v,ylab="RMSE",xlab="Lambda",type="l",main = "Hepatocyte")
  # add a vertical line at the optimal lambda
  abline(v=la_hep,lty=2,col="red")
  dev.off()
  
  pdf("OptimizationCurve-Chol.pdf",width=5,height=3)
  par(mar=c(4,4,2,1))
  plot(rev(seqs),rmseC.v,ylab="RMSE",xlab="Lambda",type="l",main = "Cholangiocyte")
  # add a vertical line at the optimal lambda
  abline(v=la_chol,lty=2,col="red")
  dev.off()
  age_chol <- predC.m[which.min(rmseC.v),]
  
  # calculate the correlation between predicted age and chronological age
  r <- cor(age_tr,age_hep,method = "pearson")
  # calculate the p-value
  p <- cor.test(age_tr,age_hep,method = "pearson")$p.value
  # Calculate the Median Absolute Error between predicted age and real age
  mae <- median(abs(age_tr - age_hep))
  # draw scatter plot of predicted age vs. chronological age
  
  pdf("pre_vs_chro.pdf",width=12,height=6)
  par(mfrow=c(1,2))
  plot(age_tr,age_hep,pch = 16,col = "blue",ylab = "Predicted age (years)",xlab = "Chronological age (years)",main = "Hepatocyte",xlim = c(10,95),ylim = c(10,95))
  abline(0,1,col="red")
  text(20,90,paste0("r = ",format(r, digits = 3),"\n","p = ",format(p, digits = 1),"\n","MedAE: ",format(mae, digits = 2),"years"))
  
  r <- cor(age_tr,age_chol,method = "pearson")
  # calculate the p-value
  p <- cor.test(age_tr,age_chol,method = "pearson")$p.value
  # Calculate the Median Absolute Error between predicted age and real age
  mae <- median(abs(age_tr - age_chol))
  # draw scatter plot of predicted age vs. chronological age
  plot(age_tr,age_chol,pch = 16,col = "blue",ylab = "Predicted age (years)",xlab = "Chronological age (years)",main = "Cholangiocyte",xlim = c(10,95),ylim = c(10,95))
  abline(0,1,col="red")
  text(20,90,paste0("r = ",format(r, digits = 3),"\n","p = ",format(p, digits = 1),"\n","MedAE: ",format(mae, digits = 2),"years"))
  dev.off()
  # clock
  dmct_all <- CTSAC(obj = obj,covs = covs)[,c(2,4)]
  id_h <- match(rownames(dmct_all)[which(dmct_all[,2] != 0)],rownames(m_tr))
  id_c <- match(rownames(dmct_all)[which(dmct_all[,1] != 0)],rownames(m_tr))
  clocks <- list()
  clocks[["Hep"]] <- glmnet(t(m_tr[id_h,]),age_tr,family="gaussian",alpha = 1,lambda = la_hep)
  clocks[["Chol"]] <- glmnet(t(m_tr[id_c,]),age_tr,family="gaussian",alpha = 1,lambda = la_chol)
  data<-list("dmct_all"=dmct_all,"age_pre"=list(age_hep,age_chol),"dmct_cell" = list(hepDMCT.lm,cholDMCT.lm))
  return(list(clocks,data))
  return(clocks)
})

setMethod("clock","pp",function(obj,covs,seqs,nP=10,...){
  print("Start building liver clock...")
  library(EpiDISH)
  m_tr <- obj@m
  age_tr <- obj@s$age
  ctf <- obj@ctf.o[[1]]
  idx <- seq(1,ncol(m_tr),1)
  n.idx <- idx
  bags <- list();
  for(p in 1:nP){
    bags[[p]] <- sample(idx,ncol(m_tr)/nP,replace=FALSE)
    idx <- setdiff(idx,bags[[p]])
  }
  #for each 9 partitions run CellDMC, record hepatocyte & chol. age-DMCs
  hepDMCT.lm <- list()
  cholDMCT.lm <- list()
  test.li <- list()
  train.li <- list()
  dmcts <- list()
  
  for(p in 1:nP){
    test.idx <- bags[[p]]
    train.idx <- setdiff(n.idx,test.idx)
    test.li[[p]] <- test.idx
    train.li[[p]] <- train.idx
    dmcts[[p]] <- CTSAC(obj = obj,idx = train.idx,covs = covs)[,c(2,4)]
    hepDMCT.lm[[p]] <- rownames(dmcts[[p]])[which(dmcts[[p]][,2] != 0)]
    cholDMCT.lm[[p]] <- rownames(dmcts[[p]])[which(dmcts[[p]][,1] != 0)]
    cat("Got: ",length(hepDMCT.lm[[p]])," hep age-DMCs","\n","Got: ",length(cholDMCT.lm[[p]])," chol age-DMCs","\n")
    
  }
  predH.m <- matrix(NA,nrow=length(seqs),ncol=ncol(m_tr))
  predC.m <- matrix(NA,nrow=length(seqs),ncol=ncol(m_tr))
  
  for(p in 1:nP){
    hep.idx <- match(hepDMCT.lm[[p]],rownames(m_tr))
    chol.idx <- match(cholDMCT.lm[[p]],rownames(m_tr))
    
    model_chol <- glmnet(t(m_tr[chol.idx,train.li[[p]]]),y=age_tr[train.li[[p]]],alpha=1,standardize=TRUE,lambda=seqs)
    model_hep <- glmnet(t(m_tr[hep.idx,train.li[[p]]]),y=age_tr[train.li[[p]]],alpha=1,standardize=TRUE,lambda=seqs)
    
    match(test.li[[p]],n.idx) -> maptoN.idx
    predH.m[,maptoN.idx] <- t(predict(model_hep,newx=t(m_tr[hep.idx,test.li[[p]]])))
    match(test.li[[p]],n.idx) -> maptoN.idx
    predC.m[,maptoN.idx] <- t(predict(model_chol,newx=t(m_tr[chol.idx,test.li[[p]]])))
    
  }
  
  
  rmseH.v <- vector();
  rmseC.v <- vector();
  for(li in 1:nrow(predH.m)){
    rmseH.v[li] <- sqrt(mean((predH.m[li,] - age_tr)^2))
  }
  for(li in 1:nrow(predC.m)){
    rmseC.v[li] <- sqrt(mean((predC.m[li,] - age_tr)^2))
  } 
  la_chol <- rev(seqs)[which.min(rmseC.v)]
  la_hep <- rev(seqs)[which.min(rmseH.v)]
  if(la_chol == min(seqs) | la_hep == min(seqs)){
    stop("The training has failed this time, please run it again!")
  }
  
  age_hep <- predC.m[which.min(rmseH.v),]
  pdf("OptimizationCurve-Hep.pdf",width=5,height=3);
  par(mar=c(4,4,2,1))
  plot(rev(seqs),rmseH.v,ylab="RMSE",xlab="Lambda",type="l",main = "Hepatocyte")
  # add a vertical line at the optimal lambda
  abline(v=la_hep,lty=2,col="red")
  dev.off()

  pdf("OptimizationCurve-Chol.pdf",width=5,height=3)
  par(mar=c(4,4,2,1))
  plot(rev(seqs),rmseC.v,ylab="RMSE",xlab="Lambda",type="l",main = "Cholangiocyte")
  # add a vertical line at the optimal lambda
  abline(v=la_chol,lty=2,col="red")
  dev.off()
  age_chol <- predC.m[which.min(rmseC.v),]
  
  # calculate the correlation between predicted age and chronological age
  r <- cor(age_tr,age_hep,method = "pearson")
  # calculate the p-value
  p <- cor.test(age_tr,age_hep,method = "pearson")$p.value
  # Calculate the Median Absolute Error between predicted age and real age
  mae <- median(abs(age_tr - age_hep))
  # draw scatter plot of predicted age vs. chronological age
  
  pdf("pre_vs_chro.pdf",width=12,height=6)
  par(mfrow=c(1,2))
  plot(age_tr,age_hep,pch = 16,col = "blue",ylab = "Predicted age (years)",xlab = "Chronological age (years)",main = "Hepatocyte",xlim = c(10,95),ylim = c(10,95))
  abline(0,1,col="red")
  text(20,90,paste0("r = ",format(r, digits = 3),"\n","p = ",format(p, digits = 1),"\n","MedAE: ",format(mae, digits = 2),"years"))
  
  r <- cor(age_tr,age_chol,method = "pearson")
  # calculate the p-value
  p <- cor.test(age_tr,age_chol,method = "pearson")$p.value
  # Calculate the Median Absolute Error between predicted age and real age
  mae <- median(abs(age_tr - age_chol))
  # draw scatter plot of predicted age vs. chronological age
  plot(age_tr,age_chol,pch = 16,col = "blue",ylab = "Predicted age (years)",xlab = "Chronological age (years)",main = "Cholangiocyte",xlim = c(10,95),ylim = c(10,95))
  abline(0,1,col="red")
  text(20,90,paste0("r = ",format(r, digits = 3),"\n","p = ",format(p, digits = 1),"\n","MedAE: ",format(mae, digits = 2),"years"))
  dev.off()
  
  dmct_all <- CTSAC(obj = obj,covs = covs)[,c(2,4)]
  id_h <- match(rownames(dmct_all)[which(dmct_all[,2] != 0)],rownames(m_tr))
  id_c <- match(rownames(dmct_all)[which(dmct_all[,1] != 0)],rownames(m_tr))
  clocks <- list()
  clocks[["Hep"]] <- glmnet(t(m_tr[id_h,]),age_tr,family="gaussian",alpha = 1,lambda = la_hep)
  clocks[["Chol"]] <- glmnet(t(m_tr[id_c,]),age_tr,family="gaussian",alpha = 1,lambda = la_chol)
  data<-list("dmct_all"=dmct_all,"age_pre"=list(age_hep,age_chol),"dmct_cell" = list(hepDMCT.lm,cholDMCT.lm))
  return(list(clocks,data))
  
})

#4# Using the clock to predict the age of a new bulk sample or sorted for validation
setGeneric("validate",function(obj,...){
  standardGeneric("validate")
})
setMethod("validate","pp",function(obj,clocks,gse,...){
  m_te <- obj@m
  age_te <- obj@s$age
  age_pre <- list()
  # predict the age of the test set
  #m_pre <- model.matrix(~.,data = data.frame(t(m_te)))


  for(i in 1:length(clocks)){
    cell <- names(clocks)[i]
    rowna <- clocks[[cell]]$beta@Dimnames[[1]]
    m1_new <- matrix(data = NA,ncol = ncol(m_te),nrow = length(rowna))
    rownames(m1_new) <- rowna
    colnames(m1_new) <- colnames(m_te)
    
    common_features <- intersect(rownames(m_te),rowna)
    idx_e <- match(common_features,rowna)
    idx_t <- match(common_features,rownames(m_te))
    m1_new[idx_e,] <- m_te[idx_t,]
    missing_features <- setdiff(rowna,common_features)
    varies_t<-rowna[clocks[[cell]]$beta@i]
    miss_cpg<-intersect(varies_t,missing_features)
    l <- length(miss_cpg)
    cat("The number of missing CpGs: ",l,"\n")
    idx_na <- setdiff(seq(1,length(rowna),1),idx_e)
    idx_miss <- match(miss_cpg,rowna)
    idx_0 <- setdiff(idx_na,idx_miss)
    m1_new[idx_0,] <- 0
    age_pre[[cell]] <- predict(clocks[[cell]],newx = t(m1_new),na.action=na.omit)
  }
  for (i in 1:length(age_pre)){
    cell <- names(age_pre)[i]
    r <- cor(age_pre[[cell]],age_te,method = "pearson")
    # Calculate the p value of PCC
    p <- cor.test(age_pre[[cell]],age_te,method = "pearson")$p.value
    # Calculate the Median Absolute Error between predicted age and real age
    MAE <- median(abs(age_pre[[cell]]-age_te))
    # Set the cell type name
    if(cell == "Hep"){
      cell_all <- "Hepatocyte"
    }else if (cell == "Chol") {
      cell_all <- "Cholangiocyte"
    }else {
      cell_all <- cell
    }
    pdf(paste0(cell," clock for dataset"," ",gse,".pdf"))
    # Plot the predicted age and real age with a 45 degree dashed line and a mark of correlation and RMSE
    plot(age_te,age_pre[[cell]],pch = 16,ylab = paste0(cell_all," age (years)"),xlab = "Chronological age (years)",main = paste0(cell_all," clock for dataset"," ",gse),xlim = c(10,95),ylim = c(10,95))
    abline(0,1,lty = 2)
    text(20,90,paste0("r = ",format(r, digits = 3),"\n","p = ",format(p, digits = 1),"\n","MedAE: ",format(MAE, digits = 2),"years"))
    # Add a linear model line, black and solid.
    #abline(fit,col = "black",lwd = 2)
    dev.off()
  }
  
  return(age_pre)
  
})





#8#
setGeneric("bench",function(obj,...){
  standardGeneric("bench")
})
setMethod("bench","train",function(obj,seqs,nP=10,...){
  print("Start building all clock...")
  m_tr <- obj@m
  age_tr <- obj@age
  ctf <- obj@ctf
  idx <- seq(1,ncol(m_tr),1)
  n.idx <- idx
  bags <- list();
  for(p in 1:nP){
    bags[[p]] <- sample(idx,ncol(m_tr)/nP,replace=FALSE)
    idx <- setdiff(idx,bags[[p]])
  }
  test.li <- list()
  train.li <- list()
  for(p in 1:nP){
    test.idx <- bags[[p]]
    train.idx <- setdiff(n.idx,test.idx)
    test.li[[p]] <- test.idx
    train.li[[p]] <- train.idx
  }
  pred.m <- matrix(NA,nrow=length(seqs),ncol=ncol(m_tr))
  for(p in 1:nP){
    model <- glmnet(t(m_tr[,train.li[[p]]]),y=age_tr[train.li[[p]]],alpha=1,standardize=TRUE,lambda=seqs)
    match(test.li[[p]],n.idx) -> maptoN.idx
    pred.m[,maptoN.idx] <- t(predict(model,newx=t(m_tr[,test.li[[p]]])))
  }
  rmse.v <- vector();
  for(li in 1:nrow(pred.m)){
    rmse.v[li] <- sqrt(mean((pred.m[li,] - age_tr)^2));
  }
  la <- rev(seqs)[which.min(rmse.v)]
  clocks <- list()
  pdf("OptimizationCurve-All.pdf",width=5,height=3);
  par(mar=c(4,4,2,1));
  plot(rev(seqs),rmse.v,ylab="RMSE",xlab="Lambda",type="l",main = "All");
  # add a vertical line at the optimal lambda
  abline(v=la,lty=2,col="red")
  dev.off()
  clocks[["All"]] <- glmnet(t(m_tr),age_tr,family="gaussian",alpha = 1,lambda = la)
  age_all <- pred.m[which.min(rmse.v),]
  pdf("all_vs_chro.pdf",width=6,height=6)
  
  r <- cor(age_tr,age_all,method = "pearson")
  # calculate the p-value
  p <- cor.test(age_tr,age_all,method = "pearson")$p.value
  # Calculate the Median Absolute Error between predicted age and real age
  mae <- median(abs(age_tr - age_all))
  # draw scatter plot of predicted age vs. chronological age
  plot(age_tr,age_all,pch = 16,col = "blue",ylab = "Predicted age (years)",xlab = "Chronological age (years)",main = "benchmarking",xlim = c(10,95),ylim = c(10,95))
  abline(0,1,col="red")
  text(20,90,paste0("r = ",format(r, digits = 3),"\n","p = ",format(p, digits = 1),"\n","MedAE: ",format(mae, digits = 2),"years"))
  dev.off()
  
  return(clocks)
})

setMethod("bench","pp",function(obj,seqs,nP=10,...){
  print("Start building all clock...")
  m_tr <- obj@m
  age_tr <- obj@s$age
  ctf <- obj@ctf.o[[1]]
  idx <- seq(1,ncol(m_tr),1)
  n.idx <- idx
  bags <- list()
  for(p in 1:nP){
    bags[[p]] <- sample(idx,ncol(m_tr)/nP,replace=FALSE)
    idx <- setdiff(idx,bags[[p]])
  }
  test.li <- list()
  train.li <- list()
  for(p in 1:nP){
    test.idx <- bags[[p]]
    train.idx <- setdiff(n.idx,test.idx)
    test.li[[p]] <- test.idx
    train.li[[p]] <- train.idx
  }
  pred.m <- matrix(NA,nrow=length(seqs),ncol=ncol(m_tr))
  for(p in 1:nP){
    model <- glmnet(t(m_tr[,train.li[[p]]]),y=age_tr[train.li[[p]]],alpha=1,standardize=TRUE,lambda=seqs)
    match(test.li[[p]],n.idx) -> maptoN.idx
    pred.m[,maptoN.idx] <- t(predict(model,newx=t(m_tr[,test.li[[p]]])))
  }
  if(any(is.na(pred.m))){
    idx_na <- which(apply(pred.m, 2, function(x) any(is.na(x))))
    pred.m <- pred.m[,-idx]
    age_tr <- age_tr[-idx]
    m_tr <- m_tr[,-idx]
  }
  rmse.v <- vector();
  for(li in 1:nrow(pred.m)){
    rmse.v[li] <- sqrt(mean((pred.m[li,] - age_tr)^2));
  }
  la <- rev(seqs)[which.min(rmse.v)]
  clocks <- list()
  pdf("OptimizationCurve-All.pdf",width=5,height=3);
  par(mar=c(4,4,2,1));
  plot(rev(seqs),rmse.v,ylab="RMSE",xlab="Lambda",type="l",main = "All");
  # add a vertical line at the optimal lambda
  abline(v=la,lty=2,col="red")
  dev.off()
  clocks[["All"]] <- glmnet(t(m_tr),age_tr,family="gaussian",alpha = 1,lambda = la)
  age_all <- pred.m[which.min(rmse.v),]
  pdf("all_vs_chro.pdf",width=6,height=6)

  r <- cor(age_tr,age_all,method = "pearson")
  # calculate the p-value
  p <- cor.test(age_tr,age_all,method = "pearson")$p.value
  # Calculate the Median Absolute Error between predicted age and real age
  mae <- median(abs(age_tr - age_all))
  # draw scatter plot of predicted age vs. chronological age
  plot(age_tr,age_all,pch = 16,col = "blue",ylab = "Predicted age (years)",xlab = "Chronological age (years)",main = "benchmarking",xlim = c(10,95),ylim = c(10,95))
  abline(0,1,col="red")
  text(20,90,paste0("r = ",format(r, digits = 3),"\n","p = ",format(p, digits = 1),"\n","MedAE: ",format(mae, digits = 2),"years"))
  dev.off()
  
  return(list(clocks,age_all))
})

setGeneric("colon_clock",function(obj,...){
  standardGeneric("colon_clock")
})
setMethod("colon_clock","pp",function(obj,covs,seqs,nP=10,...){
  print("Start building colon clock...")
  library(EpiDISH)
  m_tr <- obj@m
  age_tr <- obj@s$age
  ctf <- obj@ctf.o[[1]]
  idx <- seq(1,ncol(m_tr),1)
  n.idx <- idx
  bags <- list();
  for(p in 1:nP){
    bags[[p]] <- sample(idx,ncol(m_tr)/nP,replace=FALSE)
    idx <- setdiff(idx,bags[[p]])
  }
  #for each 9 partitions run CellDMC, record lm. age-DMCs
  lmDMCT.lm <- list()
  test.li <- list()
  train.li <- list()
  dmcts <- list()
  
  for(p in 1:nP){
    test.idx <- bags[[p]]
    train.idx <- setdiff(n.idx,test.idx)
    test.li[[p]] <- test.idx
    train.li[[p]] <- train.idx
    dmcts[[p]] <- CTSAC(obj = obj,idx = train.idx,covs = covs)[,4]
    lmDMCT.lm[[p]] <- names(dmcts[[p]])[which(dmcts[[p]] != 0)]
    cat("Got: ",length(lmDMCT.lm[[p]])," lm age-DMCs","\n")
    
  }
  pred.m <- matrix(NA,nrow=length(seqs),ncol=ncol(m_tr))
  for(p in 1:nP){
    lm.idx <- match(lmDMCT.lm[[p]],rownames(m_tr))
    model_lm <- glmnet(t(m_tr[lm.idx,train.li[[p]]]),y=age_tr[train.li[[p]]],alpha=1,standardize=TRUE,lambda=seqs)
    match(test.li[[p]],n.idx) -> maptoN.idx
    pred.m[,maptoN.idx] <- t(predict(model_lm,newx=t(m_tr[lm.idx,test.li[[p]]])))
  }
  if(any(is.na(pred.m))){
    idx_na <- which(apply(pred.m, 2, function(x) any(is.na(x))))
    pred.m <- pred.m[,-idx]
    age_tr <- age_tr[-idx]
    m_tr <- m_tr[,-idx]
  }
  
  rmse.v <- vector();
  for(li in 1:nrow(pred.m)){
    rmse.v[li] <- sqrt(mean((pred.m[li,] - age_tr)^2))
  }
  la_lm<- rev(seqs)[which.min(rmse.v)]
 
  if(la_lm == min(seqs)){
    stop("The training has failed this time, please run it again!")
  }
  
  age_lm <- pred.m[which.min(rmse.v),]
  pdf("OptimizationCurve-lm.pdf",width=5,height=3);
  par(mar=c(4,4,2,1))
  plot(rev(seqs),rmse.v,ylab="RMSE",xlab="Lambda",type="l",main = "Lymphocyte")
  # add a vertical line at the optimal lambda
  abline(v=la_lm,lty=2,col="red")
  dev.off()
  
  
  # calculate the correlation between predicted age and chronological age
  r <- cor(age_tr,age_lm,method = "pearson")
  # calculate the p-value
  p <- cor.test(age_tr,age_lm,method = "pearson")$p.value
  # Calculate the Median Absolute Error between predicted age and real age
  mae <- median(abs(age_tr - age_lm))
  # draw scatter plot of predicted age vs. chronological age
  
  pdf("pre_vs_chro.pdf",width=6,height=6)
  plot(age_tr,age_lm,pch = 16,col = "blue",ylab = "Predicted age (years)",xlab = "Chronological age (years)",main = "Lymphocyte",xlim = c(10,95),ylim = c(10,95))
  abline(0,1,col="red")
  text(20,90,paste0("r = ",format(r, digits = 3),"\n","p = ",format(p, digits = 1),"\n","MedAE: ",format(mae, digits = 2),"years"))
  dev.off()
  
  dmct_all <- CTSAC(obj = obj,covs = covs)[,4]
  id <- match(names(dmct_all)[which(dmct_all != 0)],rownames(m_tr))
  clocks <- list()
  clocks[["Lym"]] <- glmnet(t(m_tr[id,]),age_tr,family="gaussian",alpha = 1,lambda = la_lm)
  data<-list("dmct_all"=dmct_all,"age_pre"= age_lm,"dmct_cell" = lmDMCT.lm)
  return(list(clocks,data))

  
  
  })