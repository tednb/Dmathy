source("~/code/Dmathy/Damthy.R")

# Build a S4 class for training set
setClass("train",slots = list(m="matrix",age="numeric",GSE="character",age_his="list"))
# Build a S4 class for testing set
setClass("test",slots = list(m="matrix",age="numeric",GSE="character",age_his="list"))


#1# Extracting needed data from pp S4 object
# pp is an object from Dmathy, seed is a random seed, dmct is a data matrix with rows as age-CpGs and columns as cell types
setGeneric("extract",function(obj,...){
  standardGeneric("extract")
})
setMethod("extract","pp",function(obj,seed,...){
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
    n <- round(sum(booler)*0.7)
    cat("sample:",n," ","total:",sum(booler),"\n","\n")
    set.seed(seed)
    m_se <- sample(m_idx[which(booler)],n)
    m_s <- append(m_s,m_se)
  }
  m_t <- setdiff(m_idx,m_s)
  # training set and its age
  m_tr <- m[,m_s]
  age_tr <- age[m_s]
  # Divide the histogram by a 5-year interval
  hist_tr <- hist(age_tr,xlab = "age",main = "histogram of training set", breaks = age_seq)
  train_obj <- new("train", m = m_tr, age = age_tr, GSE = obj@name, age_his = list(hist_tr))
  # testing set and its age
  m_te <- m[,m_t]
  age_te <- age[m_t]
  # break the testing set into 5 age groups
  hist_te <- hist(age_te,xlab = "age",main = "histogram of testing set", breaks = age_seq)
  test_obj <- new("test", m = m_te, age = age_te, GSE = obj@name, age_his = list(hist_te))
  return(list(train_obj,test_obj))
})



#2# Using the training set to build a clock by elastic net regression
setGeneric("clock",function(obj,...){
  standardGeneric("clock")
})
setMethod("clock","train",function(obj,dmct,...){
  m_tr <- obj@m
  age_tr <- obj@age
  m_trains <- list()
  for (i in 1:ncol(dmct)) {
    cell <- colnames(dmct)[i]  # cell types
    idx <- match(rownames(dmct)[which(dmct[,i] != 0)],rownames(m_tr))
    idx <- idx[!is.na(idx)]
    m_trains[[cell]] <- m_tr[idx,]
  }
  # Build a clock for each cell type
  clocks <- list()
  for (i in 1:length(m_trains)){
    cell <- names(m_trains)[i]
    cat("cell:",cell,"\n")
    m <- t(m_trains[[i]])
    # Build a clock by elastic net regression using glmnet package
    require(glmnet)
    model <- glmnet(m,age_tr,family="gaussian",alpha = 0.5)
    # Extract all the lambda values and its corresponding RMSEs
    lambda <- model$lambda
    # Calculate RMSE for each lambda value
    rmse <- c()
    for (j in 1:length(lambda)){
      age_pre <- predict(model,newx = m,s = lambda[j])
      rmse <- append(rmse,sqrt(mean((age_pre-age_tr)^2)))
    }
    # Draw a plot of RMSEs against lambda values and output a PDF file
    pdf(paste0(cell,"_lambda_rmse.pdf"))
    plot(lambda,rmse,type = "l",xlab = "lambda",ylab = "RMSE")
    dev.off()
    # Find the lambda value with the minimum RMSE
    lambda_min <- lambda[which.min(rmse)]
    # Build a clock by elastic net regression using the lambda value with the minimum RMSE
    clocks[[cell]] <- glmnet(m,age_tr,family="gaussian",alpha = 0.5,lambda = lambda_min)
  }
  return(clocks)
})
#3# Using the clock to predict the age of a new bulk sample or sorted for validation
setGeneric("validate",function(obj,...){
  standardGeneric("validate")
})
setMethod("validate","test",function(obj,clocks,...){
  m_te <- obj@m
  age_te <- obj@age
  if (grepl("GSE",obj@GSE)){
    m_tests <- list()
    for (i in 1:length(clocks)) {
      cell <- names(clocks)[i]  # cell types
      print(cell)
      rowns<-clocks[[i]]$beta@Dimnames[[1]]
      idx <- match(rowns,rownames(m_te))
      idx <- idx[!is.na(idx)]
      m_tests[[cell]] <- m_te[idx,]
    }
    age_pre <- list()
    for (i in 1:length(m_tests)){
      cell <- names(m_tests)[i]
      m <- t(m_tests[[i]])
      # Predict the age of each cell type's testing set.
      age_pre[[cell]] <- predict(clocks[[cell]],newx = m)
      # Calculate the PCC between predicted age and real age
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
      # Fit a linear model to the predicted age and real age
      #fit <- lm(age_pre[[cell]]~age_te)
      # Output PDF files for each cell type
      pdf(paste0(cell," clock for dataset"," ",obj@GSE,".pdf"))
      # Plot the predicted age and real age with a 45 degree dashed line and a mark of correlation and RMSE
      plot(age_pre[[cell]],age_te,pch = 16,col = "blue",ylab = paste0(cell_all," age (years)"),xlab = "Chronological age (years)",main = paste0(cell_all," clock for dataset"," ",obj@GSE),xlim = c(10,80),ylim = c(10,80))
      abline(0,1,lty = 2)
      text(20,75,paste0("r = ",format(r, digits = 3),"\n","p = ",format(p, digits = 1),"\n","MedAE: ",format(MAE, digits = 2),"years"))
      # Add a linear model line, black and solid.
      #abline(fit,col = "black",lwd = 2)
      dev.off()
    }
    return(age_pre)
  }else if(grepl("sorted",obj@name)){
    # Using specific clock to predict the age of this kind of sorted sample
    if (grepl("hep",obj@name)){
      cell_type <- "Hep"}
    if(cell_type == "Hep"){
      cell_all <- "Hepatocyte"
    }else if (cell_type == "Chol") {
      cell_all <- "Cholangiocyte"
    }else {
      cell_all <- cell_type
    }
    idx <- match(clocks$cell_type$beta@Dimnames[[1]],rownames(m_te))
    idx <- idx[!is.na(idx)]
    age_pre <- predict(clocks[[cell_type]],newx = t(m_te[idx,]))
    cor <- cor(age_pre,age,method = "pearson")
    p <- cor.test(age_pre,age)$p.value
    MAE <- median(abs(age_pre-age))
    # Output PDF files for this cell type
    pdf(paste0(cell_type,"_sorted_predict.pdf"))
    # Plot the predicted age and real age with a 45 degree dashed line and a mark of correlation and RMSE
    plot(age_pre,age,pch = 16,col = "black",ylab = paste0(cell_all," age (years)"),xlab = "Chronological age (years)",main = paste("Predict",cell_all,"age","by",cell_all,"clock",sep = " "),xlim = c(10,80),ylim = c(10,80))
    abline(0,1,lty = 2)
    text(20,75,paste0("r = ",format(cor, digits = 3),"\n","p = ",format(p, digits = 1),"\n","MedAE: ",format(MAE, digits = 2),"years"))
    dev.off()
  }
})

#4# Gnerate a test object using a validation set
setGeneric("f_testset",function(obj,...){
  standardGeneric("f_testset")
})
setMethod("f_testset","pp",function(obj,...){
  m_te <- obj@m
  age_te <- obj@s$age
  GSE <- obj@name
  hist_te <- hist(obj@s$age,xlab = "age",main = "histogram of validation set")
  test_obj <- new("test", m = m_te, age = age_te, GSE = obj@name, age_his = list(hist_te))
  return(test_obj)
})

#5# Generate a trian object using a validation set
setGeneric("f_trainset",function(obj,...){
  standardGeneric("f_trainset")
})
setMethod("f_trainset","train",function(obj,test,...){
  # idx of common CpGs between training set and validation set
  idx <- match(rownames(test@m),rownames(obj@m))
  idx <- idx[!is.na(idx)]
  # Build a new training set
  m_tr <- obj@m[idx,]
  train_obj <- new("train", m = m_tr, age = obj@age, GSE = obj@GSE, age_his = obj@age_his)
  return(train_obj)
})

#6# pipeline for a validation set
pip_pre <- function(qc.o,trainset,dmct,...){
  gse <- readline(prompt = "What is your set's GSE number?")
  qc.o@name <- gse
  vali <- f_testset(qc.o) #new predicting set
  trainss <- f_trainset(trainset,vali)# new training set
  clockss <-clock(trainss,dmct) # new clock
  age_pre <- validate(vali ,clockss)
  return(age_pre)
}