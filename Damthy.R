## ----packages-------------------------------------------------------
library(ggplot2)
library(factoextra)
library(dplyr)
library(pryr)
library(magrittr)
library(broom)
library(parallel)
library(data.table)
library(pbapply)

## ----intialize------------------------------------------------------
GSEnumber <- readline(prompt = "What's the GSE number you are analyzing? ")


## ----create class structure-----------------------------------------
setClass("Damthy", slots = list(name = "character"),
         prototype = list(name = "For infinuim"))
  setClass("rawD",
           contains = "Damthy",
           slots = list(idat = "matrix",
                        m_pval = "data.frame",
                        GPL = "data.frame",
                        series = "data.frame"))
  setClass("pp",
           contains = "Damthy",
           slots = list(name = "character"), 
           prototype = list(name = "processed data"))
    setClass("raw.o",
             contains = "pp",
             slots = list(raw.m = "matrix",
                          raw.p = "matrix",
                          raw.s = "data.frame",
                          raw.g = "data.frame",
                          beta_dtr = "list"))
    setClass("qc.o",
             contains = "pp",
             slots = list(name = "character",
                          m = "matrix",
                          s = "data.frame",
                          beta_dtr = "list",
                          svd.o = "list",
                          pca.o = "list",
                          ctf.o = "list"),
             prototype = list(name = "quality control"))
    setClass("he.o",
               contains = "pp",
               slots = list(name = "character",
                            m = "matrix",
                            s = "data.frame",
                            svd.o = "list",
                            ctf.o = "list",
                            pca.o = "list"))


## ----instantiation function-----------------------------------------
f.rawD <- function(idat=0,GPL,m_pval=0,series,...){
  if (idat == 0){
    print(paste0("you got the processed data with p values of ",GSEnumber))
    rawD <- new("rawD", m_pval = m_pval,GPL = GPL, series = series)
  }else if(m_pval == 0) {
    rawD <- new("rawD", idat = idat, GPL = GPL, series = series)
    print(paste0("you got the idat data of ",GSEnumber))
  }
  return(rawD)
}

f.raw.o <- function(rawD,...){
  raw.o <- new("raw.o", raw.s = rawD@series)
  print("Tidy matrix with p values, you need to check the first column of the matrix (CpGs)!")
  DNAm.beta <- rawD@m_pval[, seq(from=2, to=ncol(rawD@m_pval)-1, by= 2)]
  rownames(DNAm.beta) <- rawD@m_pval[, 1]
  DNAm.Pval <- rawD@m_pval[, seq(from=3, to=ncol(rawD@m_pval), by= 2)]
  colnames(DNAm.Pval) <- colnames(DNAm.beta)
  rownames(DNAm.Pval) <- rawD@m_pval[, 1]
  raw.o@raw.p <- as.matrix(DNAm.Pval)
  # cpg qc
  cpg.Pval <- pbapply(DNAm.Pval, 1, function(x) (sum(x > 0.01)/length(x))) #计算每一个CpG p>0.01的比例
  raw.o@raw.m <- as.matrix(DNAm.beta[(cpg.Pval < (1/ncol(DNAm.beta))), ]) #比例大于样本数倒数的要去掉（有一个p>0.001就要去掉）
  del <- nrow(DNAm.beta)-nrow(raw.o@raw.m)
  cat("There are",del,"CpGs removed","\n")
  raw.o@raw.g <- rawD@GPL[match(rownames(raw.o@raw.m),rawD@GPL$ID),] #match GPL to raw.m
  print(paste0("you got the processed data of ",GSEnumber))
  raw.o@beta_dtr <- list(beta_dtp(raw.o))
  gc()
  return(raw.o)
}

f.qc.o <- function(raw.o,...){
qc.o <- new("qc.o", name = paste0(GSEnumber,"_all"))
choice <- readline(prompt = "Do you need to adjust type 2 probe bias? y/n")
if (choice == "y"){
  result <-bmiq(raw.o)
  qc.o@m <- result[[1]]
  design.v <- result[[2]]
}else if (choice == "n"){
  design.v <- tryCatch({
    raw.o@raw.g$Infinium_Design_Type %>% gsub("II",2,.) %>% gsub("I",1,.) %>% as.numeric()
  }, error = function(e) { 
    print("Check the GPL file:GPL$Infinium_Design_Type,II,I") 
  })
  qc.o@m <-raw.o@raw.m
}else{
  print("please input a right word:y/n")
}
qc.o@beta_dtr <- list(beta_dtp(qc.o,design.v))
return(qc.o)
}

f.he.o <- function(qc.o,...){
  he.o <- new("he.o", name = paste0(GSEnumber,"_healthy"))
  idx <- which(qc.o@s$desease == 0)
  he.o@m <- qc.o@m[,idx]
  he.o@s <- qc.o@s[idx,-which(colnames(qc.o@s) == "desease")]
  return(he.o)
}

## ----coverage of probe--------------------------------------------------
setGeneric("coverage",function(obj,...){
  standardGeneric("coverage")
})
setMethod("coverage","raw.o",function(obj,...){
  library(ggplot2)
  library(pbapply)
  library(gridExtra)
  # coverage in each probe
  cg <- pbapply(obj@raw.p,1,function(x) {sum(x < 0.01)/length(x)})
  pc.df <- data.frame(rownames(obj@raw.p),cg)
  colnames(pc.df)<-c("probe","coverage")
  # coverage in each sample
  sg <- pbapply(obj@raw.p,2,function(x) {sum(x < 0.01)/length(x)})
  ps.df <- data.frame(colnames(obj@raw.p),sg)
  colnames(ps.df)<-c("sample","coverage")
  plot1 <- ggplot(pc.df, aes(x = probe, y = coverage)) + 
    geom_point() +
    theme(axis.text.x = element_blank())+
    xlab("Probe") + ylab("Coverage")+
    ggtitle("Probe Coverage Plot")
  plot2 <- ggplot(ps.df, aes(x = sample, y = coverage)) + 
    geom_point() +
    theme(axis.text.x = element_blank())+
    xlab("Sample") + ylab("Coverage") + 
    ggtitle("Sample Coverage Plot")
  my_grid <- grid.arrange(plot1, plot2, layout_matrix=rbind(c(1,1),c(2,2)))
  ggsave("coverage.pdf",my_grid,width = 8,height = 8)
  
  })


## ----beta_dtp-----------------------------------------------------------
setGeneric("beta_dtp",function(obj,...){
  standardGeneric("beta_dtp")
})
setMethod("beta_dtp","raw.o",function(obj,...){
  density2 <- density(obj@raw.m[which(obj@raw.g$Infinium_Design_Type == "II"),1])
  density1 <- density(obj@raw.m[which(obj@raw.g$Infinium_Design_Type == "I"),1])
  density_data <- data.frame(
  beta1 = density1$x,
  beta2 = density2$x,
  density1 = density1$y,
  density2 = density2$y
  )
  fplot <- ggplot(density_data, aes(x = beta1, y = density1)) + 
  geom_line(aes(x = beta1, y = density1, color = "Type I")) + 
  geom_line(aes(x = beta2, y = density2, color = "Type II")) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
  labs(x = "Beta", y = "Density", title = "Density plot of two types of probe") +
  scale_color_manual(name = "Probe Type", values = c("Type I" = "blue", "Type II" = "red"), labels = c("Type I", "Type II")) +
  theme(legend.position = c(0.95, 0.95),
        legend.justification = c("right", "top"),
        legend.margin = margin(-10, 0, 0, -10),
        legend.box.margin = margin(0, 0, 10, 10),
        legend.direction = "vertical",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1))
  return(fplot)
})

setMethod("beta_dtp","qc.o",function(obj,design.v,...){
  density2 <- density(obj@m[which(design.v == 2),1])
  density_data <- data.frame(
    beta = density2$x,
    density = density2$y
  )
  fplot <- raw.o@beta_dtr[[1]] + 
    geom_line(data = density_data,aes(x = beta, y = density, color = "Type II-BMIQ")) + 
    scale_color_manual(name = "Probe Type", values = c("Type I" = "blue", "Type II" = "red", "Type II-BMIQ" = "black"), labels = c("Type I", "Type II", "Type II-BMIQ")) + labs(x = "Beta", y = "Density", title = "Density plot of two types of probe after BMIQ")
  return(fplot)
})
#1# ----Dealing with the duplicates------------------------------------
setGeneric("Del_dup",function(obj,...){
  standardGeneric("Del_dup")
})
setMethod("Del_dup","raw.o",function(obj,...){
  library(dendextend)
  c <- readline(prompt = "Have you dealt the duplicates? T/F")
  library(broom)
  library(isva)
  # outliers
  identify_outliers <- function(x) {
    # 计算所有数的绝对偏差
    abs_dev <- abs(x - median(x))
    # 计算绝对中位差
    mad <- median(abs_dev)
    # 根据绝对中位差，将可能的离群值定义为离群区域的2倍以上的值
    # 输出所有可能的离群值
    outliers <- x[abs_dev / mad > 2]
    return(outliers)
  }
  ----------------------------------
    m <- obj@raw.m - rowMeans(obj@raw.m)
  print("Estimating dimensionality")
  n <- EstDimRMT(m)$dim
  print("Run SVD")
  svd.m <- svd(m)
  v_top <- svd.m$v[,1:n]
  #  rownames(v_top) <- colnames(obj@m)
  cor_m <- as.dist(1-cor(t(v_top)))
  clu.o <- hclust(cor_m, method = "complete")
  nm <- obj@raw.m
  ns <- obj@raw.s
  duplicates <- table(factor(obj@raw.s$`individual id:ch1`))[table(factor(obj@raw.s$`individual id:ch1`)) > 1]
  
  par(mfrow=c(1,length(duplicates)),mar=c(4,4,4,4))
  pdf("dup_clust.pdf", width=10,height=5)
  for (i in names(duplicates)){
    idx <- which(obj@raw.s$`individual id:ch1` %in% i)
    pos <- as.numeric(clu.o$order) %in% idx
    #delete and average----------
    if (length(idx) == 2){
      dis <- abs(match(idx[1],clu.o$order)-match(idx[2],clu.o$order))
      if(dis == 1){
        mean_col <- colMeans(obj@raw.m[,idx])
        nm <- cbind(nm,mean_col)
        colnames(nm)[ncol(nm)] <- rownames(obj@raw.m)[idx[1]]
        ns <- rbind(ns,obj@raw.s[idx[1],])
      }else{
        s1<-sum(obj@raw.p[,ids[1]]<0.01)
        s2<-sum(obj@raw.p[,ids[2]]<0.01)
        if (s1>s2){
          nm <- cbind(nm,obj@raw.m[,ids[1]])
          ns <- rbind(ns,obj@raw.s[idx[1],])
        }else {
          nm <- cbind(nm,obj@raw.m[,ids[2]])
          ns <- rbind(ns,obj@raw.s[idx[2],])
        }
      }
    }else {
      a <- identify_outliers(sort(match(idx,clu.o$order)))
      dx <- clu.o$order[a]
      px <- setdiff(idx,dx)
      mean_col <- colMeans(obj@raw.m[,px])
      nm <- cbind(nm,mean_col)
      colnames(nm)[ncol(nm)] <- rownames(obj@raw.m)[px[1]]
      ns <- rbind(ns,obj@raw.s[px[1],])
    }
    #---------------------------------------
    labels <- rep("", length(clu.o$order))
    labels[pos] <- as.character(idx)
    dend <- as.dendrogram(clu.o)
    dend <- color_branches(dend,k=13)
    dend %>% set("labels", labels) %>%  set("labels_cex",0.7) %>% plot
    
  }
  ddx <- which(obj@raw.s$`individual id:ch1` %in% names(duplicates))
  nm <- nm[,-ddx]
  ns <- ns[-ddx,]
  dev.off()
  return(list(nm,ns))
  gc()
})

#2# ----bmiq-----------------------------------------------------------
setGeneric("bmiq",function(obj,...){
  standardGeneric("bmiq")
})
setMethod("bmiq","raw.o",function(obj,...){
source("~/code/packages/BMIQ_1.6.R")
tmp.m <- obj@raw.m
num_cpus <- detectCores()-2
tmp.lv <- split(tmp.m, col(tmp.m))
design.v <- tryCatch({
  obj@raw.g$Infinium_Design_Type %>% gsub("II",2,.) %>% gsub("I",1,.) %>% as.numeric()
}, error = function(e) { 
  print("Check the GPL file:GPL$Infinium_Design_Type,II,I") 
})
print("Start BMIQ")
tmp.m <-tryCatch({
result_list <- mclapply(tmp.lv, function(x) {
  bmiq.o <- BMIQ(x, design.v)
  return(bmiq.o$nbeta)
}, mc.cores = num_cpus)
do.call(cbind, result_list)
}, warning = function(w) {
  message("Warning:", w$message)
  return(NULL)
}, error = function(e) {
  message("Error:", e$message)
  return(NULL)
}
)
# for(s in 1:ncol(tmp.m)){
#   beta.v <- obj@raw.m[,s]
#   bmiq.o <- BMIQ(beta.v, design.v, sampleID=s)
#   tmp.m[,s] <- bmiq.o$nbeta
#   print(paste("Done BMIQ for sample ",s,sep=""))
# }
rownames(tmp.m) <- rownames(obj@raw.m)
colnames(tmp.m) <- colnames(obj@raw.m)
print("BMIQ is over")
lv <- list(tmp.m,design.v)
return(lv)
})

#3# ----CTF------------------------------------------------------------
setGeneric("CTF",function(obj,...){
  standardGeneric("CTF")
})
setMethod("CTF","pp",function(obj,ref,wth=0.4,useW=TRUE,type="850k",...){
library(EpiSCORE)
avSIM.m <- constAvBetaTSS(obj@m, type=type)
estF.o <- wRPC(avSIM.m, ref=ref, useW=useW, wth=wth, maxit=500)
fplot <- boxplot(estF.o$estF, main = "Estimation of cell-type fractions ", xlab = "Cell types", ylab = "Fraction")
lv <- list(estF.o$estF,fplot)
return(lv)
})



#4# ----svd_lm---------------------------------------------------------
setGeneric("lm_svd",function(obj,...){
  standardGeneric("lm_svd")
})
setMethod("lm_svd","pp",function(obj,...){
  if (readline(prompt = "Are you ready to rum lm with a neat tibble of sample features? T/F")){
  library(broom)
  library(isva)
  m <- obj@m - rowMeans(obj@m)
  print("Estimating dimensionality")
  n <- EstDimRMT(m)$dim
  sam_var <- bind_cols(obj@ctf.o[[1]],obj@s)
  print("Run SVD")
  if (length(obj@svd.o)==0){
  svd.m <- svd(m)
  }else {
  svd.m <- obj@svd.o[[2]] 
  }
  v_top <- svd.m$v[,1:n]
  pval.m <- matrix(ncol = n,nrow = ncol(sam_var))
  print("Run lm")
  for (i in 1:ncol(sam_var)){
    print(i)
    if (class(sam_var[[i]]) != "factor"){
      pval <- apply(v_top, 2,function(x) {tidy(summary(lm(x ~ sam_var[[i]])))[[5]][2]})
    }else {
      pval <- apply(v_top, 2,function(x) {glance(summary(lm(x ~ sam_var[[i]])))[[5]][1]})
    }
    pval.m[i,] <- pval
  }
  colnames(pval.m) <- paste0("PC-",1:n)
  rownames(pval.m) <- colnames(sam_var)
  list <- list(pval.m,svd.m)
  return(list)
  }
})


#5# ----p_heatmap------------------------------------------------------
setGeneric("figure",function(obj,...){
  standardGeneric("figure")
})
setMethod("figure","pp",function(obj,...){
  c <- readline(prompt = "Have you dealt the duplicates? T/F")
  if (c){
  fV.v <- obj@svd.o[[2]]$d^2/sum(obj@svd.o[[2]]$d^2)
  plot.new()
  name <- paste0("SVDsummary_",obj@name)
  pdf(paste0(name,".pdf"),width=8,height=8)
  layout(matrix(1:2,nrow=2),heights=c(1,1))
  par(mar=c(4,6,2,1))
  plot(fV.v[1:ncol(obj@svd.o[[1]])],pch=23,type="b",col="red",ylab="fracV",xlab="Top-PC",main="")
  my_colors <- c("#3C2E2B", "#C92742", "#ECF659", "#FAC7EE", "#FFFFFF")
  # Define breaks for color legend
  my_breaks <- c(-300,-50,-15,-5,log10(0.05),0); #EPiSCORE  chosed by plot(m)
  # breaks.v <- c(-82,-25,-5,-3,log10(0.05),0); #EpiDISH
  
  # Create heat map without clustering
#  try(heatmap.2(log10(obj@svd.o[[1]]), dendrogram = "none", Rowv = FALSE, Colv = FALSE, col = my_colors, breaks = my_breaks, key = FALSE, key.title = NA, margins = c(5,10),trace = "none",labRow = "",labCol = ""))
  par(mar = c(4, 22, 2, 1))
  image(x=1:ncol(obj@svd.o[[1]]),y=1:nrow(obj@svd.o[[1]]),z=log10(t(obj@svd.o[[1]])),col=my_colors,breaks=my_breaks,xlab="",ylab="",axes=FALSE,asp=2);
  axis(1,at=1:ncol(obj@svd.o[[1]]),labels=paste("PC-",1:ncol(obj@svd.o[[1]]),sep=""),las=2)
  axis(2,at=1:nrow(obj@svd.o[[1]]),labels=rownames(obj@svd.o[[1]]),las=2)
  # Add legend
  legend("bottomright", legend = c("p<1e-50", "p<1e-15", "p<1e-5", "p<0.05", "p=ns"), fill = my_colors, bty = "o",
         cex = 0.8,
         pt.cex = 1)
  dev.off()}
})


#6# ------PCA ----------------------------------------------------------

setGeneric("PCA",function(obj,...){
  standardGeneric("PCA")
})

setMethod("PCA", "pp", function(obj, s, ...) {
  if (length(obj@pca.o)==0){
#    matrix <- scale(t(obj@m), center = TRUE, scale = TRUE)
    matrix <- obj@m
    pca_o <- prcomp(matrix, scale = FALSE)
    pca_data <- as.data.frame(pca_o$rotation[, 1:s])
  } else {
    pca_o <- obj@pca.o[[1]]
    pca_data <- as.data.frame(obj@pca.o[[1]]$rotation[, 1:s])
  }
  n <- readline(prompt = "Which feature do you want to mark?")
  idx <- NULL
  fill_values <- NULL
  test <- 0
  par(mfrow=c(1,(s-1)),mar=c(4,4,4,4))
  name <- paste("PCAsummary_", obj@name, "_", n, sep = "")
  pdf(paste0(name, ".pdf"), width=6,height=6)
  if (n == "sex") {
    idx <- which(obj@s$sex == 0)
    fill_values <- c("male" = "gray", "Female" = "brown")
  } else if (n == "disease") {
    idx <- which(obj@s$disease == 0)
    fill_values <- c("Bad" = "gray", "Healthy" = "blue")
  } else if (n %in% colnames(obj@ctf.o[[1]])) {
    test <- 1
    nf <- obj@ctf.o[[1]][,n]
    for (i in 1:(ncol(pca_data))) {
    # Pearson 相关性系数
    r <- cor(nf, pca_data[,i], method = "pearson")
    p_value <- cor.test(nf, pca_data[,i])$p.value
    plot(nf,pca_data[,i], xlab = paste0("%",n), ylab = paste0("PC", i), pch = 16, col = "black", cex = 0.5)
    fit <- lm(pca_data[,i] ~ nf)
    abline(fit, col = "red")
    y_hat <- predict(fit) # 进行拟合预测
    y_upper <- predict(fit, interval = "confidence", level = 0.95)[, 2] # 计算上界
    y_lower <- predict(fit, interval = "confidence", level = 0.95)[, 3] # 计算下界
    polygon(c(nf, rev(nf)), c(y_upper, rev(y_lower)), col = "#BDBDBD", border = NA) # 绘制误差带
    text(x = min(nf)+0.1*min(nf), y = max(pca_data[,i]),labels = paste0("r=",format(r, digits = 3)),cex = 0.8)
    text(x = min(nf)+0.1*min(nf), y = max(pca_data[,i]) - 0.1*(max(pca_data[,i]) - min(pca_data[,i])),cex = 0.8, labels = paste0("p=",format.pval(p_value,digits = 3)))
    }
    
  }else {
    stop("Please input a valid feature!")
    test = 1
  }
if (test == 0){
  
  
  for (i in 1:(ncol(pca_data)-1)) {
    
    g <- ggplot(data = pca_data, aes_string(x = paste0("PC", i), y = paste0("PC", i+1))) +
      geom_point(shape = 21, aes(fill = names(fill_values)[1]), color = "black", size = 2) +
      geom_point(data = pca_data[idx, ], shape = 21, aes(fill = names(fill_values)[2]), color = "black", size = 2) +
      labs(x = colnames(pca_data)[i], y = colnames(pca_data)[i+1], fill = "") +
      scale_fill_manual(values = fill_values) +
      theme_classic() +
      theme(legend.position = "top", legend.justification = "right")
    
    print(g)
  }
}
  dev.off()
  list(pca_o)
  })  




































































