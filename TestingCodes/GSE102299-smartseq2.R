set.seed(5000)
library(rpart)
library(rpart.plot)
library(AUC)
library(pROC)
library(party)
library(partykit)
library(ROCR)
library(evtree)
library(tree)
library(caret)
library(C50)
library(ggplot2)
library(BiocManager)
library(svMisc)

suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(MultiAssayExperiment))

setwd("~/Documents/R")

dataset <- readRDS("datasets/GSE102299-smartseq2.rds")
dataset <- updateObject(dataset)
colData(dataset)[,14]
head(colData(dataset))
table(colData(dataset)[,14])

(dataset <- experiments(dataset)[["gene"]])
dataset <- assays(dataset)[["TPM"]]
dataset <- t(dataset)
dataset[1:3,1:3]
dim(dataset)


#based on characteristics_ch1.4 factor (treatment: IL25 & treatment: IL33) #
dataset <- dataset[1:470,]
dim(dataset)
dataset[1:3,1:10]
row.names(dataset)[1:188] <- "1"
row.names(dataset)[189:470] <- "0"
dataset[1:3,1:10]

write.csv(dataset, "Test.csv")
dataset <- read.csv("Test.csv")
dataset<- data.frame(dataset)
#-------------------------------------------------------------------------

#Find out the lowset p_value among the genes in the training_set
gene_id <- c()
p_values <- c()
for(n in 2:ncol(dataset)) {
  wilcox <- wilcox.test(dataset[[n]] ~ dataset[[1]], data = dataset, exact = T, paired = FALSE)
  p_values[n-1] <- wilcox$p.value
  cat("\n", n)
}

gene_id <- c(1:(ncol(dataset)-1))
head(gene_id)
head(p_values)
lowest_p <- data.frame(gene_id, p_values)
lowest_p <- lowest_p[order(p_values),]
head(lowest_p)
tail(lowest_p)
candidate_genes <- lowest_p$gene_id
candidate_genes <- candidate_genes[1:1000]
head(candidate_genes)
dataset <- dataset[ , c(1, candidate_genes)]


rpart_speed <- c()
rpart_precision <- c()
rpart_recall <- c()
rpart_F1score <- c()
rpart_auc <- c()
rpart_complexity <- c()

ctree_speed <- c()
ctree_precision <- c()
ctree_recall <- c()
ctree_F1score <- c()
ctree_auc <- c()
ctree_complexity <- c()

evtree_speed <- c()
evtree_precision <- c()
evtree_recall <- c()
evtree_F1score <- c()
evtree_auc <- c()
evtree_complexity <- c()

tree_speed <- c()
tree_precision <- c()
tree_recall <- c()
tree_F1score <- c()
tree_auc <- c()
tree_complexity <- c()

C5.0_speed <- c()
C5.0_precision <- c()
C5.0_recall <- c()
C5.0_F1score <- c()
C5.0_auc <- c()
C5.0_complexity <- c()

for(k in 1:100) {
  cat("\014")
  cat(k)
  cat("\n")
  
  shuffled <- dataset[sample(nrow(dataset), replace=FALSE), ] # Shuffling dataset
  dim(dataset)
  dim(shuffled)
  rpart_speed_avearge <- c()
  rpart_precision_avearge <- c()
  rpart_recall_avearge <- c()
  rpart_F1score_avearge <- c()
  rpart_auc_avearge <- c()
  rpart_complexity_avearge <- c()
  
  ctree_speed_avearge <- c()
  ctree_precision_avearge <- c()
  ctree_recall_avearge <- c()
  ctree_F1score_avearge <- c()
  ctree_auc_avearge <- c()
  ctree_complexity_avearge <- c()
  
  evtree_speed_avearge <- c()
  evtree_precision_avearge <- c()
  evtree_recall_avearge <- c()
  evtree_F1score_avearge <- c()
  evtree_auc_avearge <- c()
  evtree_complexity_avearge <- c()
  
  tree_speed_avearge <- c()
  tree_precision_avearge <- c()
  tree_recall_avearge <- c()
  tree_F1score_avearge <- c()
  tree_auc_avearge <- c()
  tree_complexity_avearge <- c()
  
  C5.0_speed_avearge <- c()
  C5.0_precision_avearge <- c()
  C5.0_recall_avearge <- c()
  C5.0_F1score_avearge <- c()
  C5.0_auc_avearge <- c()
  C5.0_complexity_avearge <- c()
  
  i=1 
  j=47
  for(x in 1:10) {
    test_set <- shuffled[i:j, ]
    train_set <- shuffled[-c(i:j), ]
    #--------------------------------------------------------------------------
    #rpart
    start_time <- Sys.time()
    tree1 <- rpart(X ~., data = train_set, method = "class")
    end_time <- Sys.time()
    rpart_speed_avearge <- append(rpart_speed_avearge, (as.numeric(difftime(end_time,start_time, units = c("secs")))))
    result <- predict(tree1, test_set)[,2]
    confu_matrix <- table(factor(test_set$X, levels = c(0,1)), factor(result >= 0.5,levels = c(F,T)))
    
    stopifnot((ncol(confu_matrix))&&(nrow(confu_matrix)) == 2)
    
    rpart_precision_avearge <- append(rpart_precision_avearge, confu_matrix[2,2]/(confu_matrix[2,2]+confu_matrix[1,2]))
    rpart_recall_avearge <- append(rpart_recall_avearge, confu_matrix[2,2]/(confu_matrix[2,2]+confu_matrix[2,1]))
    rpart_F1score_avearge <- append(rpart_F1score_avearge, (2*((confu_matrix[2,2]/(confu_matrix[2,2]+confu_matrix[1,2]))*(confu_matrix[2,2]/(confu_matrix[2,2]+confu_matrix[2,1])))/((confu_matrix[2,2]/(confu_matrix[2,2]+confu_matrix[1,2]))+(confu_matrix[2,2]/(confu_matrix[2,2]+confu_matrix[2,1]))))
    )
    
    test_error <- try({roc <- roc(result, test_set$X)}, silent = T)
    if(class(test_error)=="try-error") {
      rpart_auc_avearge <- append(rpart_auc_avearge, 1)
    } else {
      roc <- roc(result, test_set$X)
      rpart_auc_avearge <- append(rpart_auc_avearge, roc$auc[1])
    }
    test_error <- NULL
    rpart_complexity_avearge <- append(rpart_complexity_avearge, c(sum(tree1$frame$var != "<leaf>")+sum(tree1$frame$var == "<leaf>"))) #count number of leaves and branches
    #--------------------------------------------------------------------------
    cat("\014")
    cat(k)
    cat("\n")
    
    #ctree
    start_time <- Sys.time()
    tree2 <- ctree(X ~ ., data = train_set)
    end_time <- Sys.time()
    ctree_speed_avearge <- append(ctree_speed_avearge, (as.numeric(difftime(end_time,start_time, units = c("secs")))))
    
    result <- predict(tree2, test_set)
    confu_matrix <- table(factor(test_set$X, levels = c(0,1)), factor(result >= 0.5,levels = c(F,T)))
    
    stopifnot((ncol(confu_matrix))&&(nrow(confu_matrix)) == 2)
    
    ctree_precision_avearge <- append(ctree_precision_avearge, confu_matrix[2,2]/(confu_matrix[2,2]+confu_matrix[1,2]))
    ctree_recall_avearge <- append(ctree_recall_avearge, confu_matrix[2,2]/(confu_matrix[2,2]+confu_matrix[2,1]))
    ctree_F1score_avearge <- append(ctree_F1score_avearge, (2*((confu_matrix[2,2]/(confu_matrix[2,2]+confu_matrix[1,2]))*(confu_matrix[2,2]/(confu_matrix[2,2]+confu_matrix[2,1])))/((confu_matrix[2,2]/(confu_matrix[2,2]+confu_matrix[1,2]))+(confu_matrix[2,2]/(confu_matrix[2,2]+confu_matrix[2,1]))))
    )
    
    test_error <- try({roc <- roc(result, test_set$X)}, silent = T)
    if(class(test_error)=="try-error") {
      ctree_auc_avearge <- append(ctree_auc_avearge, 1)
    } else {
      roc <- roc(result, test_set$X)
      #roc
      ctree_auc_avearge <- append(ctree_auc_avearge, roc$auc[1])
    }
    test_error <- NULL
    ctree_complexity_avearge <- append(ctree_complexity_avearge, length(tree2)) #count number of leaves and branches
    #plot(tree2, type = "simple")  
    
    #--------------------------------------------------------------------------
    cat("\014")
    cat(k)
    cat("\n")
    #evtree
    start_time <- Sys.time()
    tree3 <- evtree(factor(X) ~., data = train_set) 
    end_time <- Sys.time()
    evtree_speed_avearge <- append(evtree_speed_avearge, (as.numeric(difftime(end_time,start_time, units = c("secs")))))
    
    result <- predict(tree3, test_set)
    confu_matrix <- table(factor(result, levels = c(0,1)), factor(test_set$X, levels = c(0,1)))
    
    
    stopifnot((ncol(confu_matrix))&&(nrow(confu_matrix)) == 2)
    
    evtree_precision_avearge <- append(evtree_precision_avearge, confu_matrix[2,2]/(confu_matrix[2,2]+confu_matrix[2,1]))
    evtree_recall_avearge <- append(evtree_recall_avearge, confu_matrix[2,2]/(confu_matrix[2,2]+confu_matrix[2,1]))
    evtree_F1score_avearge <- append(evtree_F1score_avearge, (2*((confu_matrix[2,2]/(confu_matrix[2,2]+confu_matrix[1,2]))*(confu_matrix[2,2]/(confu_matrix[2,2]+confu_matrix[2,1])))/((confu_matrix[2,2]/(confu_matrix[2,2]+confu_matrix[1,2]))+(confu_matrix[2,2]/(confu_matrix[2,2]+confu_matrix[2,1]))))
    )
    
    test_error <- try({roc <- roc(result, test_set$X)}, silent = T)
    if(class(test_error)=="try-error") {
      evtree_auc_avearge <- append(evtree_auc_avearge, 1)
    } else {
      roc <- roc(result, test_set$X)
      #roc
      evtree_auc_avearge <- append(evtree_auc_avearge, roc$auc[1])
    }
    test_error <- NULL
    
    evtree_complexity_avearge <- append(evtree_complexity_avearge, length(tree3)) #count number of leaves and branches
    #plot(tree3, type = "simple") 
    #--------------------------------------------------------------------------
    cat("\014")
    cat(k)
    cat("\n")
    #tree
    start_time <- Sys.time()
    tree4 <- tree(factor(X)~., data = train_set)
    
    end_time <- Sys.time()
    tree_speed_avearge <- append(tree_speed_avearge, (as.numeric(difftime(end_time,start_time, units = c("secs")))))
    
    result <- predict(tree4, test_set)[,2]
    
    confu_matrix <- table(factor(test_set$X, levels = c(0,1)), factor(result >= 0.5,levels = c(F,T)))
    
    stopifnot((ncol(confu_matrix))&&(nrow(confu_matrix)) == 2)
    
    tree_precision_avearge <- append(tree_precision_avearge, confu_matrix[2,2]/(confu_matrix[2,2]+confu_matrix[1,2]))
    tree_recall_avearge <- append(tree_recall_avearge, confu_matrix[2,2]/(confu_matrix[2,2]+confu_matrix[2,1]))
    tree_F1score_avearge <- append(tree_F1score_avearge, (2*((confu_matrix[2,2]/(confu_matrix[2,2]+confu_matrix[1,2]))*(confu_matrix[2,2]/(confu_matrix[2,2]+confu_matrix[2,1])))/((confu_matrix[2,2]/(confu_matrix[2,2]+confu_matrix[1,2]))+(confu_matrix[2,2]/(confu_matrix[2,2]+confu_matrix[2,1]))))
    )
    
    test_error <- try({roc <- roc(result, test_set$X)}, silent = T)
    if(class(test_error)=="try-error") {
      tree_auc_avearge <- append(tree_auc_avearge, 1)
    } else {
      roc <- roc(result, test_set$X)
      #roc
      tree_auc_avearge <- append(tree_auc_avearge, roc$auc[1])
    }
    test_error <- NULL
    
    tree_complexity_avearge <- append(tree_complexity_avearge,c(sum(tree4$frame$var != "<leaf>")+sum(tree4$frame$var == "<leaf>"))) #count number of leaves and branches
    #plot(tree4)
    
    #--------------------------------------------------------------------------
    cat("\014")
    cat(k)
    cat("\n")
    #C5.0
    start_time <- Sys.time()
    tree5 <- C5.0(factor(X)~., data = train_set)
    
    end_time <- Sys.time()
    C5.0_speed_avearge <- append(C5.0_speed_avearge, (as.numeric(difftime(end_time,start_time, units = c("secs")))))
    result <- predict(tree5, test_set)
    
    confu_matrix <- table(factor(result, levels = c(0,1)), factor(test_set$X, levels = c(0,1)))
    
    stopifnot((ncol(confu_matrix))&&(nrow(confu_matrix)) == 2)
    
    C5.0_precision_avearge <- append(C5.0_precision_avearge, confu_matrix[2,2]/(confu_matrix[2,2]+confu_matrix[2,1]))
    C5.0_recall_avearge <- append(C5.0_recall_avearge, confu_matrix[2,2]/(confu_matrix[2,2]+confu_matrix[2,1]))
    C5.0_F1score_avearge <- append(C5.0_F1score_avearge, (2*((confu_matrix[2,2]/(confu_matrix[2,2]+confu_matrix[1,2]))*(confu_matrix[2,2]/(confu_matrix[2,2]+confu_matrix[2,1])))/((confu_matrix[2,2]/(confu_matrix[2,2]+confu_matrix[1,2]))+(confu_matrix[2,2]/(confu_matrix[2,2]+confu_matrix[2,1]))))
    )
    
    test_error <- try({roc <- roc(result, test_set$X)}, silent = T)
    if(class(test_error)=="try-error") {
      C5.0_auc_avearge <- append(C5.0_auc_avearge, 1)
    } else {
      roc <- roc(result, test_set$X)
      #roc
      C5.0_auc_avearge <- append(C5.0_auc_avearge, roc$auc[1])
    }
    test_error <- NULL
    
    C5.0_complexity_avearge <- append(C5.0_complexity_avearge,length(as.party(tree5)))
    #plot(tree5) 
    
    #--------------------------------------------------------------------------
    
    i=j+1
    if(i == 424) {
      j=j+47
    } else {
      j=j+47
    }
  }
  #rpart
  rpart_speed <- append(rpart_speed, mean(rpart_speed_avearge))
  rpart_precision <- append(rpart_precision, mean(rpart_precision_avearge))
  rpart_recall <- append(rpart_recall, mean(rpart_recall_avearge))
  rpart_F1score <- append(rpart_F1score, mean(rpart_F1score_avearge))
  rpart_auc <- append(rpart_auc, mean(rpart_auc_avearge))
  rpart_complexity <- append(rpart_complexity, mean(rpart_complexity_avearge))
  
  #ctree
  ctree_speed <- append(ctree_speed, mean(ctree_speed_avearge))
  ctree_precision <- append(ctree_precision, mean(ctree_precision_avearge))
  ctree_recall <- append(ctree_recall, mean(ctree_recall_avearge))
  ctree_F1score <- append(ctree_F1score, mean(ctree_F1score_avearge))
  ctree_auc <- append(ctree_auc, mean(ctree_auc_avearge))
  ctree_complexity <- append(ctree_complexity, mean(ctree_complexity_avearge))
  
  #evtree
  evtree_speed <- append(evtree_speed, mean(evtree_speed_avearge))
  evtree_precision <- append(evtree_precision, mean(evtree_precision_avearge))
  evtree_recall <- append(evtree_recall, mean(evtree_recall_avearge))
  evtree_F1score <- append(evtree_F1score, mean(evtree_F1score_avearge))
  evtree_auc <- append(evtree_auc, mean(evtree_auc_avearge))
  evtree_complexity <- append(evtree_complexity, mean(evtree_complexity_avearge))
  
  #tree
  tree_speed <- append(tree_speed, mean(tree_speed_avearge))
  tree_precision <- append(tree_precision, mean(tree_precision_avearge))
  tree_recall <- append(tree_recall, mean(tree_recall_avearge))
  tree_F1score <- append(tree_F1score, mean(tree_F1score_avearge))
  tree_auc <- append(tree_auc, mean(tree_auc_avearge))
  tree_complexity <- append(tree_complexity, mean(tree_complexity_avearge))
  
  #C5.0
  C5.0_speed <- append(C5.0_speed, mean(C5.0_speed_avearge))
  C5.0_precision <- append(C5.0_precision, mean(C5.0_precision_avearge))
  C5.0_recall <- append(C5.0_recall, mean(C5.0_recall_avearge))
  C5.0_F1score <- append(C5.0_F1score, mean(C5.0_F1score_avearge))
  C5.0_auc <- append(C5.0_auc, mean(C5.0_auc_avearge))
  C5.0_complexity <- append(C5.0_complexity, mean(C5.0_complexity_avearge))
  
}


speed_table <- data.frame("rpart" = rpart_speed, "ctree" = ctree_speed, "evtree" = evtree_speed, "tree" = tree_speed, "C5.0" = C5.0_speed)
speed_table
write.csv(speed_table, "~/Documents/R/GSE102299-smartseq2/speed_table.csv")

precision_table <- data.frame("rpart" = rpart_precision, "ctree" = ctree_precision, "evtree" = evtree_precision, "tree" = tree_precision, "C5.0" = C5.0_precision)
precision_table
write.csv(precision_table, "~/Documents/R/GSE102299-smartseq2/precision_table.csv")

recall_table <- data.frame("rpart" = rpart_recall, "ctree" = ctree_recall, "evtree" = evtree_recall, "tree" = tree_recall, "C5.0" = C5.0_recall)
recall_table
write.csv(recall_table, "~/Documents/R/GSE102299-smartseq2/recall_table.csv")

F1score_table <- data.frame("rpart" = rpart_F1score, "ctree" = ctree_F1score, "evtree" = evtree_F1score, "tree" = tree_F1score, "C5.0" = C5.0_F1score)
F1score_table
write.csv(F1score_table, "~/Documents/R/GSE102299-smartseq2/F1score_table.csv")

auc_table <- data.frame("rpart" = rpart_auc, "ctree" = ctree_auc, "evtree" = evtree_auc, "tree" = tree_auc, "C5.0" = C5.0_auc)
auc_table
write.csv(auc_table, "~/Documents/R/GSE102299-smartseq2/auc_table.csv")

complexity_table <- data.frame("rpart" = rpart_complexity, "ctree" = ctree_complexity, "evtree" = evtree_complexity, "tree" = tree_complexity, "C5.0" = C5.0_complexity)
complexity_table
write.csv(complexity_table, "~/Documents/R/GSE102299-smartseq2/complexity_table.csv")


pdf("~/Documents/R/GSE102299-smartseq2/Speed.pdf")
boxplot(speed_table, main = "Speed" , xlab = "Methods", ylab = "Speed", col=c("DarkSalmon", "CadetBlue", "Darkkhaki", "DarkSeaGreen", "Thistle"))
dev.off()

pdf("~/Documents/R/GSE102299-smartseq2/Precision.pdf")
boxplot(precision_table, main = "Precision" , xlab = "Methods", ylab = "Precision", col=c("DarkSalmon", "CadetBlue", "Darkkhaki", "DarkSeaGreen", "Thistle"))
dev.off()

pdf("~/Documents/R/GSE102299-smartseq2/Recall.pdf")
boxplot(recall_table, main = "Recall" , xlab = "Methods", ylab = "Recall", col=c("DarkSalmon", "CadetBlue", "Darkkhaki", "DarkSeaGreen", "Thistle"))
dev.off()

pdf("~/Documents/R/GSE102299-smartseq2/F1score.pdf")
boxplot(F1score_table, main = "F1score" , xlab = "Methods", ylab = "F1score", col=c("DarkSalmon", "CadetBlue", "Darkkhaki", "DarkSeaGreen", "Thistle"))
dev.off()

pdf("~/Documents/R/GSE102299-smartseq2/AUC.pdf")
boxplot(auc_table, main = "Area Under the Curve" , xlab = "Methods", ylab = "AUC(%)", col=c("DarkSalmon", "CadetBlue", "Darkkhaki", "DarkSeaGreen", "Thistle"))
dev.off()

pdf("~/Documents/R/GSE102299-smartseq2/Complexity.pdf")
boxplot(complexity_table, main = "Complexity" , xlab = "Methods", ylab = "Complexity", col=c("DarkSalmon", "CadetBlue", "Darkkhaki", "DarkSeaGreen", "Thistle"))
dev.off()

