require(beanplot)
require(tidyr)
#install.packages(c("FactoMineR", "factoextra"))
#install.packages("dplyr", dependencies = TRUE)
library(FactoMineR)
library(factoextra)
library(yarrr)
library(tidyr)
library(tidyverse)
library(table1)
library(sjlabelled)
library(dplyr)
library(utils)
library(glmnet)
#### Make Luksza Input Files ####
setwd("~/3.0 Hasting Research/Novel Model/Cohen Melanoma/Penalized_Regression_Data")
samples <- matrix(c("3466"))
for (i in samples){
  # Read in data
  valid_data <- as.data.frame(read.table(paste("Validated_",i,".txt", sep=""), header=TRUE))
  data_MT <- as.data.frame(read.table(paste(i, "_netctl.out", sep=""), header=TRUE))
  data_WT <- as.data.frame(read.table(paste(i, ".WTmatch_netCTL.out", sep=""), header=TRUE))
  # Match data
  valid_data <- valid_data[!duplicated(valid_data[,3]),]
  data_MT <- data_MT[data_MT$Peptide %in% valid_data$Peptide,]
  data_MT <- data_MT[!duplicated(data_MT[,1]),]
  data_WT <- data_WT[data_WT$Peptide %in% data_MT$Peptide,]
  data_WT <- data_WT[!duplicated(data_WT[,1]),]
  data_WT <- data_WT[match(data_MT$Peptide, data_WT$Peptide),]
  valid_data <- valid_data[match(data_MT$Peptide, valid_data$Peptide),]
  # Isolate data
  MT_peptide <- data_MT[,1]
  MT_MHC <- data_MT[,3]
  MT_MHC <- 50^MT_MHC
  WT_peptide <- data_WT[,2]
  WT_MHC <- data_WT[,3]
  WT_MHC <- 50^WT_MHC
  genename <- valid_data[,2]
  l <- nrow(valid_data)
  id <- matrix(c("3903"),nrow=l, ncol=1)
  mut <- matrix(c("MUT"),nrow=l, ncol=1)
  forluksza_total <- cbind(genename, id, WT_peptide, MT_peptide, mut, WT_MHC, MT_MHC)
  forluksza_total <- as.data.frame(forluksza_total)
  # Print data to file
  sink(file=paste(i,"_lukszaprogram.in", sep=""))
  print(forluksza_total)
  sink(NULL)
}
#### Read in Data ####
samples <- matrix(c("3466", "3702", "3703h11.1", "3703h2.1", "3713", 
                    "3879h1.1", "3879h2.1", "3919", "3926h1.1",
                    "3926h1.1"))
sink(file="Overall_data_for_regression.txt")
for (i in samples){
  valid_data <- as.data.frame(read.table(paste("Validated_",i,".txt", sep=""), header=TRUE))
  data_MT <- as.data.frame(read.table(paste(i, "_netctl.out", sep=""), header=TRUE))
  data_WT <- as.data.frame(read.table(paste(i, ".WTmatch_netCTL.out", sep=""), header=TRUE))
  sequence_alignment <- as.data.frame(read.table(paste(i, ".epitopes.annotated.tsv", sep=""), header=TRUE, row.names = NULL))
  TCR_binding <- as.data.frame(read.table(paste("neoantigen_fitness_", i, ".txt", sep=""), header=TRUE))
  quants <- as.data.frame(read.table(paste(i,"_quant.sf", sep=""), header=TRUE))
  
  #### Isolate Overlap and Match Order ####
  valid_data <- valid_data[!duplicated(valid_data[,3]),]
  data_MT <- data_MT[data_MT$Peptide %in% valid_data$Peptide,]
  data_MT <- data_MT[!duplicated(data_MT[,1]),]
  data_WT <- data_WT[data_WT$Peptide %in% data_MT$Peptide,]
  data_WT <- data_WT[!duplicated(data_WT[,1]),]
  data_WT <- data_WT[match(data_MT$Peptide, data_WT$Peptide),]
  valid_data <- valid_data[match(data_MT$Peptide, valid_data$Peptide),]
  sequence_alignment <- sequence_alignment[match(data_MT$Peptide, sequence_alignment$row.names),]
  TCR_binding <- TCR_binding[match(data_MT$Peptide, TCR_binding$MutantPeptide),]
  
  #### Isolating RNA Expression ####
  quants[,1]=gsub(".18","",quants[,1], fixed=TRUE)
  quants[,1]=gsub(".17","",quants[,1], fixed=TRUE)
  quants[,1]=gsub(".16","",quants[,1], fixed=TRUE)
  quants[,1]=gsub(".15","",quants[,1], fixed=TRUE)
  quants[,1]=gsub(".14","",quants[,1], fixed=TRUE)
  quants[,1]=gsub(".13","",quants[,1], fixed=TRUE)
  quants[,1]=gsub(".12","",quants[,1], fixed=TRUE)
  quants[,1]=gsub(".11","",quants[,1], fixed=TRUE)
  quants[,1]=gsub(".10","",quants[,1], fixed=TRUE)
  quants[,1]=gsub(".9","",quants[,1], fixed=TRUE)
  quants[,1]=gsub(".8","",quants[,1], fixed=TRUE)
  quants[,1]=gsub(".7","",quants[,1], fixed=TRUE)
  quants[,1]=gsub(".6","",quants[,1], fixed=TRUE)
  quants[,1]=gsub(".5","",quants[,1], fixed=TRUE)
  quants[,1]=gsub(".4","",quants[,1], fixed=TRUE)
  quants[,1]=gsub(".3","",quants[,1], fixed=TRUE)
  quants[,1]=gsub(".2","",quants[,1], fixed=TRUE)
  quants[,1]=gsub(".1","",quants[,1], fixed=TRUE)
  Gencode <- as.data.frame(read.table("ensemblTranscripts_GTF_new.txt"))
  Gencode <- Gencode[Gencode$V7 %in% quants$Name,]
  quants <- quants[match(Gencode$V7, quants$Name),]
  total <- cbind(Gencode, quants)
  total <- total[total$V8 %in% TCR_binding$Mutation,]
  total_truncated <- total[,c(8,12)]
  total_truncated <- aggregate(TPM ~ V8, data=total_truncated, FUN=sum)
  total_truncated <- total_truncated[match(TCR_binding$Mutation, total_truncated$V8),]
  #### Isolate data for regression ####
  for_regression_MT <- data_MT[,1:4]
  WT_MHC <- data_WT[,3]
  Seqalign <- sequence_alignment[,4]
  TCR <- TCR_binding[,7]
  expression <- total_truncated[,2]
  for_regression_outcome <- valid_data[,4]
  for_regression_outcome <- as.data.frame(for_regression_outcome)
  sink(file="For_regression_outocme.txt",append=TRUE)
  print(for_regression_outcome)
  sink(NULL)
  for_regression <- cbind(for_regression_MT, WT_MHC, Seqalign, TCR, expression)
  print(for_regression)
}
sink(NULL)

#### Penalized Regression ####
# Good time to clean environment
setwd("~/3.0 Hasting Research/Novel Model/Cohen Melanoma/Penalized_Regression_Data")
data <- as.data.frame(read.table("Overall_data_for_regression.txt", row.names = NULL))
outcomes <- as.data.frame(read.table("For_regression_outocme.txt", row.names =NULL, header =TRUE))
data <- data[,-1]
data[is.na(data)] <- 0 
#data <- data[,2:7]
outcomes <- outcomes[,-1]
outcomes <- as.data.frame(outcomes)
ids <- matrix(c("no"), nrow=370)
colnames(ids) <- "immunogenic"
y <- ids
for (i in 1:370){
  if (outcomes[i,]>0){
    y[i,]="yes"
  }
}
fraction_yes <- rep(1-sum(y=="yes")/nrow(y), sum(y =="yes"))
fraction_no <- rep(1-sum(y=="no")/nrow(y), sum(y=="no"))
weights <- numeric(nrow(y))
weights[y=="yes"] <- fraction_yes
weights[y=="no"] <- fraction_no
x <- model.matrix( ~ MHC  + TAP + Cle + WT_MHC + Seqalign + TCR + expression - 1, data) 
glm_fit <- glmnet(x=x, y=y, family="binomial", weights=weights)
plot(glm_fit, xvar="dev", label=TRUE, cex.lab=1.5, cex.axis=1.5, lwd =2, cex=1.5, cex.main=1.5)
summary(glm_fit)
cv.fit <- cv.glmnet(x=x, y=y, family="binomial", type.measure = "class", weights=weights)
summary(cv.fit)
plot(cv.fit, cex.lab=1.5, cex.axis=1.5, cex=1.5, cex.main=1.5)
coef(cv.fit,s="lambda.min")


#### Fusion format output ####
setwd("~/3.0 Hasting Research/Novel Model/Cohen Melanoma/Fusions/")
samples <- matrix(c("fusion"))
i="fusion"
sink(file="Fusion_Overall_data.txt")
data_MT <- as.data.frame(read.table(paste(i, "_netctl.out", sep=""), header=TRUE))
data_WT <- as.data.frame(read.table(paste(i, ".WTmatch_netCTL.out", sep=""), header=TRUE))
sequence_alignment <- as.data.frame(read.table(paste(i, ".epitopes.annotated.tsv", sep=""), header=TRUE, row.names = NULL))
TCR_binding <- as.data.frame(read.table(paste("neoantigen_fitness_", i, ".txt", sep=""), header=TRUE))
#### Isolate Overlap and Match Order ####
data_MT <- data_MT[!duplicated(data_MT[,1]),]
data_WT <- data_WT[data_WT$Peptide %in% data_MT$Peptide,]
data_WT <- data_WT[!duplicated(data_WT[,1]),]
data_WT <- data_WT[match(data_MT$Peptide, data_WT$Peptide),]
sequence_alignment <- sequence_alignment[match(data_MT$Peptide, sequence_alignment$row.names),]
TCR_binding <- TCR_binding[match(data_MT$Peptide, TCR_binding$MutantPeptide),]
#### Isolate data for regression ####
for_regression_MT <- data_MT[,1:4]
WT_MHC <- data_WT[,3]
Seqalign <- sequence_alignment[,4]
TCR <- TCR_binding[,7]
for_regression <- cbind(for_regression_MT, WT_MHC, Seqalign, TCR)
print(for_regression)
sink(NULL)
