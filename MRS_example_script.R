############################### Example script to showcase methylation risk scores (MRS)
############################### Michael Thompson, Brian Hill, Noah Zaitlen, Eran Halperin
############################### Los Angeles, August 2nd, 2022
# See academic license
## Load relevant packages
library(data.table)
library(dplyr)
require(pROC)
## First, download an example dataset
setwd("~/.")
system("wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE20nnn/GSE20067/matrix/GSE20067_series_matrix.txt.gz")
# Unzip it
system("gunzip GSE20067_series_matrix.txt.gz")
# Get the case control status of individuals
casectrl=system("grep case/control GSE20067_series_matrix.txt", intern=T)
casectrl=unlist(strsplit(casectrl, split='\t'))[-1]
casectrl=sapply(strsplit(casectrl, split=': '), function(x) x[2])
casectrl=sapply(strsplit(casectrl, split='\"'), function(x) x[1])
casectrl.bin=abs(as.numeric(as.factor(casectrl))-2)
# Get the sample names
inds=system("grep ID_REF GSE20067_series_matrix.txt",intern=T)
inds=unlist(strsplit(inds, split='\t'))[-1]
inds=sapply(strsplit(inds, split='\"'), function(x) x[2])
casedf=data.frame(ID=inds, Case=casectrl.bin)
rownames(casedf)=inds
df=fread("GSE20067_series_matrix.txt", skip=72,data.table=F)
rownames(df)=df$ID_REF
df=as.matrix(df[,-1])
storage.mode(df)<-"numeric"
sumna=sapply(1:nrow(df), function(i) sum(is.na(df[i,])))
## Perform some standard preprocessing of the data
# remove sites with < 10% of individuals
df=df[sumna <= 19,]
sumna=sapply(1:ncol(df), function(i) sum(is.na(df[,i])))
df=df[,sumna <= nrow(df)/10]
# mean impute the remaining sites
for(i in 1:ncol(df)){
  miss=which(is.na(df[,i]))
  if(length(miss) > 0){
    df[miss,i]=mean(df[,i], na.rm = T)
  }
}
casedf=casedf[intersect(rownames(casedf), colnames(df)),]

## Now download a weight you wish to use
system("wget https://raw.githubusercontent.com/cozygene/EHR_MRS_UCLA/main/weights450k/labs/CREATININE.txt")
weight=fread("CREATININE.txt")
wgt.matrix=matrix(data=weight$beta,ncol=1,dimnames=list(weight$cpg, "beta"))
## Use the weights to impute an MRS for creatinine into this downloaded dataset
# First make sure the CpGs (methylation sites, or features) overlap
# One way of doing this is by taking the intersection of features
wgt.matrix=wgt.matrix[intersect(rownames(wgt.matrix), rownames(df)),,drop=F]
# To perform imputation, take a weighted combination of CpG values per individual
# weighted by the weight specified in the MRS weight file (CREATININE.txt)
# store this quantity as "pred"
pred=t(df[rownames(wgt.matrix),rownames(casedf)]) %*% wgt.matrix
casedf$pred=c(pred)
# check if this imputed MRS for creatinine is significantly
# associated with the case control status
summary(glm(Case ~ pred, family="binomial", data=casedf))
# evaluate how predictive this imputed creatinine value is
# of case contorl status using AUC
pROC::auc(predictor=casedf$pred, response=casedf$Case)
pROC::ci.auc(predictor=casedf$pred, response=casedf$Case, boot.n=1000,
             method="b")
