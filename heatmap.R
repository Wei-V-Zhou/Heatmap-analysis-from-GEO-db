

##################
# Load libraries #
##################
# clear objectives and garbage collection
rm(list=ls())
gc()
options(stringsAsFactors=FALSE)
# load packages
library(plyr)
library(limma)
library(Biobase)
library(ggplot2)
library(stringr)
library(reshape2)
library(GEOquery)
library(pheatmap)
library(ggfortify)
library(hugene10sttranscriptcluster.db)

# load data
if(!file.exists("GSE2603_gSet.Rdata")){
  gSet<-getGEO("GSE2603",destdir=".",GSEMatrix=TRUE,
               AnnotGPL=FALSE,getGPL=FALSE)
  if (length(gSet)>1) idx<-grep("GPL96",attr(gSet,"names")) else idx<-1
  gSet<-gSet[[idx]]
  save(gSet,file="GSE2603_gSet.Rdata")
}
load("GSE2603_gSet.Rdata")
# log2 transform
ex<-exprs(gSet)
qx<-as.numeric(quantile(ex,c(0.,0.25,0.5,0.75,0.99,1.0),na.rm=T))
LogC<-(qx[5]>100)||
  (qx[6]-qx[1]>50 && qx[2]>0)||
  (qx[2]>0 && qx[2]<1 && qx[4]>1 && qx[4]<2)
if(LogC){ex[which(ex<=0)]<-NaN
exprs(gSet)<-log2(ex)}
# load the expression data
exprSet<-exprs(gSet) 
# load the sample name
samples<-sampleNames(gSet)
# load the annotation
pdata <- pData(gSet)

# load group name
group_list<-as.character(pdata[,8])
{
  b_expr=exprSet[,grep("Breast cancer cell line",group_list)]
  ba_expr=as.matrix(exprSet[,grep("breast cancer cell line",group_list)])
  s_expr=exprSet[,grep("Subpopulation",group_list)]
  sa_expr=as.matrix(exprSet[,grep("subpopulation",group_list)])
  p_expr=exprSet[,grep("Primary Breast Cancer",group_list)]
}
{
  exprSet=cbind(b_expr,ba_expr,s_expr,sa_expr,p_expr)
  group_list=c(rep("Breast_cancer_cell_line",ncol(ba_expr)+1),
               rep("Subpopulation",ncol(sa_expr)),
               rep("Subpopulation",ncol(s_expr)),
               rep("Breast_cancer_cell_line",ncol(b_expr)-1),
               rep("Primary_Breast_Cancer",ncol(p_expr)))
}

# load the GPL data
if(T){
  gpl <- getGEO("GPL96", destdir = ".")
  probe2gene<-Table(gpl)[,c(1,11)]
}
# filter the probe with null genes
ids=probe2gene[probe2gene[,2]!='',]
# filter the probe with many genes
a<-strsplit(as.character(ids[,2])," /// ")
tmp<-mapply(cbind,ids[,1],a)
df<-ldply(tmp,data.frame)
probe2gene=df[,2:3]
save(probe2gene,file="probe2gene.Rdata")

# remove the probe without gene annotation
if(T){
  exprSet=exprSet[rownames(exprSet) %in% probe2gene[,1],]
  probe2gene=probe2gene[match(rownames(exprSet),probe2gene[,1]),]
}
dim(exprSet)
dim(probe2gene)
tail(sort(table(probe2gene[,2])),n=12L)

# get the maximum expression of the same gene
{
  MAX=by(exprSet,probe2gene[,2],
         function(x) rownames(x)[which.max(rowMeans(x))])
  MAX=as.character(MAX)
  exprSet=exprSet[rownames(exprSet) %in% MAX,]
  rownames(exprSet)=probe2gene[match(rownames(exprSet),probe2gene[,1]),2]
}
dim(exprSet)
exprSet[1:5,1:5]
save(exprSet,group_list,file="final_exprSet.Rdata")

# test the gene expression
exprSet_L<-melt(exprSet)
colnames(exprSet_L)<-c('probe','sample','value')
dat1=exprSet_L[exprSet_L$probe=="GAPDH",]
dat2=exprSet_L[exprSet_L$probe=="ACTB",]
p<-ggplot(dat2,aes(x=sample,y=value))+geom_boxplot()
print(p)

# hcluster plot
{
  colnames(exprSet)=paste(group_list,1:ncol(exprSet))
  nodepar<-list(lab.cex=0.6,pch=c(NA,19),
                cex=0.7,col="blue")
  hc<-hclust(dist(t(exprSet)))
  par(mar=c(5,5,5,10))
  plot(as.dendrogram(hc),nodePar=nodepar,horiz=F)
}



# PCA plot
df<-as.data.frame(t(exprSet))
df$group<-group_list
autoplot(prcomp(df[,1:(ncol(df)-1)]),data=df,colour="group")

# heatmap plot
rm(list=ls())
load("final_exprSet.Rdata")
cg=names(tail(sort(apply(exprSet,1,sd)),1000))
pheatmap(exprSet[cg,],show_colnames=F,show_rownames=F) 
{
  n=t(scale(t(exprSet[cg,])))
  n[n>2]=2 
  n[n< -2]= -2
  n[1:4,1:4]
}
pheatmap(n,show_colnames =F,show_rownames=F)
ac=data.frame(g=group_list)
colnames(n)=rownames(ac)
pheatmap(n,show_colnames=F,show_rownames=F,annotation_col=ac)

# differential expression genes
{
  design<-model.matrix(~0+factor(group_list))
  colnames(design)<-levels(factor(group_list))
  rownames(design)<-colnames(exprSet)
}

contrast.matrix<-makeContrasts(paste0(unique(group_list),collapse="-"),levels=design)
load("./probe2gene.Rdata")
{
  fit<-lmFit(exprSet,design)
  fit2<-contrasts.fit(fit,contrast.matrix)
  fit2<-eBayes(fit2,0.01)
  options(digits=4)
}
tempOutput=topTable(fit2,coef=1,n=Inf)
nrDEG=na.omit(tempOutput)
save(nrDEG,file="nrDEG.Rdata")

palette(c("#dfeaf4","#f4dfdf","#f2cb98", "#AABBCC"))
dev.new(width=4+dim(get)[[2]]/5, height=6)
par(mar=c(2+round(max(nchar(sampleNames(gSet)))/2),4,2,1))
title <- paste ("GSE2603", '/', annotation(gSet), " selected samples", sep ='')
boxplot(ex, boxwex=0.6, notch=T, main=title, outline=FALSE, las=2)
legend("topleft", labels, fill=palette(), bty="n")

