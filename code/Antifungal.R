options(stringsAsFactors = FALSE)
library(PharmacoGx)
library(VennDiagram)
library(RColorBrewer)
library(R.utils)
library(downloader)
library(gdata)
library(AnnotationDbi)
library(hgu133a.db)
library(gridExtra)
library(biomaRt)
library(rJava)
library(xlsx)
library(ggplot2)
library(RankAggreg)
nbcore <- 4

drug.perturbation<-"CMAP_signatures.RData"
if(!file.exists(drug.perturbation)){
  drug.perturbation<-PharmacoGx::downloadPertSig("CMAP")
}else{
  load(drug.perturbation)
}

ket<-drug.perturbation[,'ketoconazole',]
ket<-ket[order(ket[,1],decreasing = T),]
ket_up<-ket[1:10,c('fdr','pvalue')]
ket<-ket[order(ket[,1],decreasing = FALSE),]
ket_down<-ket[1:10,c('fdr','pvalue')]
ketoconazole<-c(rep(1,times=10),rep(-1,times=10))
names(ketoconazole)<-c(rownames(ket_up),rownames(ket_down))
ket_myfn<-"CMAP_genes_ketoconazole_connectivity.RData"
if(!file.exists(ket_myfn)){
  message("Be aware that computing sensitivity will take some time...")
  cl<-parallel::makeCluster(nbcore)
  ket_re<-parApply(drug.perturbation[ , ,c("tstat","pvalue")],2,function(x,ketoconazole){return(PharmacoGx::connectivityScore(x=x,y=ketoconazole,method = "gsea",nperm = 100))},cl=cl,ketoconazole=ketoconazole)
  stopCluster(cl)
  rownames(ket_re)<-c("Connectivity","P_Value")
  ket_re<-t(ket_re)
  save(ket_re,file=ket_myfn)
}else{
  load(ket_myfn)
}
ket_re1<-ket_re[order(ket_re[,1],decreasing = T),]
ket_re1<-ket_re1[1:600,c("Connectivity","P_Value")]
ket_top<-ket_re1[1:10,c("Connectivity","P_Value")]
ket_top<-round(ket_top,digits = 3)

pdf('ketoconazole.pdf')
drug<-rownames(ket_top)
score<-ket_top[,1]
df<-data.frame(x = drug,y = score)
ggplot(df,aes(x=reorder(drug,score),y=score))+geom_bar(stat = 'identity',width = 0.5)+coord_flip()
dev.off()
write.csv(ket_top,file = 'ketoconazole.csv')
# jpeg('ketoconazole.jpg')
# drug<-rownames(ket_top)
# score<-ket_top[,1]
# df<-data.frame(x = drug,y = score)
# ggplot(df,aes(x=reorder(drug,score),y=score))+geom_bar(stat = 'identity',width = 0.5)+coord_flip()
# dev.off()
###############################################################

mic<-drug.perturbation[,'miconazole',]
mic<-mic[order(mic[,1],decreasing = T),]
mic_up<-mic[1:10,c('fdr','pvalue')]
mic<-mic[order(mic[,1],decreasing = FALSE),]
mic_down<-mic[1:10,c('fdr','pvalue')]
miconazole<-c(rep(1,times=10),rep(-1,times=10))
names(miconazole)<-c(rownames(mic_up),rownames(mic_down))
mic_myfn<-"CMAP_genes_miconazole_connectivity.RData"
if(!file.exists(mic_myfn)){
  message("Be aware that computing sensitivity will take some time...")
  cl<-parallel::makeCluster(nbcore)
  mic_re<-parApply(drug.perturbation[ , ,c("tstat","pvalue")],2,function(x,miconazole){return(PharmacoGx::connectivityScore(x=x,y=miconazole,method = "gsea",nperm = 100))},cl=cl,miconazole=miconazole)
  stopCluster(cl)
  rownames(mic_re)<-c("Connectivity","P_Value")
  mic_re<-t(mic_re)
  save(mic_re,file=mic_myfn)
}else{
  load(mic_myfn)
}
mic_re1<-mic_re[order(mic_re[,1],decreasing = T),]
mic_re1<-mic_re1[1:600,c("Connectivity","P_Value")]
mic_top<-mic_re1[1:10,c("Connectivity","P_Value")]
mic_top<-round(mic_top,digits = 3)

pdf('miconazole.pdf')
drug<-rownames(mic_top)
score<-mic_top[,1]
df<-data.frame(x = drug,y = score)
ggplot(df,aes(x=reorder(drug,score),y=score))+geom_bar(stat = 'identity',width = 0.5)+coord_flip()
dev.off()
write.csv(mic_top,file = 'miconazole.csv')
##########################################################
amB<-drug.perturbation[,'amphotericin B',]
amB<-amB[order(amB[,1],decreasing = T),]
amB_up<-amB[1:10,c('fdr','pvalue')]
amB<-amB[order(amB[,1],decreasing = FALSE),]
amB_down<-amB[1:10,c('fdr','pvalue')]
amphotericin<-c(rep(1,times=10),rep(-1,times=10))
names(amphotericin)<-c(rownames(amB_up),rownames(amB_down))
myfn3<-"CMAP_amB_connectivity.RData"
if(!file.exists(myfn3)){
  message("Be aware that computing sensitivity will take some time...")
  cl<-parallel::makeCluster(nbcore)
  amB_re<-parApply(drug.perturbation[ , ,c("tstat","pvalue")],2,function(x,amphotericin){return(PharmacoGx::connectivityScore(x=x,y=amphotericin,method = "gsea",nperm = 100))},cl=cl,amphotericin=amphotericin)
  stopCluster(cl)
  rownames(amB_re)<-c("Connectivity","P_Value")
  amB_re<-t(amB_re)
  save(amB_re,file=myfn3)
}else{
  load(myfn3)
}
amB_re1<-amB_re[order(amB_re[,1],decreasing = T),]
amB_re1<-amB_re1[1:600,c("Connectivity","P_Value")]
amB_top<-amB_re1[1:10,c("Connectivity","P_Value")]
amB_top<-round(amB_top,digits = 3)

pdf('amphotericin_B.pdf')
drug<-rownames(amB_top)
score<-amB_top[,1]
df<-data.frame(x = drug,y = score)
ggplot(df,aes(x=reorder(drug,score),y=score))+geom_bar(stat = 'identity',width = 0.5)+coord_flip()
dev.off()


#######################################################
nys<-drug.perturbation[,'nystatin',]
nys<-nys[order(nys[,1],decreasing = T),]
nys_up<-nys[1:10,c('fdr','pvalue')]
nys<-nys[order(nys[,1],decreasing = FALSE),]
nys_down<-nys[1:10,c('fdr','pvalue')]
nystatin<-c(rep(1,times=10),rep(-1,times=10))
names(nystatin)<-c(rownames(nys_up),rownames(nys_down))
myfn4<-"CMAP_nystatin_connectivity.RData"
if(!file.exists(myfn4)){
  message("Be aware that computing sensitivity will take some time...")
  cl<-parallel::makeCluster(nbcore)
  nys_re<-parApply(drug.perturbation[ , ,c("tstat","pvalue")],2,function(x,nystatin){return(PharmacoGx::connectivityScore(x=x,y=nystatin,method = "gsea",nperm = 100))},cl=cl,nystatin=nystatin)
  stopCluster(cl)
  rownames(nys_re)<-c("Connectivity","P_Value")
  nys_re<-t(nys_re)
  save(nys_re,file=myfn4)
}else{
  load(myfn4)
}
nys_re1<-nys_re[order(nys_re[,1],decreasing = T),]
nys_re1<-nys_re1[1:600,c("Connectivity","P_Value")]
nys_top<-nys_re1[1:10,c("Connectivity","P_Value")]
nys_top<-round(nys_top,digits = 3)

pdf('nystatin.pdf')
drug<-rownames(nys_top)
score<-nys_top[,1]
df<-data.frame(x = drug,y = score)
ggplot(df,aes(x=reorder(drug,score),y=score))+geom_bar(stat = 'identity',width = 0.5)+coord_flip()
dev.off()

###############################################################
cop<-drug.perturbation[,'copper sulfate',]
cop<-cop[order(cop[,1],decreasing = T),]
cop_up<-cop[1:10,c('fdr','pvalue')]
cop<-cop[order(cop[,1],decreasing = FALSE),]
cop_down<-cop[1:10,c('fdr','pvalue')]
copper<-c(rep(1,times=10),rep(-1,times=10))
names(copper)<-c(rownames(cop_up),rownames(cop_down))
myfn5<-"CMAP_cop_connectivity.RData"
if(!file.exists(myfn5)){
  message("Be aware that computing sensitivity will take some time...")
  cl<-parallel::makeCluster(nbcore)
  cop_re<-parApply(drug.perturbation[ , ,c("tstat","pvalue")],2,function(x,copper){return(PharmacoGx::connectivityScore(x=x,y=copper,method = "gsea",nperm = 100))},cl=cl,copper=copper)
  stopCluster(cl)
  rownames(cop_re)<-c("Connectivity","P_Value")
  cop_re<-t(cop_re)
  save(cop_re,file=myfn5)
}else{
  load(myfn5)
}
cop_re<-cop_re[order(cop_re[,1],decreasing = T),]
cop_top<-cop_re[1:10,c("Connectivity","P_Value")]
cop_top<-round(cop_top,digits = 6)
cop_re<-cop_re[order(rownames(cop_re),decreasing = T),]

pdf('copper.pdf')
drug<-rownames(cop_top)
score<-cop_top[,1]
df<-data.frame(x = drug,y = score)
ggplot(df,aes(x=reorder(drug,score),y=score))+geom_bar(stat = 'identity',width = 0.5)+coord_flip()
dev.off()

#######################################################
flu<-drug.perturbation[,'flucytosine',]
flu<-flu[order(flu[,1],decreasing = T),]
flu_up<-flu[1:10,c('fdr','pvalue')]
flu<-flu[order(flu[,1],decreasing = FALSE),]
flu_down<-flu[1:10,c('fdr','pvalue')]
flucytosine<-c(rep(1,times=10),rep(-1,times=10))
names(flucytosine)<-c(rownames(flu_up),rownames(flu_down))
myfn6<-"CMAP_flu_connectivity.RData"
if(!file.exists(myfn6)){
  message("Be aware that computing sensitivity will take some time...")
  cl<-parallel::makeCluster(nbcore)
  flu_re<-parApply(drug.perturbation[ , ,c("tstat","pvalue")],2,function(x,flucytosine){return(PharmacoGx::connectivityScore(x=x,y=flucytosine,method = "gsea",nperm = 100))},cl=cl,flucytosine=flucytosine)
  stopCluster(cl)
  rownames(flu_re)<-c("Connectivity","P_Value")
  flu_re<-t(flu_re)
  save(flu_re,file=myfn6)
}else{
  load(myfn6)
}
flu_re1<-flu_re[order(flu_re[,1],decreasing = T),]
flu_re1<-flu_re1[1:600,c("Connectivity","P_Value")]
flu_top<-flu_re1[1:10,c("Connectivity","P_Value")]
flu_top<-round(flu_top,digits = 3)

pdf('flucytosine.pdf')
drug<-rownames(flu_top)
score<-flu_top[,1]
df<-data.frame(x = drug,y = score)
ggplot(df,aes(x=reorder(drug,score),y=score))+geom_bar(stat = 'identity',width = 0.5)+coord_flip()
dev.off()

pdf("ketoconazole_table.pdf")

grid.table(ket_top)
dev.off()

pdf("miconazole_table.pdf")
grid.table(mic_top)
dev.off()

pdf("amphotericin B_table.pdf")
grid.table(amB_top)
dev.off()

pdf("nystatin_table.pdf")
grid.table(nys_top)
dev.off()

pdf("copper_table.pdf")
grid.table(cop_top)
dev.off()

pdf("flucytosine_table.pdf")
grid.table(flu_top)
dev.off()


total<-matrix(NA,ncol = dim(ket_re1)[1],nrow = 5,byrow = T)
rownames(total)<-c('ketoconazole','miconazole','amphotericin B','nystatin','flucytosine')
total[1,]<-rownames(ket_re1)
total[2,]<-rownames(mic_re1)
total[3,]<-rownames(amB_re1)
total[4,]<-rownames(nys_re1)
total[5,]<-rownames(flu_re1)
weight<-matrix(NA,ncol = dim(ket_re1)[1],nrow = 5,byrow = T)
rownames(weight)<-c('ketoconazole','miconazole','amphotericin B','nystatin','flucytosine')
weight[1,]<-ket_re1[,1]
weight[2,]<-mic_re1[,1]
weight[3,]<-amB_re1[,1]
weight[4,]<-nys_re1[,1]
weight[5,]<-flu_re1[,1]
r1<-RankAggreg(total,30,weight,method="CE",distance="Kendall",rho=.01, N=10000)
r2<-RankAggreg(total,15,weight,method="CE",distance="Spearman",rho=.01, N=10000)
plot(r1)
plot(r2)
plot(r1["Optimal List"])
ket_re["sertaconazole",]
