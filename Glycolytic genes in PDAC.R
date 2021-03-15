#1
setwd("D:/Yang")
library("org.Hs.eg.db")
rt=read.table("id.txt",sep="\t",check.names=F,header=T)
View(rt)
genes=as.vector(rt[,1])
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs <- as.character(entrezIDs)
out=cbind(rt,entrezID=entrezIDs)
write.table(out,file="id1.txt",sep="\t",quote=F,row.names=F)
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
rt=read.table("id1.txt",sep="\t",header=T,check.names=F)
rt=rt[is.na(rt[,"entrezID"])==F,]
gene=rt$entrezID
kk <- enrichGO(gene = gene,
               OrgDb = org.Hs.eg.db, 
               pvalueCutoff =0.05, 
               qvalueCutoff = 0.05,
               ont="all",
               readable =T)
write.table(kk,file="GO.txt",sep="\t",quote=F,row.names = F)
pdf(file="barplotGO.pdf",width = 10,height = 6)
barplot(kk, drop = TRUE, showCategory =10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()
pdf(file="bubbleGO.pdf",width = 7,height = 6)
dotplot(kk,showCategory = 10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()
rt=read.table("id1.txt",sep="\t",header=T,check.names=F)
rt=rt[is.na(rt[,"entrezID"])==F,]
gene=rt$entrezID
kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =0.05, qvalueCutoff =0.05)
write.table(kk,file="KEGG.txt",sep="\t",quote=F,row.names = F)
pdf(file="barplotKEGG.pdf",width = 10,height = 6)
barplot(kk, drop = TRUE, showCategory = 30)
dev.off()
pdf(file="bubbleKEGG.pdf",width = 7,height = 6)
dotplot(kk, showCategory = 30)
dev.off()
install.packages("GOplot")
install.packages("ggdendro")
library(GOplot)
library(ggplot2)
library(ggdendro)
library(gridExtra)
library(RColorBrewer)
library(stringr)
library(org.Hs.eg.db)
library(clusterProfiler)
gene_symbol<-read.table(file="gene_list.txt",header=T)$SYMBOL
entrez_id <- bitr(gene_symbol, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb=org.Hs.eg.db)
ego_all <- enrichGO(gene =entrez_id$ENTREZID,
                    OrgDb= org.Hs.eg.db,
                    ont = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    readable = TRUE,
                    pool = TRUE)
GO=ego_all[1:6,c(1,2,3,9,7)]
GO$geneID=str_replace_all(GO$geneID,"/",",")
names(GO)=c("Category","ID","term","Genes","adj_pval")
gene=data.frame(ID=gene_symbol,logFC=rnorm(length(gene_symbol),mean=0,sd=2))
gene=data.frame(ID=gene_symbol,logFC=rnorm(length(gene_symbol),mean=0,sd=2))
circ <- circle_dat(GO,gene)
chord <- chord_dat(data = circ, genes = gene, process = GO$term)
GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 5)
out_name="binfohome_goplot_all1"
cairo_pdf(file=paste(out_name,"_GOChord.pdf",seq=""),onefile=TRUE,width=15,height=17)
GOChord(chord,space=0.02,gene.order="logFC",gene.space=0.25,gene.size=5)
dev.off()
#2
library(dplyr)
setwd("D:/Yang")
library(NMF)
setwd("D:/Yang")
exp <- read.csv("D:/Yang/glyexp.csv",header = TRUE,sep = ",",row.names = 1,dec = ".",na.strings = "NA",stringsAsFactors=FALSE,check.names = FALSE)
rownames(exp) <- exp[,1]
exp <- exp[,-1]
exp1<-exp[,!colnames(exp) %in% c("sample20")]
dat=exp1
d=dat
gene_e=319 
mads=apply(d,1,mad)
d=d[rev(order(mads))[1:gene_e],]
d = sweep(d,1, apply(d,1,median,na.rm=T))
library(ConsensusClusterPlus)
#title=tempdir()
capabilities()
options(bitmapType='cairo')
tiff()
results = ConsensusClusterPlus(as.matrix(d),maxK=6,reps=1000,pItem=0.8,pFeature=1,
                               title="D:/Yang",
                               clusterAlg="hc",distance="pearson",seed=12621,plot="tiff",writeTable=F)
#consensusClass - the sample classifications
clu2=as.data.frame(results[[2]][["consensusClass"]])
colnames(clu2)="cluster"
setwd("D:/Yang")
table(clu2)
clu2$ID=row.names(clu2)
clu22 <- clu2[which(clu2$cluster==2),]
clu21 <- clu2[which(clu2$cluster==1),]
write.csv(clu22,file="clu2.csv")
write.csv(clu21,file="clu1.csv")
write.csv(clu2,file="clu.csv")
#3
setwd("D:/Yang")  
time=read.table("time.txt",header=T,sep="\t",check.names=F)
gene <- read.csv("glyexp.csv",header = TRUE,sep = ",",row.names = 1)
row.names(gene) <- gene[,1]
gene <- gene[,-1]
as.data.frame(gene)
colnames(gene)<-substr(colnames(gene),0,12)
gene1 <- t(gene)
gene2=cbind(id=row.names(gene1), gene1)
row.names(gene2)=NULL
colnames(gene2)[1]<-"sample"
sur <- merge(time,gene2,by="sample")
write.table(sur,file="survival.txt",sep="\t",row.names=F,quote=F)
install.packages("survival")
picDir="picture"                                               
dir.create(picDir)
library(survival)
rt=read.table("survival.txt",header=T,sep="\t",check.names=F)
rt$futime=rt$futime/365                                        
outTab=data.frame()
for(gene in colnames(rt[,4:ncol(rt)])){
  a=rt[,gene]<=median(rt[,gene])
  diff=survdiff(Surv(futime, fustat) ~a,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  outTab=rbind(outTab,cbind(gene=gene,pvalue=pValue))
  pValue=round(pValue,3)
  #pValue=format(pValue, scientific = TRUE)
  fit <- survfit(Surv(futime, fustat) ~ a, data = rt)
  summary(fit)
  tiff(file=paste(picDir,"\\",gene,".survival.tiff",sep=""),
       width = 14,            
       height =14,           
       units ="cm",
       compression="lzw",
       bg="white",
       res=600)
  plot(fit, 
       lwd=2,
       col=c("red","blue"),
       xlab="Time (year)",
       mark.time=T,
       ylab="Survival rate",
       main=paste(gene,"(p=", pValue ,")",sep="") )
  legend("topright", 
         c("High","Low"), 
         lwd=2, 
         col=c("red","blue"))
  dev.off()
}
write.table(outTab,file="survival.xls",sep="\t",row.names=F,quote=F)
#4
setwd("D:/Yang")
exp = read.csv(file="tumor.csv",sep = ',',header=F,row.names = 1)
exp1 <- t(exp)
exp1[,1]<-substr(exp1[,1],0,12)
group=read.table("clu.csv",sep = ',',header=T)
colnames(group) <- c("ID","cluster")
group$ID <-substr(group$ID,0,12)
group1 <- paste("cluster", group$cluster, sep ="" ,collapse = NULL)
group$cluster <- group1
time=read.table("time.txt",header=T,sep="\t") 
colnames(time) <- c("ID","fustat","futime")
mer <- merge(group,exp1,by="ID")
mer1 <- merge(time,mer,by="ID")
mer2 <- mer1[order(mer1[,4]),]
mer2$ID=="TCGA.HZ.A9TJ"
mer3 <- mer2[-c(61:63),]
unique(mer3$cluster)
length(which(mer3$cluster=="cluster1"))
length(which(mer3$cluster=="cluster2"))
write.csv(mer3,file="mer3.csv")
exp <- mer3[,-(2:4)]
#colnames(exp)[1] <- NA
rownames(exp) <- exp[,1]
exp <- exp[,-1]
exp_t <- t(exp)
write.csv(exp_t,file="exp_t.csv")
library("dplyr")
library("DESeq2")
library("limma")
library("edgeR")
setwd("D:/Yang")
library("limma")
conNum=103  #cluster1                                                    
treatNum=70 #cluster2
fdrFilter=0.05                                                      
logFCfilter=2                                                     
outTab=data.frame()
grade=c(rep(1,conNum),rep(2,treatNum))
rt=exp_t
exp=rt[,1:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]
data[data<0]=0
for(i in row.names(data)){
  geneName=unlist(strsplit(i,"\\|",))[1]
  geneName=gsub("\\/", "_", geneName)
  rt=rbind(expression=data[i,],grade=grade)
  rt=as.matrix(t(rt))
  wilcoxTest<-wilcox.test(expression ~ grade, data=rt)
  conGeneMeans=mean(data[i,1:conNum])
  treatGeneMeans=mean(data[i,(conNum+1):ncol(data)])
  logFC=log2(treatGeneMeans)-log2(conGeneMeans)
  pvalue=wilcoxTest$p.value
  conMed=median(data[i,1:conNum])
  treatMed=median(data[i,(conNum+1):ncol(data)])
  diffMed=treatMed-conMed
  if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){  
    outTab=rbind(outTab,cbind(gene=i,conMean=conGeneMeans,treatMean=treatGeneMeans,logFC=logFC,pValue=pvalue))
  }
}
pValue=outTab[,"pValue"]
fdr=p.adjust(as.numeric(as.vector(pValue)),method="fdr")
outTab=cbind(outTab,fdr=fdr)
write.csv(outTab,row.names=F,file="all.csv")
outTab <- read.csv("all.csv",header=T,check.names=FALSE)
class(outTab$logFC)
class(outTab$fdr)
outTab$sig[outTab$fdr <= fdrFilter & outTab$logFC >= logFCfilter] <- "up"
outTab$sig[outTab$fdr <= fdrFilter & outTab$logFC <= -logFCfilter] <- "down"
outTab$sig[outTab$fdr > fdrFilter | abs(outTab$logFC) < logFCfilter] <- "no"
write.table(outTab,file="all.xls",sep="\t",row.names=F,quote=F)
outDiff=outTab[( abs(as.numeric(as.vector(outTab$logFC)))>logFCfilter & as.numeric(as.vector(outTab$fdr))<fdrFilter),]
write.table(outDiff,file="diff.xls",sep="\t",row.names=F,quote=F)
heatmap=rbind(ID=colnames(data[as.vector(outDiff[,1]),]),data[as.vector(outDiff[,1]),])
write.table(heatmap,file="tcgaDiffMetabExp.txt",sep="\t",col.names=F,quote=F)
library(ggplot2)
outTab <- outTab[,-8]
p<-ggplot(outTab,aes(x=outTab$logFC,y=-log10(outTab$fdr),colour=sig))+xlab("log2 Fold Change")+ylab("-log10P-Value")+
  geom_point(size=5,alpha=0.6,cex.axis=5, cex.lab=5)+
  theme(axis.title.x = element_text(size=25),
        axis.text.x=element_text(size=25),
        axis.title.y = element_text(size=25),
        axis.text.y = element_text(size=25),
        legend.title=element_text(family="Times",size=25),
        legend.text=element_text(family="Times",size=25))+
  scale_color_manual(values =c("#0072B5","black","#BC3C28"))+
  geom_vline(aes(xintercept=2), colour="gray",size=1.2 ,linetype=2)+
  geom_vline(aes(xintercept=-2), colour="gray",size=1.2 ,linetype=2)+
  geom_hline(aes(yintercept=-log10(0.05)),colour="gray",size=1.2 ,linetype=2) 
guides(fill=guide_legend(title=NULL))
print(p)
ggsave(("vol.tiff"),plot = print(p),width = 10,height = 8,units = "in")
library("pheatmap")
View(data)
outDiff <- outDiff[,-8]
View(outDiff)
outDiff1 <- outDiff[order(outDiff[,7]),]
hmExp=data[as.vector(outDiff1[,1]),]
hmExp=log2(hmExp+0.1)
Type=c(rep("cluster1",conNum),rep("cluster2",treatNum))
names(Type)=colnames(data)
Type=as.data.frame(Type)
pdf(file="heatmap.pdf",height=7,width=11)
pheatmap(hmExp, 
         annotation=Type, 
         color = colorRampPalette(c("green", "black", "red"))(50),
         cluster_cols =F,
         show_colnames = F,
         show_rownames = F,
         fontsize = 12,
         fontsize_row=3,
         fontsize_col=10)
dev.off()
km <- mer3[ , c(1:4)]
write.csv(km,file="km.csv")
rt = read.csv(file="km.csv",sep = ',',header=T,row.names = 1)
install.packages("pheatmap")
library(pheatmap)
library(survival)
library("survminer")
rt$futime=rt$futime/365
diff=survdiff(Surv(futime, fustat) ~cluster,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)
fit <- survfit(Surv(futime, fustat) ~ cluster, data = rt)
pdf(file="survival.pdf",onefile = FALSE,
    width = 5.5,            
    height =5)             
ggsurvplot(fit, 
           data=rt,
           conf.int=TRUE,
           pval=paste0("p=",pValue),
           pval.size=4,
           risk.table=TRUE,
           legend.labs=c("cluster1", "cluster2"),
           legend.title="Cluster",
           xlab="Time(years)",
           break.time.by = 1,
           risk.table.title="",
           palette=c("red", "blue"),
           risk.table.height=.25)
dev.off()
summary(fit)  
#5
library(dplyr)
setwd("D:/Yang") 
gene <- read.table("gene_sig1.txt",header=T,sep="\t",check.names=F) 
exp <- read.table("survival.txt",header=F,sep="\t",check.names=F) 
exp <- t(exp)
colnames(exp) <- exp[1,]
exp <- exp[-1,]
colnames(gene) <- "sample"
as.character(exp[,1])
as.character(gene$sample)
exp <- as.data.frame(exp)
class(exp)
uni <- semi_join(exp,gene,by="sample")
View(uni)
uni <- t(uni)
write.csv(uni,file="uniexp1.csv")
exp <- read.csv("uniexp1.csv")
time=read.table("time.txt",header=T,sep="\t",check.names=F)
rt <- merge(time,exp,by="sample")
write.csv(rt,file="uniexptime1.csv",row.names = F)
library(survival)
rt=read.csv("uniexptime1.csv",header=T,check.names=F,row.names=1)    
pFilter=0.05 
outTab=data.frame()
sigGenes=c("futime","fustat")
for(i in colnames(rt[,3:ncol(rt)])){
  cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  if(coxP<pFilter){
    sigGenes=c(sigGenes,i)
    outTab=rbind(outTab,
                 cbind(id=i,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
    )
  }
}
write.table(outTab,file="UniCox1.txt",sep="\t",row.names=F,quote=F)
uniSigExp=rt[,sigGenes]
uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
write.table(uniSigExp,file="UniSigExp1.txt",sep="\t",row.names=F,quote=F)
rt <- read.table("UniCox1.txt",header=T,sep="\t",row.names=1,check.names=F)
gene <- rownames(rt)
hr <- sprintf("%.3f",rt$"HR")
hrLow  <- sprintf("%.3f",rt$"HR.95L")
hrHigh <- sprintf("%.3f",rt$"HR.95H")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
pdf(file="forest1.pdf", width = 6,height = 10)
n <- nrow(rt)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2))
xlim = c(0,3)
par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)
par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, 'red', 'green')
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
axis(1)
dev.off()
#6
library("glmnet")
library("survival")
setwd("D:/Yang") 
rt=read.table("UniSigExp1.txt",header=T,sep="\t",row.names=1)            
rt$futime=rt$futime/365
x=as.matrix(rt[,c(3:ncol(rt))])
y=data.matrix(Surv(rt$futime,rt$fustat))
fit=glmnet(x, y, family = "cox")
cvfit=cv.glmnet(x, y, family="cox")
pdf(file="lasso1.pdf", width = 8,height = 5)
plot(cvfit) 
dev.off()
coef=coef(fit, s = cvfit$lambda.min)
index=which(coef != 0)
actCoef=coef[index]
lassoGene=row.names(coef)[index]
geneCoef=cbind(Gene=lassoGene,Coef=actCoef)
write.table(geneCoef,file="geneCoef1.txt",sep="\t",quote=F,row.names=F)
trainFinalGeneExp=rt[,lassoGene]
myFun=function(x){crossprod(as.numeric(x),actCoef)}
trainScore=apply(trainFinalGeneExp,1,myFun)
outCol=c("futime","fustat",lassoGene)
risk=as.vector(ifelse(trainScore>median(trainScore),"high","low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(trainScore),risk)
write.table(cbind(id=rownames(outTab),outTab),file="Risk1.txt",sep="\t",quote=F,row.names=F)
#7
library(survival)
library("survminer")
setwd("D:/Yang")              
bioSurvival=function(inputFile=null,outFile=null){
  rt=read.table(inputFile,header=T,sep="\t")                   
  #?Ƚϸߵͷ????????????죬?õ???????pֵ
  diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  pValue=signif(pValue,4)
  pValue=format(pValue, scientific = TRUE)
  fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
  surPlot=ggsurvplot(fit, 
                     data=rt,
                     conf.int=TRUE,
                     pval=paste0("p=",pValue),
                     pval.size=4,
                     risk.table=TRUE,
                     legend.labs=c("High risk", "Low risk"),
                     legend.title="Risk",
                     xlab="Time(years)",
                     break.time.by = 1,
                     risk.table.title="",
                     palette=c("red", "blue"),
                     risk.table.height=.25)
  pdf(file=outFile,onefile = FALSE,width = 6.5,height =5.5)
  print(surPlot)
  dev.off()
}
bioSurvival(inputFile="Risk1.txt",outFile="Risk1.pdf")
library(pheatmap)
bioRiskPlot=function(inputFile=null,riskScoreFile=null,survStatFile=null,heatmapFile=null){
  rt=read.table(inputFile,sep="\t",header=T,row.names=1,check.names=F)  
  rt=rt[order(rt$riskScore),]                                           
  riskClass=rt[,"risk"]
  lowLength=length(riskClass[riskClass=="low"])
  highLength=length(riskClass[riskClass=="high"])
  line=rt[,"riskScore"]
  line[line>10]=10
  pdf(file=riskScoreFile,width = 10,height = 3.5)
  plot(line, type="p", pch=20,
       xlab="Patients (increasing risk socre)", ylab="Risk score",
       col=c(rep("green",lowLength),rep("red",highLength)) )
  abline(h=median(rt$riskScore),v=lowLength,lty=2)
  legend("topleft", c("High risk", "low Risk"),bty="n",pch=19,col=c("red","green"),cex=1.2)
  dev.off()
  color=as.vector(rt$fustat)
  color[color==1]="red"
  color[color==0]="green"
  pdf(file=survStatFile,width = 10,height = 3.5)
  plot(rt$futime, pch=19,
       xlab="Patients (increasing risk socre)", ylab="Survival time (years)",
       col=color)
  legend("topleft", c("Dead", "Alive"),bty="n",pch=19,col=c("red","green"),cex=1.2)
  abline(v=lowLength,lty=2)
  dev.off()
  rt1=rt[c(3:(ncol(rt)-2))]
  rt1=log2(rt1+1)
  rt1=t(rt1)
  annotation=data.frame(type=rt[,ncol(rt)])
  rownames(annotation)=rownames(rt)
  pdf(file=heatmapFile,width = 10,height = 3.5)
  pheatmap(rt1, 
           annotation=annotation, 
           cluster_cols = FALSE,
           fontsize_row=11,
           show_colnames = F,
           fontsize_col=3,
           color = colorRampPalette(c("green", "black", "red"))(50) )
  dev.off()
}
bioRiskPlot(inputFile="Risk1.txt",riskScoreFile="riskScore1.pdf",survStatFile="survStat1.pdf",heatmapFile="heatmap1.pdf")
install.packages("survivalROC")
library(survivalROC)
rt=read.table("Risk1.txt",header=T,sep="\t",check.names=F,row.names=1)    
pdf(file="ROC1.pdf",width=6,height=6)
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, 
                predict.time =1, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
     xlab="False positive rate", ylab="True positive rate",
     main=paste("ROC curve (", "AUC = ",sprintf("%.3f",roc$AUC),")"),
     lwd = 3, cex.main=1.5, cex.lab=1.5, cex.axis=1.5, font=1.5)
abline(0,1)
dev.off()
#8
library(survival)
setwd("D:/Yang")                     
risk=read.table("Risk1.txt",header=T,sep="\t",check.names=F,row.names=1)       
cli=read.table("Clinical2.txt",sep="\t",check.names=F,header=T,row.names=1)     
sameSample=intersect(row.names(cli),row.names(risk))
risk=risk[sameSample,]
cli=cli[sameSample,]
rt=cbind(futime=risk[,1],fustat=risk[,2],cli,risk=risk[,(ncol(risk)-1)])
write.csv(rt,file="mer.csv")

uniTab=data.frame()
for(i in colnames(rt[,3:ncol(rt)])){
  cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  uniTab=rbind(uniTab,
               cbind(id=i,
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
}
write.table(uniTab,file="uniCox1.txt",sep="\t",row.names=F,quote=F)
multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
multiCoxSum=summary(multiCox)
multiTab=data.frame()
multiTab=cbind(
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
multiTab=cbind(id=row.names(multiTab),multiTab)
write.table(multiTab,file="multiCox1.txt",sep="\t",row.names=F,quote=F)
bioForest=function(coxFile=null,forestFile=null,forestCol=null){
  rt <- read.table(coxFile,header=T,sep="\t",row.names=1,check.names=F)
  gene <- rownames(rt)
  hr <- sprintf("%.3f",rt$"HR")
  hrLow  <- sprintf("%.3f",rt$"HR.95L")
  hrHigh <- sprintf("%.3f",rt$"HR.95H")
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
  pdf(file=forestFile, width = 10,height = 4.5)
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2))
  xlim = c(0,3)
  par(mar=c(3,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=1
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
  text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)
  par(mar=c(3,1,2,1),mgp=c(2,0.5,0))
  xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
  abline(v=1,col="black",lty=2,lwd=2)
  boxcolor = ifelse(as.numeric(hr) > 1, forestCol, forestCol)
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
  axis(1)
  dev.off()
}
bioForest(coxFile="uniCox1.txt",forestFile="uniForest1.pdf",forestCol="green")
bioForest(coxFile="multiCox1.txt",forestFile="multiForest1.pdf",forestCol="red")
#9
library(ggplot2)
library(ggpubr)
setwd("D:/Yang")
rt <- read.table(file="rt.txt",header=T,sep="\t",check.names=F)
pdf(file="11gene.pdf",onefile = FALSE,width = 10,height =6) 
ggplot(rt, aes(x = gene, y = exp))+ 
  geom_boxplot(aes(fill = risk),position=position_dodge(0.8),width=0.7)+
  theme_bw() +
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+
  labs(x='Gene', y='Expression', fill='Risk') + 
  theme(axis.text.x = element_text(size=12,color="black",angle = 45,hjust = 0.3, vjust = 0.5), 
        axis.title.x = element_text(size=14,color="black",face = "bold"), 
        axis.title.y = element_text(size=14,color="black",face = "bold"),
        legend.title = element_text(size=14,face = "bold"),
        legend.text = element_text(size=12))+
  stat_compare_means(method = "wilcox.test",aes(group=risk),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
dev.off()
#10
library(ggplot2)
library(plyr)
library(ggthemes)
setwd("D:/Yang")
group=read.table("group.csv",sep = ',',header=T)
cli <- read.table("clinical3.csv",sep = ',',header=T)
mer <- merge(group,cli,by="ID")
write.csv(mer,file="mer.csv")
g <- read.table("grade.csv",sep = ',',header=T,row.names = 1)
levels(g$Grade) <- c("G4","G3","G2","G1")
levels(g$Grade)
pdf(file="grade.pdf",onefile = FALSE,
    width = 3.5,            
    height =6) 
ggplot(g,aes(risk,Percent,fill=Grade))+
  geom_bar(stat="identity",position="stack")+
  scale_fill_manual(values=c("#4169E1","#1E90FF","#00BFFF","#87CEFA"),breaks = c('G4','G3','G2','G1'))+
  ggtitle("")+
  #coord_flip()+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'),
        axis.title.x=element_text(size=14,color="black",hjust=0.5),
        axis.title.y=element_text(size=14,color="black",hjust=0.5),
        axis.text.x=element_text(size=12,color="black"),
        axis.text.y=element_text(size=12,color="black"),
        legend.position="top",
        legend.text = element_text(size = 10, color = "black"),
        legend.title = element_text(size = 11, color = "black"))+
  guides(fill=guide_legend(title="Grade"))
dev.off()
s <- read.table("stage.csv",sep = ',',header=T,row.names = 1)
levels(s$Stage) <- c("IV","III","II","I")
levels(s$Stage)
pdf(file="stage.pdf",onefile = FALSE,
    width = 3.2,            
    height =6) 
ggplot(s,aes(risk,Percent,fill=Stage))+
  geom_bar(stat="identity",position="stack")+
  scale_fill_manual(values=c("#4169E1","#1E90FF","#00BFFF","#87CEFA"),breaks = c('IV','III','II','I'))+
  ggtitle("")+
  #coord_flip()+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'),
        axis.title.x=element_text(size=14,color="black",hjust=0.5),
        axis.title.y=element_text(size=14,color="black",hjust=0.5),
        axis.text.x=element_text(size=12,color="black"),
        axis.text.y=element_text(size=12,color="black"),
        legend.position="top",
        legend.text = element_text(size = 10, color = "black"),
        legend.title = element_text(size = 11, color = "black"))+
  guides(fill=guide_legend(title="Stage"))
dev.off()
t <- read.table("T.csv",sep = ',',header=T,row.names = 1)
levels(t$T.stage) <- c("T4","T3","T2","T1")
levels(t$T.stage)
pdf(file="T.pdf",onefile = FALSE,
    width = 3.5,            
    height =6) 
ggplot(t,aes(risk,Percent,fill=T.stage))+
  geom_bar(stat="identity",position="stack")+
  scale_fill_manual(values=c("#4169E1","#1E90FF","#00BFFF","#87CEFA"),breaks = c('T4','T3','T2','T1'))+
  ggtitle("")+
  #coord_flip()+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'),
        axis.title.x=element_text(size=14,color="black",hjust=0.5),
        axis.title.y=element_text(size=14,color="black",hjust=0.5),
        axis.text.x=element_text(size=12,color="black"),
        axis.text.y=element_text(size=12,color="black"),
        legend.position="top",
        legend.text = element_text(size = 10, color = "black"),
        legend.title = element_text(size = 11, color = "black"))+
  guides(fill=guide_legend(title="T.stage"))
dev.off()
n <- read.table("N.csv",sep = ',',header=T,row.names = 1)
levels(n$N.stage) <- c("N1","N0")
levels(n$N.stage)
pdf(file="N.pdf",onefile = FALSE,
    width = 3,            
    height =6) 
ggplot(n,aes(risk,Percent,fill=N.stage))+
  geom_bar(stat="identity",position="stack")+
  scale_fill_manual(values=c("#4169E1","#00BFFF"),breaks = c('N1','N0'))+
  ggtitle("")+
  #coord_flip()+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'),
        axis.title.x=element_text(size=14,color="black",hjust=0.5),
        axis.title.y=element_text(size=14,color="black",hjust=0.5),
        axis.text.x=element_text(size=12,color="black"),
        axis.text.y=element_text(size=12,color="black"),
        legend.position="top",
        legend.text = element_text(size = 10, color = "black"),
        legend.title = element_text(size = 11, color = "black"))+
  guides(fill=guide_legend(title="N.stage"))
dev.off()
m <- read.table("M.csv",sep = ',',header=T,row.names = 1)
levels(m$M.stage) <- c("M1","M0")
levels(m$M.stage)
pdf(file="M.pdf",onefile = FALSE,
    width = 3,            
    height =6) 
ggplot(m,aes(risk,Percent,fill=M.stage))+
  geom_bar(stat="identity",position="stack")+
  scale_fill_manual(values=c("#4169E1","#00BFFF"),breaks = c('M1','M0'))+
  ggtitle("")+
  #coord_flip()+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'),
        axis.title.x=element_text(size=14,color="black",hjust=0.5),
        axis.title.y=element_text(size=14,color="black",hjust=0.5),
        axis.text.x=element_text(size=12,color="black"),
        axis.text.y=element_text(size=12,color="black"),
        legend.position="top",
        legend.text = element_text(size = 10, color = "black"),
        legend.title = element_text(size = 11, color = "black"))+
  guides(fill=guide_legend(title="M.stage"))
dev.off()
setwd("D:/Yang")
group=read.table("group.csv",sep = ',',header=T)
clu=read.table("clu.csv",sep = ',',header=T)
clu[,1]<-substr(clu[,1],0,12)
mer <- merge(clu,group,by="ID")
write.csv(mer,file="mer.csv")
#11
setwd("D:/Yang")
exp = read.csv(file="tumor.csv",row.names=1,sep = ',',header=T) 
exp1 <- t(exp)
rownames(exp1)<-substr(rownames(exp1),0,12)
write.csv(exp1,file="tumor1.csv")
exp = read.csv(file="tumor1.csv",sep = ',',header=T) 
colnames(exp)[1]<-"id"
group=read.table("group.csv",sep = ',',header=T)
mer <- merge(group,exp,by="id")
mer1 <- mer[order(mer[,2]),]
unique(mer$risk)
length(which(mer$risk=="high"))
length(which(mer$risk=="low"))
write.csv(mer1,file="mer1.csv")
exp <- mer1[,-2]
colnames(exp)[1] <- NA
rownames(exp) <- exp[,1]
exp <- exp[,-1]
exp_t <- t(exp)
write.csv(exp_t,file="exp_t.csv")
library("dplyr")
library("DESeq2")
library("limma")
library("edgeR")
setwd("D:/Yang")
library("limma")
conNum=86   #high                                                    
treatNum=87 #low
fdrFilter=0.05                                                      
logFCfilter=2                                                     
outTab=data.frame()
grade=c(rep(1,conNum),rep(2,treatNum))
rt=exp_t
exp=rt[,1:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]
data[data<0]=0
for(i in row.names(data)){
  geneName=unlist(strsplit(i,"\\|",))[1]
  geneName=gsub("\\/", "_", geneName)
  rt=rbind(expression=data[i,],grade=grade)
  rt=as.matrix(t(rt))
  wilcoxTest<-wilcox.test(expression ~ grade, data=rt)
  conGeneMeans=mean(data[i,1:conNum])
  treatGeneMeans=mean(data[i,(conNum+1):ncol(data)])
  logFC=log2(treatGeneMeans)-log2(conGeneMeans)
  pvalue=wilcoxTest$p.value
  conMed=median(data[i,1:conNum])
  treatMed=median(data[i,(conNum+1):ncol(data)])
  diffMed=treatMed-conMed
  if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){  
    outTab=rbind(outTab,cbind(gene=i,conMean=conGeneMeans,treatMean=treatGeneMeans,logFC=logFC,pValue=pvalue))
  }
}
pValue=outTab[,"pValue"]
fdr=p.adjust(as.numeric(as.vector(pValue)),method="fdr")
outTab=cbind(outTab,fdr=fdr)
write.csv(outTab,row.names=F,file="all1.csv")
outTab <- read.csv("all.csv",header=T,check.names=FALSE)
class(outTab$logFC)
class(outTab$fdr)
outTab$sig[outTab$fdr <= fdrFilter & outTab$logFC >= logFCfilter] <- "up"
outTab$sig[outTab$fdr <= fdrFilter & outTab$logFC <= -logFCfilter] <- "down"
outTab$sig[outTab$fdr > fdrFilter | abs(outTab$logFC) < logFCfilter] <- "no"
write.table(outTab,file="all.xls",sep="\t",row.names=F,quote=F)
outDiff=outTab[( abs(as.numeric(as.vector(outTab$logFC)))>logFCfilter & as.numeric(as.vector(outTab$fdr))<fdrFilter),]
write.table(outDiff,file="diff.xls",sep="\t",row.names=F,quote=F)
heatmap=rbind(ID=colnames(data[as.vector(outDiff[,1]),]),data[as.vector(outDiff[,1]),])
write.table(heatmap,file="tcgaDiffMetabExp.txt",sep="\t",col.names=F,quote=F)
library(ggplot2)
p<-ggplot(outTab,aes(x=outTab$logFC,y=-log10(outTab$fdr),colour=sig))+xlab("log2 Fold Change")+ylab("-log10P-Value")+
  geom_point(size=5,alpha=0.6,cex.axis=5, cex.lab=5)+
  theme(axis.title.x = element_text(size=25),
        axis.text.x=element_text(size=25),
        axis.title.y = element_text(size=25),
        axis.text.y = element_text(size=25),
        legend.title=element_text(family="Times",size=25),
        legend.text=element_text(family="Times",size=25))+
  scale_color_manual(values =c("#0072B5","black","#BC3C28"))+
  geom_vline(aes(xintercept=2), colour="gray",size=1.2 ,linetype=2)+
  geom_vline(aes(xintercept=-2), colour="gray",size=1.2 ,linetype=2)+
  geom_hline(aes(yintercept=-log10(0.05)),colour="gray",size=1.2 ,linetype=2) 
guides(fill=guide_legend(title=NULL))
print(p)
ggsave(("vol1.tiff"),plot = print(p),width = 10,height = 8,units = "in")
library("pheatmap")
View(data)
View(outDiff)
outDiff1 <- outDiff[order(outDiff[,7]),]
hmExp=data[as.vector(outDiff1[,1]),]
hmExp=log2(hmExp+0.1)
Type=c(rep("high",conNum),rep("low",treatNum))
names(Type)=colnames(data)
Type=as.data.frame(Type)
pdf(file="heatmap.pdf",height=7,width=12)
pheatmap(hmExp, 
         annotation=Type, 
         color = colorRampPalette(c("green", "black", "red"))(50),
         cluster_cols =F,
         show_colnames = F,
         show_rownames = F,
         fontsize = 12,
         fontsize_row=3,
         fontsize_col=10)
dev.off()
#12
setwd("D:/Yang")
library("org.Hs.eg.db")
rt=read.table("id.txt",sep="\t",check.names=F,header=T)
View(rt)
genes=as.vector(rt[,1])
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs <- as.character(entrezIDs)
out=cbind(rt,entrezID=entrezIDs)
write.table(out,file="id1.txt",sep="\t",quote=F,row.names=F)
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
rt=read.table("id1.txt",sep="\t",header=T,check.names=F)
rt=rt[is.na(rt[,"entrezID"])==F,]
gene=rt$entrezID
kk <- enrichGO(gene = gene,
               OrgDb = org.Hs.eg.db, 
               pvalueCutoff =0.05, 
               qvalueCutoff = 0.05,
               ont="all",
               readable =T)
write.table(kk,file="GO.txt",sep="\t",quote=F,row.names = F)
pdf(file="barplotGO.pdf",width = 12,height = 7)
barplot(kk, drop = TRUE, showCategory =10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()
pdf(file="bubbleGO.pdf",width = 12,height = 7)
dotplot(kk,showCategory = 10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()
rt=read.table("id1.txt",sep="\t",header=T,check.names=F)
rt=rt[is.na(rt[,"entrezID"])==F,]
gene=rt$entrezID
kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =0.05, qvalueCutoff =0.05)
write.table(kk,file="KEGG.txt",sep="\t",quote=F,row.names = F)
pdf(file="barplotKEGG.pdf",width = 7,height = 6)
barplot(kk, drop = TRUE, showCategory = 30)
dev.off()
pdf(file="bubbleKEGG.pdf",width = 7,height = 6)
dotplot(kk, showCategory = 30)
dev.off()
install.packages("GOplot")
install.packages("ggdendro")
library(GOplot)
library(ggplot2)
library(ggdendro)
library(gridExtra)
library(RColorBrewer)
library(stringr)
library(org.Hs.eg.db)
library(clusterProfiler)
gene_symbol<-read.table(file="gene_list.txt",header=T)$SYMBOL
entrez_id <- bitr(gene_symbol, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb=org.Hs.eg.db)
ego_all <- enrichGO(gene =entrez_id$ENTREZID,
                    OrgDb= org.Hs.eg.db,
                    ont = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    readable = TRUE,
                    pool = TRUE)
GO=ego_all[1:6,c(1,2,3,9,7)]
GO$geneID=str_replace_all(GO$geneID,"/",",")
names(GO)=c("Category","ID","term","Genes","adj_pval")
gene=data.frame(ID=gene_symbol,logFC=rnorm(length(gene_symbol),mean=0,sd=2))
gene=data.frame(ID=gene_symbol,logFC=rnorm(length(gene_symbol),mean=0,sd=2))
circ <- circle_dat(GO,gene)
chord <- chord_dat(data = circ, genes = gene, process = GO$term)
GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 5)
out_name="binfohome_goplot_all1"
cairo_pdf(file=paste(out_name,"_GOChord.pdf",seq=""),onefile=TRUE,width=15,height=17)
GOChord(chord,space=0.02,gene.order="logFC",gene.space=0.25,gene.size=5)
dev.off()
#13
library(dplyr)
setwd("D:/Yang")  
gene=read.table("gene1.txt",header=T,sep="\t",check.names=F)
exp <- read.csv("glyexp.csv",header = TRUE,sep = ",",row.names = 1)
a <- semi_join(exp,gene,by="ID")
rownames(a) <- a[,1]
a <- a[,-1]
b <- t(a)
write.csv(b,file="cor1.csv")
M<-cor(b)
head(M)
write.csv(M,file="cor2.csv")
library(corrplot)
library(RColorBrewer)
pdf(file="corplot1.pdf", width = 6,height = 6)
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
# matrix of the p-value of the correlation
p.mat <- cor.mtest(b)
head(p.mat[, 1:5])
write.csv(p.mat,file="cor3.csv")
corrplot(M, type="upper", order="hclust", 
         p.mat = p.mat, sig.level = 0.05,
         col=brewer.pal(n=8, name="RdBu"),
         tl.col="black", tl.srt=90)
dev.off()
#14
install.packages("estimate")
library(limma)
library(estimate)
setwd("D:/Yang") 
rt=read.csv("tumor.csv",header=T,check.names=F)
colnames(rt)<-substr(colnames(rt),0,12)
write.table(rt,file="uniq.symbol.txt",sep="\t",quote=F,col.names=T,row.names = F)
filterCommonGenes(input.f="uniq.symbol.txt", 
                  output.f="commonGenes.gct", 
                  id="GeneSymbol")
estimateScore(input.ds = "commonGenes.gct",
              output.ds="estimateScore.gct")
scores=read.table("estimateScore.gct",skip = 2,header = T)
rownames(scores)=scores[,1]
scores=t(scores[,3:ncol(scores)])
rownames(scores)=gsub("\\.","\\-",rownames(scores))
out=rbind(ID=colnames(scores),scores)
write.table(out,file="scores.txt",sep="\t",quote=F,col.names=F)
risk <- read.table(file="risk1.txt",header=T,sep="\t",check.names=F)
risk1 <- risk[,-c(2:15)]
score <- read.table(file="scores.txt",header=T,sep="\t",check.names=F)
colnames(risk1) <- c("ID","risk")
mer <- merge(risk1,score,by="ID")
write.csv(mer,file="mer.csv")
pdf(file="StromalScore.pdf",onefile = FALSE,width = 5,height =5) 
boxplot(StromalScore~risk,mer,width=c(1,1), 
        outline=F,
        lwd=2,
        col = c("#E69F00","#56B4E9"),
        border=c("black","black"),
        cex = 1.5,
        cex.axis=1.3,
        cex.lab=1.3)
dev.off()
pdf(file="ImmuneScore.pdf",onefile = FALSE,width = 5,height =5) 
boxplot(ImmuneScore~risk,mer,width=c(1,1), 
        outline=F,
        lwd=2,
        col = c("#E69F00", "#56B4E9"),
        border=c("black","black"),
        cex = 1.5,
        cex.axis=1.3,
        cex.lab=1.3)
dev.off()
pdf(file="ESTIMATEScore.pdf",onefile = FALSE,width = 5,height =5) 
boxplot(ESTIMATEScore~risk,mer,width=c(1,1), 
        outline=F,
        lwd=2,
        col = c("#E69F00", "#56B4E9"),
        border=c("black","black"),
        cex = 1.5,
        cex.axis=1.3,
        cex.lab=1.3)
dev.off()
pdf(file="TumorPurity.pdf",onefile = FALSE,width = 5,height =5) 
boxplot(TumorPurity~risk,mer,width=c(1,1), 
        outline=F,
        lwd=2,
        col = c("#E69F00", "#56B4E9"),
        border=c("black","black"),
        cex = 1.5,
        cex.axis=1.3,
        cex.lab=1.3)
dev.off()
#15
install.packages('e1071')
library("limma")                                                    
setwd("D:/Yang") 
rt=read.csv("tumor.csv",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]
v <-voom(data, plot = F, save.plot = F)
out=v$E
out=rbind(ID=colnames(out),out)
write.table(out,file="uniq.symbol.txt",sep="\t",quote=F,col.names=F)       
source("ssGSEA18.CIBERSORT.R")
results=CIBERSORT("ref.txt", "uniq.symbol.txt", perm=100, QN=TRUE)
a <- read.table(file="risk1.txt",header=T,sep="\t",check.names=F)
b <- a[,-(2:15)]
colnames(b) <- c("id","risk")
x <- read.table(file="CIBERSORT-Results.txt",header=T,sep="\t",check.names=F,row.names=1)
rownames(x)<-substr(rownames(x),0,12)
id <- rownames(x)
z <- cbind(id,x)
t <- merge(z,b,by="id")
rownames(t) <- t[,1]
t <- t[,-1]
t=t[order(t[,26]),]
p=t[t[,"P-value"]<0.05,]
p=t[,1:(ncol(t)-4)]
m <- t(p)
rt <- m
group <- t[,-c(1:25)]
class(group)
group <- as.data.frame(group)
rownames(group) <- rownames(t)
library(pheatmap)
pdf("heatmap.pdf",height=5,width=9)
pheatmap(rt, annotation=group, 
         color = colorRampPalette(c("blue",  "red"))(50),
         cluster_cols =F,
         fontsize=8,
         fontsize_row=8,
         scale="row",
         show_colnames=F,
         fontsize_col=3)
dev.off()
install.packages("ggpubr")
library(ggpubr)
pFilter=0.05
f <- group[,1]
group1 <- cbind(f,group)
colnames(group1)=c("cluster","Group")
outTab=data.frame()
data=cbind(p,group1)
for(i in colnames(data[,1:(ncol(data)-2)])){
  rt1=data[,c(i,"Group")]
  colnames(rt1)=c("expression","Group")
  ksTest<-wilcox.test(expression ~ Group, data = rt1)
  pValue=ksTest$p.value
  if(pValue<pFilter){
    outTab=rbind(outTab,cbind(rt1,gene=i))
    print(pValue)
  }
}
write.table(outTab,file="data.txt",sep="\t",row.names=F,quote=F)
data=read.table("data.txt",sep="\t",header=T,check.names=F)       
p=ggboxplot(data, x="gene", y="expression", color = "Group",
            ylab="Fraction",
            xlab="",
            palette = c("blue","red") )
p=p+rotate_x_text(45)
pdf(file="boxplot.pdf",width=7.5,height=5.5)                         
p+stat_compare_means(aes(group=Group),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
dev.off()
#16
library("limma")  
library(dplyr)
setwd("D:/Yang")
gpl <- read.table(file="GPL.txt",header=T,sep="\t",check.names=F)
class(gpl)
x <- as.character(gpl$gene_assignment)
ID <- sapply(strsplit(x, split = "//"),"[",2)
library(stringr)
ID <- str_trim(ID,'right')
data<-cbind(ID,gpl)
data1 <- data[,-3]
colnames(data1) <- c("ID","id")
write.table(data1,file="GPL1.txt",sep="\t",row.names = F,na = "")
exp <- read.table(file="exp.txt",header=T,sep="\t",check.names=F)
mer <- merge(data1,exp,by="id")
mer1 <- mer[,-1]
mer2 <- mer1[complete.cases(mer1[,1]),]
write.table(mer2,file="exp_id.txt",sep="\t",row.names = F,na = "")
rt <- read.table(file="exp_id.txt",header=T,sep="\t",check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
write.table(data,file="exp_id.txt",sep="\t",row.names = T,na = "")
exp <- read.table(file="exp_id.txt",header=T,sep="\t",check.names=F)
gene <- read.table(file="gene1.txt",header=T,sep="\t",check.names=F)
exp11 <- semi_join(exp,gene,by="ID")
write.table(exp11,file="exp11.txt",sep="\t",row.names = F)
exp <- read.table(file="exp11_t.txt",header=T,sep="\t",check.names=F)
cli <- read.table(file="clinical_t.txt",header=T,sep="\t",check.names=F)
mer <- merge(cli,exp,by="ID")
mer$futime <- mer$futime/12
write.table(mer,file="sur.txt",sep="\t",row.names = F)
library(survival)
library("survminer")
bioSurvival=function(inputFile=null,outFile=null){
  rt=read.table(inputFile,header=T,sep="\t")                   
  diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  pValue=signif(pValue,4)
  pValue=format(pValue, scientific = TRUE)
  fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
  surPlot=ggsurvplot(fit, 
                     data=rt,
                     conf.int=TRUE,
                     pval=paste0("p=",pValue),
                     pval.size=4,
                     risk.table=TRUE,
                     legend.labs=c("High risk", "Low risk"),
                     legend.title="Risk",
                     xlab="Time(years)",
                     break.time.by = 1,
                     risk.table.title="",
                     palette=c("red", "blue"),
                     risk.table.height=.25)
  pdf(file=outFile,onefile = FALSE,width = 6.5,height =5.5)
  print(surPlot)
  dev.off()
}
bioSurvival(inputFile="sur1.txt",outFile="Risk.pdf")
a <- read.table(file="sur1.txt",header=T,sep="\t",check.names=F)
class(a)
library(pheatmap)
bioRiskPlot=function(inputFile=null,riskScoreFile=null,survStatFile=null,heatmapFile=null){
  rt=read.table(inputFile,sep="\t",header=T,row.names=1,check.names=F)  
  rt=rt[order(rt$riskScore),]                                           
  riskClass=rt[,"risk"]
  lowLength=length(riskClass[riskClass=="low"])
  highLength=length(riskClass[riskClass=="high"])
  line=rt[,"riskScore"]
  line[line>10]=10
  pdf(file=riskScoreFile,width = 10,height = 3.5)
  plot(line, type="p", pch=20,
       xlab="Patients (increasing risk socre)", ylab="Risk score",
       col=c(rep("green",lowLength),rep("red",highLength)) )
  abline(h=median(rt$riskScore),v=lowLength,lty=2)
  legend("topleft", c("High risk", "low Risk"),bty="n",pch=19,col=c("red","green"),cex=1.2)
  dev.off()
  color=as.vector(rt$fustat)
  color[color==1]="red"
  color[color==0]="green"
  pdf(file=survStatFile,width = 10,height = 3.5)
  plot(rt$futime, pch=19,
       xlab="Patients (increasing risk socre)", ylab="Survival time (years)",
       col=color)
  legend("topleft", c("Dead", "Alive"),bty="n",pch=19,col=c("red","green"),cex=1.2)
  abline(v=lowLength,lty=2)
  dev.off()
  rt1=rt[c(3:(ncol(rt)-2))]
  rt1=log2(rt1+1)
  rt1=t(rt1)
  annotation=data.frame(type=rt[,ncol(rt)])
  rownames(annotation)=rownames(rt)
  pdf(file=heatmapFile,width = 10,height = 3.5)
  pheatmap(rt1, 
           annotation=annotation, 
           cluster_cols = FALSE,
           fontsize_row=11,
           show_colnames = F,
           fontsize_col=3,
           color = colorRampPalette(c("green", "black", "red"))(50) )
  dev.off()
}
bioRiskPlot(inputFile="sur1.txt",riskScoreFile="riskScore.pdf",survStatFile="survStat.pdf",heatmapFile="heatmap.pdf")
install.packages("survivalROC")
library(survivalROC)
rt=read.table("sur1.txt",header=T,sep="\t",check.names=F,row.names=1)    
pdf(file="ROC.pdf",width=6,height=6)
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, 
                predict.time =1, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
     xlab="False positive rate", ylab="True positive rate",
     main=paste("ROC curve (", "AUC = ",sprintf("%.3f",roc$AUC),")"),
     lwd = 3, cex.main=1.5, cex.lab=1.5, cex.axis=1.5, font=1.5)
abline(0,1)
dev.off()