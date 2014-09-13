#DSSIRDV
#=======

#Different Successional Stages' Interspecific Relationships of Degraded Vegetation

#将数据放在目录E:\RWD\4-maotai\1-article下

#设定工作空间
setwd("E:\\RWD\\4-maotai\\1-article")

#一、划分演替阶段
#载入vegan程序包
#若无则install.packages("vegan")
library(vegan)

##读取数据保存到varespec、varechem
spe <- read.csv("species.csv",header=T,row.names=1 )
env <- read.csv("site.csv",header=T ,row.names=1)
N1 <- dim(spe)

#筛选优势种
selectspe <- as.array(unlist(spe))
#排除0对分位数计算的影响
selectspe[selectspe==0]=NA
Delset <- as.array(unlist(-1*which(is.na(selectspe))))
selectspe <- selectspe[Delset]

#使用重要值的期望（均值)作为筛选优势物种的阈值
IVmean <- summary(selectspe)
#按（列)物种，将所有样地中重要值
#最大值低于重要值期望（均值)的物种数据删掉
Delset <- as.array(unlist(-1*which(apply(spe,2,max)<IVmean[4])))
spe <- spe[,Delset]
#查看筛选后的数据及行列
spe
N2 <- dim(spe)

#绘制重要值数据直方图图
dev.new(width=3,height=3)
hist(selectspe,breaks=25, freq=FALSE,col="gray",cex=3,
xlab="Important Value",main="Distribution of Important Value")
lines(density(selectspe,from=0,to=1),col="red",lwd=2)
abline(v=IVmean[4],col="blue",lty=2,lwd=2)
text(IVmean[4],7,paste("Mean=",IVmean[4]),pos=4)
savePlot(filename="Distribution of Important Value.pdf",type=c("pdf"))


#基础数据转换
#样地距离分析
spe.norm <- decostand(spe,"normalize")
spe.bray <- vegdist(spe.norm,"bray")
spe.ch <- vegdist(spe.norm, "euc")
#物种距离分析
site <- t(spe)
site.norm <- decostand(site,"normalize")
site.bray <- vegdist(site.norm,"bray")
site.ch <- vegdist(site.norm, "euc")

# 依据轮廓宽度图选择最优化的聚类簇数量（Rousseeuw质量指数)
# 绘制所有分类水平（除了k=1组的情况)轮廓宽度值（Ward 聚类)
# 首先产生一个长度等于样方数量的空向量asw
asw <- numeric(nrow(spe))
# 其次循环获得轮廓宽度值并依次填入asw向量
# 利用cluster程序包内函数plot.silhouette（)绘制轮廓宽度值k
#载入cluster程序包
#若无则install.packages("cluster")后再执行载入命令
library(cluster)

spe.bray.ward <- hclust(spe.bray, "ward.D")
site.bray.ward <- hclust(site.bray, "ward.D")

#将聚类结果重排，需要gclus程序包，若无则install.packages("gclus")
library(gclus)
spe.chwo<- reorder.hclust(spe.bray.ward,spe.ch)
site.chwo<- reorder.hclust(site.bray.ward,site.ch)


#寻找最适聚类族
for (k in 2:(nrow(spe)-1)) {
	sil <- silhouette(cutree(spe.chwo , k=k), spe.ch)
	asw[k] <- summary(sil)$avg.width
	}
# 选择最佳（最大)轮廓宽度值
k.best <- which.max(asw)
dev.new(width=3,height=3)
plot(1:nrow(spe), asw, type="h",
	xlab="k (groups)", ylab="Average Silhouette Width",
	xaxs="i",yaxs="i",cex=3,ylim=c(0,0.3),bty="L")
axis(1, k.best, paste("Optimum",k.best,sep="\n"), col="red", font=1,
  col.axis="red")
points(k.best, max(asw), pch=16, col="red", cex=0.5)
savePlot(filename="Silhouette Width.pdf",type=c("pdf"))

#分析得出最有分组数及其轮廓宽度值
cat("", "Silhouette-optimal number of clusters k =", k.best, "\n", 
	"with an average silhouette width of", max(asw), "\n")
###*******结合实际得出******适合分组数
#及其轮廓宽度值（******此处“k”自定义)，#实际分为群丛等级
k=16
asw[k]

#划分类型
spe.bw.groups <- cutree(spe.chwo, k=k)
write.csv(spe.bw.groups,file="planttype.csv")

dev.new(width=9.7,height=6)
plot(spe.chwo,hang=-1)
rect.hclust(spe.chwo, k=k, border = "blue")

#将树状图和热力图结合
dend<- as.dendrogram(spe.chwo)

#Hill双向指示种分析法
or<- vegemite(spe,spe.chwo,"Hill",zero="-")
#需要程序包RColorBrewer,若无则install.packages("RColorBrewer")
library(RColorBrewer)
dev.new(width=6,height=9.7)
heatmap(t(spe[rev(or$species)]),Rowv=NA,Colv=dend,
col=brewer.pal(9,"YlOrRd"),scale="none",margin=c(4,4),
ylab="Species name abbreviations",xlab="sites",cex=0.5,lwd=0.5,
cexRow=0.3,cexCol=0.4)
savePlot(filename="clusters.pdf",type=c("pdf"))

#计算样方聚类后物种重要值均值以命名群落
a0 <- merge(spe,spe.bw.groups,by.x="row.names",by.y="row.names",sort=FALSE)
a <- a0[,-1]
c0=aggregate(a,by=list(a$y),mean)
c=(data.frame(t(c0[,c(-1,-length(c0))])))
names(c)=unlist(c0$y)
c
#存储群丛重要值表
write.csv(c,file="planttype_importance.csv")

#二、计算种间关系

#读取种间关系分析数据（要求species作为列名排列部分与spe列名同)
#含有已划分类型的type字段
spens <- read.csv("splots-species-numberdata.csv",header=T)
(N3 <- dim(spens))
spens <- spens[,Delset]
(N4 <- dim(spens))

#按类群整理数据
Type1=subset(spens,type==1,select=-c(type,plots,splots))
Dataset1 <- as.array(unlist(-1*which(apply(Type1,2,max)==0)))
Type1 <- as.matrix(Type1[,Dataset1])

Type2=subset(spens,type==2,select=-c(type,plots,splots))
Dataset2 <- as.array(unlist(-1*which(apply(Type2,2,max)==0)))
Type2 <- as.matrix(Type2[,Dataset2])

Type3=subset(spens,type==3,select=-c(type,plots,splots))
Dataset3 <- as.array(unlist(-1*which(apply(Type3,2,max)==0)))
Type3 <- as.matrix(Type3[,Dataset3])

Type4=subset(spens,type==4,select=-c(type,plots,splots))
Dataset4 <- as.array(unlist(-1*which(apply(Type4,2,max)==0)))
Type4 <- as.matrix(Type4[,Dataset4])

Type5=subset(spens,type==5,select=-c(type,plots,splots))
Dataset5 <- as.array(unlist(-1*which(apply(Type5,2,max)==0)))
Type5 <- as.matrix(Type5[,Dataset5])

Type6=subset(spens,type==6,select=-c(type,plots,splots))
Dataset6 <- as.array(unlist(-1*which(apply(Type6,2,max)==0)))
Type6 <- as.matrix(Type6[,Dataset6])

Type7=subset(spens,type==7,select=-c(type,plots,splots))
Dataset7 <- as.array(unlist(-1*which(apply(Type7,2,max)==0)))
Type7 <- as.matrix(Type7[,Dataset7])

Type8=subset(spens,type==8,select=-c(type,plots,splots))
Dataset8 <- as.array(unlist(-1*which(apply(Type8,2,max)==0)))
Type8 <- as.matrix(Type8[,Dataset8])

Type9=subset(spens,type==9,select=-c(type,plots,splots))
Dataset9 <- as.array(unlist(-1*which(apply(Type9,2,max)==0)))
Type9 <- as.matrix(Type9[,Dataset9])

Type10=subset(spens,type==10,select=-c(type,plots,splots))
Dataset10 <- as.array(unlist(-1*which(apply(Type10,2,max)==0)))
Type10 <- as.matrix(Type10[,Dataset10])

Type11=subset(spens,type==11,select=-c(type,plots,splots))
Dataset11 <- as.array(unlist(-1*which(apply(Type11,2,max)==0)))
Type11 <- as.matrix(Type11[,Dataset11])

Type12=subset(spens,type==12,select=-c(type,plots,splots))
Dataset12 <- as.array(unlist(-1*which(apply(Type12,2,max)==0)))
Type12 <- as.matrix(Type12[,Dataset12])

Type13=subset(spens,type==13,select=-c(type,plots,splots))
Dataset13 <- as.array(unlist(-1*which(apply(Type13,2,max)==0)))
Type13 <- as.matrix(Type13[,Dataset13])

Type14=subset(spens,type==14,select=-c(type,plots,splots))
Dataset14 <- as.array(unlist(-1*which(apply(Type14,2,max)==0)))
Type14 <- as.matrix(Type14[,Dataset14])

Type15=subset(spens,type==15,select=-c(type,plots,splots))
Dataset15 <- as.array(unlist(-1*which(apply(Type15,2,max)==0)))
Type15 <- as.matrix(Type15[,Dataset15])

Type16=subset(spens,type==16,select=-c(type,plots,splots))
Dataset16 <- as.array(unlist(-1*which(apply(Type16,2,max)==0)))
Type16 <- as.matrix(Type16[,Dataset16])

#加载种间关系计算程序包spaa
#若无则加载命令install.packages("spaa")
library(spaa)

A1 <- sp.assoc(Type1)
A2 <- sp.assoc(Type2)
A3 <- sp.assoc(Type3)
A4 <- sp.assoc(Type4)
A5 <- sp.assoc(Type5)
A6 <- sp.assoc(Type6)
A7 <- sp.assoc(Type7)
A8 <- sp.assoc(Type8)
A9 <- sp.assoc(Type9)
A10 <- sp.assoc(Type10)
A11 <- sp.assoc(Type11)
A12 <- sp.assoc(Type12)
A13 <- sp.assoc(Type13)
A14 <- sp.assoc(Type14)
A15 <- sp.assoc(Type15)
A16 <- sp.assoc(Type16)


#计算种间协变
sp.pair1 <- cor(Type1, method = "spearman")

posass1 <- (sum(sp.pair1>0)-ncol(sp.pair1))/2
negass1 <- (sum(sp.pair1<0))/2
PNratio1 <- posass1/negass1

sp.pair2 <- cor(Type2, method = "spearman")

posass2 <- (sum(sp.pair2>0)-ncol(sp.pair2))/2
negass2 <- (sum(sp.pair2<0))/2
PNratio2 <- posass2/negass2

sp.pair3 <- cor(Type3, method = "spearman")

posass3 <- (sum(sp.pair3>0)-ncol(sp.pair3))/2
negass3 <- (sum(sp.pair3<0))/2
PNratio3 <- posass3/negass3

sp.pair4 <- cor(Type4, method = "spearman")

posass4 <- (sum(sp.pair4>0)-ncol(sp.pair4))/2
negass4 <- (sum(sp.pair4<0))/2
PNratio4 <- posass1/negass4

sp.pair5 <- cor(Type5, method = "spearman")

posass5 <- (sum(sp.pair5>0)-ncol(sp.pair5))/2
negass5 <- (sum(sp.pair5<0))/2
PNratio5 <- posass5/negass5

sp.pair6 <- cor(Type6, method = "spearman")

posass6 <- (sum(sp.pair6>0)-ncol(sp.pair6))/2
negass6 <- (sum(sp.pair6<0))/2
PNratio6 <- posass6/negass6

sp.pair7 <- cor(Type7, method = "spearman")

posass7 <- (sum(sp.pair7>0)-ncol(sp.pair7))/2
negass7 <- (sum(sp.pair7<0))/2
PNratio7 <- posass7/negass7

sp.pair8 <- cor(Type8, method = "spearman")

posass8 <- (sum(sp.pair8>0)-ncol(sp.pair8))/2
negass8 <- (sum(sp.pair8<0))/2
PNratio8 <- posass8/negass8

sp.pair9 <- cor(Type9, method = "spearman")

posass9 <- (sum(sp.pair9>0)-ncol(sp.pair9))/2
negass9 <- (sum(sp.pair9<0))/2
PNratio9 <- posass9/negass9

sp.pair10 <- cor(Type10, method = "spearman")

posass10 <- (sum(sp.pair10>0)-ncol(sp.pair10))/2
negass10 <- (sum(sp.pair10<0))/2
PNratio10 <- posass10/negass10

sp.pair11 <- cor(Type11, method = "spearman")

posass11 <- (sum(sp.pair11>0)-ncol(sp.pair11))/2
negass11 <- (sum(sp.pair11<0))/2
PNratio11 <- posass11/negass11

sp.pair12 <- cor(Type12, method = "spearman")

posass12 <- (sum(sp.pair12>0)-ncol(sp.pair12))/2
negass12 <- (sum(sp.pair12<0))/2
PNratio12 <- posass12/negass12

sp.pair13 <- cor(Type13, method = "spearman")

posass13 <- (sum(sp.pair13>0)-ncol(sp.pair13))/2
negass13 <- (sum(sp.pair13<0))/2
PNratio13 <- posass13/negass13

sp.pair14 <- cor(Type14, method = "spearman")

posass14 <- (sum(sp.pair14>0)-ncol(sp.pair14))/2
negass14 <- (sum(sp.pair14<0))/2
PNratio14 <- posass14/negass14

sp.pair15 <- cor(Type15, method = "spearman")

posass15 <- (sum(sp.pair15>0)-ncol(sp.pair15))/2
negass15<- (sum(sp.pair15<0))/2
PNratio15 <- posass15/negass15

sp.pair16 <- cor(Type16, method = "spearman")

posass16 <- (sum(sp.pair16>0)-ncol(sp.pair16))/2
negass16 <- (sum(sp.pair16<0))/2
PNratio16 <- posass16/negass16

assoc.type1=c( A1$S, A1$N, A1$var.ratio, posass1, negass1, PNratio1)
assoc.type2=c( A2$S, A2$N, A2$var.ratio,  posass2, negass2, PNratio2)
assoc.type3=c( A3$S, A3$N, A3$var.ratio, posass3, negass3, PNratio3)
assoc.type4=c( A4$S, A4$N, A4$var.ratio, posass4, negass4, PNratio4)
assoc.type5=c( A5$S, A5$N, A5$var.ratio, posass5, negass5, PNratio5)
assoc.type6=c( A6$S, A6$N, A6$var.ratio, posass6, negass6, PNratio6)
assoc.type7=c( A7$S, A7$N, A7$var.ratio, posass7, negass7, PNratio7)
assoc.type8=c( A8$S, A8$N, A8$var.ratio, posass8, negass8, PNratio8)
assoc.type9=c( A9$S, A9$N, A9$var.ratio, posass9, negass9, PNratio9)
assoc.type10=c( A10$S, A10$N, A10$var.ratio, posass10, negass10, PNratio10)
assoc.type11=c( A11$S, A11$N, A11$var.ratio, posass11, negass11, PNratio11)
assoc.type12=c( A12$S, A12$N, A12$var.ratio, posass12, negass12, PNratio12)
assoc.type13=c( A13$S, A13$N, A13$var.ratio, posass13, negass13, PNratio13)
assoc.type14=c( A14$S, A14$N, A14$var.ratio, posass14, negass14, PNratio14)
assoc.type15=c( A15$S, A15$N, A15$var.ratio, posass15, negass15, PNratio15)
assoc.type16=c( A16$S, A16$N, A16$var.ratio, posass16, negass16, PNratio16)

colname <- c("S", "N", "VR", "Pos_ass", "Neg_ass", "PNratio")
rowname <- c("I", "II", "III", "IV", "V", "VI", "VII", 
"VIII","IX","X","XI","XII","XIII","XIV","XV","XVI")

Association <- as.array(rbind(assoc.type1,assoc.type2, assoc.type3, assoc.type4, 
assoc.type5, assoc.type6, assoc.type7, assoc.type8, assoc.type9, assoc.type10, 
assoc.type11, assoc.type12, assoc.type13, assoc.type14, assoc.type15, assoc.type16))

colnames(Association) <- colname
rownames(Association) <- rowname

Association

#标准化指标数据
Ass.norm <- decostand(Association,"normalize",2)
##查看指标间相关性以及进行相关性检验
#需要用到psych程序包，若无则加载命令install.packages("psych")
library(psych)
Ass.cor <- corr.test(Ass.norm)

#将结果存入csv文件

write.csv(Association, file="Association_Analysis.csv")

write.csv(Ass.cor$r, file="Association_correlation.csv")

write.csv(Ass.cor$p, file="Pvalue_Association_correlation.csv")
