args=(commandArgs(TRUE))


x<-read.delim(args[1], row.names="Symbol", sep="\t")

if (args[2] == "human") {
	annotation<-read.delim("/hihg/ref/gtf/forgene/gencode.forR.txt", sep="\t")
} else {
	annotation<-read.delim("/hihg/ref/gtf/forgene/mm10.forR.txt", sep="\t")
}

if (args[3] == "ALL" || args[3] == "edgeR") {
# edgeR
library(edgeR)
z<-x[-1,]
groups<-head(x,1)
groups_matrix<-as.matrix(groups)
group<-factor(groups_matrix)
y <- DGEList(counts=z,group=group)
y <- calcNormFactors(y)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
et <- exactTest(y)
detags <- rownames(topTags(et, n=22000))
results <- topTags(et, n=60000)
ymerged_final<-merge(z, y$pseudo.counts, by.x=0, by.y=0, all.x=T, all.y=T)
colnames(ymerged_final)[1]<-"temp"
ymerged2=merge(ymerged_final, results, by.x="temp", by.y=0, all.x=T, all.y=T)
ymerged3=merge(annotation, ymerged2, by.x="Ensembl", by.y="temp", all.x=T, all.y=T)
write.csv(file="EdgeR_results.csv", ymerged3)
}

if (args[3] == "ALL" || args[3] == "DESeq") {
# DESEQ
library(DESeq)
df<-x[-1,]
groups<-head(x,1)
groups_matrix<-as.matrix(groups)
gmat<-replace(groups_matrix, groups_matrix==0, "untreated")
gmat<-replace(gmat, gmat==1, "treated")
ltype<-replace(gmat, gmat=="treated", "paired-end")
ltype<-replace(ltype, ltype=="untreated", "paired-end")
pasillaDesign=data.frame(row.names=colnames(df), condition = as.vector(gmat), libType = as.vector(ltype))
pairedSamples = pasillaDesign$libType == "paired-end"
countTable = df[ , pairedSamples ]
condition = pasillaDesign$condition[ pairedSamples ]
condition = factor(as.vector(gmat))
cds = newCountDataSet( countTable, condition )
cds = estimateSizeFactors( cds )
cds = estimateDispersions( cds )
res = nbinomTest( cds, "untreated", "treated" )
ymerged = merge(z, normalizedCounts, by.x = 0, by.y = 0, all.x=T, all.y=T)
colnames(ymerged_final)[1]<-"temp"
yfinal = merge(annotation, ymerged, by.x="Ensembl", by.y="temp", all.x=T, all.y=T)
yexport = merge(yfinal, res, by.x="Ensembl", by.y="id", all.x=T, all.y=T)
write.csv(file="DESeq_results.txt", yexport)
}

if (args[3] == "ALL" || args[3] == "baySeq") {
library(baySeq)
library(snow)
x = read.delim("../overall.txt", header=TRUE, sep="\t")
all<-x[-1,]
g2<-head(x,1)
g2<-g2[-1]
g2_mat<-as.matrix(g2)
colnames(g2_mat)<-NULL
grp<-as.vector(g2_mat)
grp1=replace(grp, grp==1, 2)
grp2=replace(grp1, grp1==0, 1)
rep1=replace(grp2, grp2==2, 1)
replicates=list(NDE=rep1, DE=grp2)
cname <- all[,1]
all <- all[,-1]
all = as.matrix(all)
CD<-new("countData", data=all, replicates=grp2, groups=replicates)
libsizes(CD) <- getLibsizes(CD)
CD@annotation <- as.data.frame(cname)
cl <- makeCluster(8)
CDP.NBML <- getPriors.NB(CD, samplesize = 1000, estimation = "QL", cl = cl)
CDPost.NBML <- getLikelihoods.NB(CDP.NBML, pET = 'BIC', cl = cl)
topCounts(CDPost.NBML, group=2)
NBML.TPs <- getTPs(CDPost.NBML, group=2, TPs = 1:100)
bayseq_de = topCounts(CDPost.NBML, group=2, number=60000)
ymerged=merge(annotation,bayseq_de, by.x="Ensembl", by.y="annotation", all.x=T, all.y=T)
write.csv(ymerged, file="BaySeq_results.txt")
}

if (args[3] == "ALL") {
library(samr)
#some of the data will be taken from edgeR
z <- x[-1,]
groups <- head(x, 1)
groups_matrix <- as.matrix(groups)
group <- factor(groups_matrix)
y <- DGEList(counts = z, group = group)
y <- calcNormFactors(y)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
et <- exactTest(y)
detags <- rownames(topTags(et, n = 22000))
results <- topTags(et, n = 60000)
# group data will be taken from bayseq
all<-x[-1,]
g2<-head(x,1)
g2<-g2[-1]
g2_mat<-as.matrix(g2)
colnames(g2_mat)<-NULL
grp<-as.vector(g2_mat)
grp1=replace(grp, grp==1, 2)
grp2=replace(grp1, grp1==0, 1)
matrix_u<-y$pseudo.counts
matrix_rounded=round(matrix_u,0)
samfit<-SAMseq(matrix_rounded, grp2, resp.type="Two class unpaired")
ymerged=merge(annotation, samfit, by.x="Ensembl", by.y="Gene Name", all.x=T, all.y=T)
write.csv(ymerged, file="SAMSeq_results.txt")
}

if (args[3] == "ALL") {
library(EBSeq)
Sizes=MedianNorm(x)
myfactor=as.factor(grp2)
EBOut=EBTest(Data=as.matrix(x), Conditions=myfactor, sizeFactors=Sizes, maxround=5)
eb_matrix<-EBOut$PPMat
ymerged=merge(annotation, eb_matrix, by.x="Ensembl", by.y=0, all.x=T, all.y=T)
write.csv(ymerged, file="EBSeq_results.txt")
}

if (args[3] == "ALL" || args[3] == "NOISeq") {
library(NOISeq)
z<-x[-1,]
mynames<-as.vector(as.matrix(read.delim("/hihg/smoke/dvanbooven/rna_comprehensive/noiseq/ensembl.txt",header=F)))
mygc<-as.vector(as.matrix(read.delim("/hihg/smoke/dvanbooven/rna_comprehensive/noiseq/gc.txt",header=F)))
mybiotype<-as.vector(as.matrix(read.delim("/hihg/smoke/dvanbooven/rna_comprehensive/noiseq/biotype.txt", header=F)))
mychrom<-read.delim("/hihg/smoke/dvanbooven/rna_comprehensive/noiseq/location.txt", row.names="Ensembl", sep="\t")
mylength<-as.vector(as.matrix(read.delim("/hihg/smoke/dvanbooven/rna_comprehensive/noiseq/length.txt", header=F)))
myfactors<-data.frame(Status=grp2)
mydata<-readData(data=z, length=mylength, gc=mygc, biotype=mybiotype, chromosome=mychrom, factors=myfactors)
QCreport(mydata, samples = NULL, factor = "Status")
mynoiseq=noiseq(mydata, k = 0.5, norm="rpkm", factor="Status", conditions=c(0,1), pnr=0.2, nss=5, v=0.02, lc=1, replicates="no")
ymerged=merge(annotation, mynoiseq@results[[1]], by.x="Ensembl", by.y=0, all.x=T, all.y=T)
yfinal=ymerged[1:10]
write.csv(file="NOISeq_Results.txt", yfinal)
}

#Remaining differential callers to add are DEGSeq and PoissonSeq.  Then all is required to merge all of the various pieces into a single annotated output file
