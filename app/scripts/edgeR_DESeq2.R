args=(commandArgs(TRUE))


x<-read.delim(args[1], row.names="Symbol", sep="\t")

if (args[2] == "human") {
	annotation<-read.delim("/hihg/ref/gtf/forgene/gencode.forR.txt", sep="\t")
} else if (args[2] == "mm10") {
	annotation<-read.delim("/hihg/ref/gtf/forgene/mm10.forR.txt", sep="\t")
} else if (args[2] == "hg38_v25") {
	annotation<-read.delim("/hihg/ref/gtf/forgene/gencode.v25.forR.txt", sep="\t")
} else if (args[2] == "hg38_v24") {
	annotation<-read.delim("/hihg/ref/gtf/forgene/gencode.v24.forR.txt", sep="\t")
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

if (args[3] == "ALL" || args[3] == "DESeq2") {
# DESEQ2
library(DESeq2)
df<-x[-1,]
groups<-head(x,1)
groups_matrix<-as.matrix(groups)
gmat<-replace(groups_matrix, groups_matrix==0, "untreated")
gmat<-replace(gmat, gmat==1, "treated")
ltype<-replace(gmat, gmat=="treated", "paired-end")
ltype<-replace(ltype, ltype=="untreated", "paired-end")
pasillaDesign=data.frame(row.names=colnames(df), condition = as.vector(gmat), lib = as.vector(ltype))
pairedSamples = pasillaDesign$libType == "paired-end"
countTable = df[ , pairedSamples ]
condition = pasillaDesign$condition[ pairedSamples ]
condition = factor(as.vector(gmat))
dds<-DESeqDataSetFromMatrix(countData=df, colData=pasillaDesign, design = ~condition)
dds<-DESeq(dds)
res = results(dds)
normalizedCounts <- counts(dds, normalized=TRUE)
ymerged = merge(df, normalizedCounts, by.x = 0, by.y = 0, all.x=T, all.y=T)
colnames(ymerged)[1]<-"temp"
yfinal = merge(annotation, ymerged, by.x="Ensembl", by.y="temp", all.x=T, all.y=T)
res2<-as.data.frame(res)
yexport = merge(yfinal, res2, by.x="Ensembl", by.y=0, all.x=T, all.y=T)
write.csv(file="DESeq2_results.txt", yexport)
}

