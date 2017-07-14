args=(commandArgs(TRUE))


if (args[2] == "human") {
	annotation<-read.delim("/hihg/ref/gtf/forgene/gencode.forR.txt", sep="\t")
	# Alternative method to download from dropbox for outside of pegasus usage.
	#download.file("https://www.dropbox.com/s/ew8h9e75f9qfczr/gencode.forR.txt?dl=0", destfile="./temp.gtf", method="wget")
	#annotation<-read.delim("temp.gtf", sep="\t")
} else {
	annotation<-read.delim("/hihg/ref/gtf/forgene/mm10.forR.txt", sep="\t")
}

if (args[3] == "ALL" || args[3] == "DEGseq") {

	library(DEGseq)
	count_matrix <- read.delim(args[1])
	count_matrix2<-count_matrix[-1,]
	g1<-head(count_matrix,1)
	j = 1
	i = length(g1)
	c1<-c()
	c2<-c()
	while (j <= i) {
		if (g1[j] == 0) { c1 <- c(c1, j) }
		if (g1[j] == 1) { c2 <- c(c2, j) }
		j = j + 1
	}
	DEGexp(geneExpMatrix1=count_matrix2, geneCol1=1, expCol1=c1, groupLabel1="group1", geneExpMatrix2=count_matrix2, geneCol2 = 1, expCol2=c2, groupLabel2="group2", method="MARS", outputDir="./")

}

