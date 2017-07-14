args=(commandArgs(TRUE))


x<-read.delim(args[1], row.names="Symbol", sep="\t")


if (args[2] == "human") {
	annotation<-read.delim("/hihg/ref/gtf/forgene/gencode.forR.txt", sep="\t")
	# Alternative method to download from dropbox for outside of pegasus usage.
	#download.file("https://www.dropbox.com/s/ew8h9e75f9qfczr/gencode.forR.txt?dl=0", destfile="./temp.gtf", method="wget")
	#annotation<-read.delim("temp.gtf", sep="\t")
} else {
	annotation<-read.delim("/hihg/ref/gtf/forgene/mm10.forR.txt", sep="\t")
}

if (args[3] == "PoissonSeq" || args[3] == "ALL") {

	library(PoissonSeq)
	hdat<-read.table(args[1], header=T, stringsAsFactors=T)
	grp1<-head(hdat,1)
	grp2<-grp1[,-1]
	z<-hdat[-1,]
	gname<-z[,1]
	n<-as.matrix(z[,-1])
	type="twoclass"
	pair<-TRUE
	y<-grp2+1
	y2<-unlist(y)
	dat<-list(n=n, y=y2, type=type, pair=pair, gname=gname)
	res<-PS.Main(dat=dat)	
	write.csv(file="PoissonSeq_results.csv", res)
}

