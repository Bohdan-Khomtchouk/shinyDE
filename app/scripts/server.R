library(shiny)
library(edgeR)
library(VennDiagram)
library(DESeq)
library(baySeq)
library(snow)
library(gplots)
library(splitstackshape)
library(goseq)
library("GenomeInfoDb")

edgerDE<-c(1,2,3)
edgerFiltered<-c(1,2,3)
edgerOverall<-c(1,2,3)
deseqDE<-c(1,2,3)
bayseqDE<-c(1,2,3)
intersectionDE<-c(1,2,3)
deseqOverall<-c(1,2,3)
bayseqOverall<-c(1,2,3)
n12=0
n13=0
n23=0
n123=0

shinyServer(function(input, output, session) {

  download.file("https://www.dropbox.com/s/ew8h9e75f9qfczr/gencode.forR.txt?dl=0", destfile="./temp1.gtf", method="wget")
  annotationH <- read.delim("temp.gtf", sep="\t")

  
  
  output$mytable1 <- renderDataTable({
    
  	# input$file1 will be NULL initially. After the user selects
    # and uploads a file, it will be a data frame with 'name',
    # 'size', 'type', and 'datapath' columns. The 'datapath'
    # column will contain the local filenames where the data can
    # be found.

    
    inFile <- input$file1
    if (is.null(inFile))
      return(NULL)
    
    cat(input$dist)

  if (input$dist == "ALL" || input$dist == "edgeR")  {
    withProgress(message='Running edgeR', detail = 'Might take a min', value=NULL, {
    
    x <- read.delim(inFile$datapath, row.names="Symbol")
    
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
	as.data.frame(ymerged_final)
	colnames(ymerged_final)[1]<-"temp"
	ymerged2=merge(ymerged_final, results, by.x='temp', by.y=0, all.x=T, all.y=T)
	ymerged3=merge(annotationH, ymerged2, by.x='Ensembl', by.y='temp', all.x=T, all.y=T)
	edgerOverall<<-ymerged3
  edgerDE<<-subset(ymerged3, FDR<0.05)
  gene_list<-as.vector(edgerDE$Ensembl)
  m<-y$pseudo.counts
  edgerFiltered<<-m[gene_list,]
  as.data.frame(ymerged3)
  })
	}
	    
  }, options = list(aLengthMenu=c(5,30,50), iDisplayLength=100))
  
  
  output$mytable2 <- renderDataTable({
  cat("whoa nelly")
  edgerDE
  inFile <- input$file1
  if (is.null(inFile))
    return(NULL)
  
  if (input$dist == "ALL" || input$dist == "DESeq")  {
  x <- read.delim(inFile$datapath, row.names="Symbol")
  
  withProgress(message='Running DESeq', detail = 'Might take a min or 2', value=NULL, { 
  
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
  resmerged=merge(annotationH, res, by.x="Ensembl", by.y="id", all.x=T, all.y=T)
  as.data.frame(resmerged)
  deseqDE<<-subset(resmerged, padj<0.05)
  
  })}
  }, options = list(aLengthMenu=c(5,30,50), iDisplayLength=100))
  
  
  output$mytable3 <- renderDataTable({
    
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, it will be a data frame with 'name',
    # 'size', 'type', and 'datapath' columns. The 'datapath'
    # column will contain the local filenames where the data can
    # be found.
    
    
    inFile <- input$file1
    if (is.null(inFile))
      return(NULL)
    
    cat(input$dist)
    
    if (input$dist == "ALL" || input$dist == "baySeq")  {
      
      x <- read.delim(inFile$datapath, header=TRUE, sep="\t")
      
      withProgress(message='Running baySeq', detail = 'Please be patient, this will take several minutes.', value=NULL, {
        
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
        CDPost.NBML <- getLikelihoods(CDP.NBML, pET = 'BIC', cl = cl)
        topCounts(CDPost.NBML, group=2)
        NBML.TPs <- getTPs(CDPost.NBML, group=2, TPs = 1:100)
        bayseq_de = topCounts(CDPost.NBML, group=2, number=60000)
        baymerged=merge(annotationH,bayseq_de, by.x="Ensembl", by.y="annotation", all.x=T, all.y=T)
        as.data.frame(baymerged)
        bayseqDE<<-subset(baymerged, FDR.DE<0.05)
        
      })
    }
    
  }, options = list(aLengthMenu=c(5,30,50), iDisplayLength=100))
  
  output$venn_it <- renderPlot ({
    
    if (input$dist == "ALL")  {
    # calculate the comparisons between differential lists.
      
    
    draw.triple.venn(area1=nrow(edgerDE), area2=nrow(deseqDE), area3=nrow(bayseqDE), n12=n12, n23 = n23, n13 = n13, n123 = n123, category=c("edgeR", "DESeq", "baySeq"))
      
    } 
    
  }, width=750, height=666)
  
  output$mytable4 <- renderDataTable ({
    
    merge1<-merge(edgerDE, deseqDE, by.x="Ensembl", by.y="Ensembl")
    merge2<-merge(deseqDE, bayseqDE, by.x="Ensembl", by.y="Ensembl")
    merge3<-merge(edgerDE, bayseqDE, by.x="Ensembl", by.y="Ensembl")
    intersectionDE<<-merge(merge1, bayseqDE, by.x="Ensembl", by.y="Ensembl")
    
    n12 <<- nrow(merge1)
    n23 <<- nrow(merge2)
    n13 <<- nrow(merge3)
    n123 <<- nrow(intersectionDE)
    
    as.data.frame(intersectionDE)
    
  })
  
  output$heatmap <- renderPlot ({
    
    withProgress(message='Drawing Heatmap', detail = 'min or so', value=NULL, {
    heatmap.2(edgerFiltered)
    })
  })
  
  output$ontologies <- renderDataTable ({
    withProgress(message='Calculating Ontology Enrichment', detail = 'min or so', value=NULL, {
    genes=as.integer(p.adjust(edgerOverall$FDR[edgerOverall$logFC!=0]<0.05))
    ensembl_for_trim=edgerOverall$Ensembl[edgerOverall$logFC!=0]
    splitEnsembl=cSplit(as.data.frame(ensembl_for_trim), "ensembl_for_trim", sep=".", type.convert=FALSE)
    names(genes)=splitEnsembl$ensembl_for_trim_1
    pwf=nullp(genes,"hg19","ensGene")
    GO.wall=goseq(pwf,"hg19", "ensGene")
    GO.wall$FDR<-p.adjust(GO.wall$over_represented_pvalue, method="BH")
    GO_overall=subset(GO.wall, FDR<0.05)
    GO_overall
    
    })
  })
  
  output$downloadData <- downloadHandler(
    
    
    filename = "rna_results.csv", 
    cat("whoa"),
    content=function(file) { 
      
      write.table(edgerOverall, file, row.names=TRUE);
    }
  )
  
})