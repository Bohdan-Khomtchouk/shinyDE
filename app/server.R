library(shiny)
library(gplots)
library(dplyr)
library(edgeR)
library(DESeq2)
library(EBSeq)
library(baySeq)
library(PoissonSeq)
library(DEGseq)

gtex_table <- c(1,1,1)
bs_table <- c(0)
final_tbl <- c(0)
annot_tbl <- c(0)
annot_tbl2 <- c(0)
globali = 0
annotationH <- c(0)

shinyServer(function(input, output, session) {

  datasetInput <- reactive({

      infile <- input$file1

      cat(file = stderr(), infile$datapath)

      if (is.null(infile$datapath)) {
        showModal(modalDialog(
          title = "Important!",
          "You need to select a valid file.",
          easyClose = TRUE
        ))
        return()
      }
      #x <- read.table(infile$datapath, header = F, col.names = c("GeneSymbol"))
      x <- read.delim(infile$datapath, row.names = "Symbol")

      annotationH <<- read.delim("temp.gtf", sep = "\t")

    #else
    #{
    #  x <- t(data.frame(do.call(rbind, strsplit(input$caption, "\n", fixed = TRUE))))
    #  colnames(x) <- "GeneSymbol"
    #  x <- mutate_each(as.data.frame(x), funs(toupper))
    #}
    as.data.frame(x)

  })
  observeEvent(input$do, {

    x <- datasetInput()
    cat(file = stderr(), length(x))
    if (length(x) == 0) {
      showModal(modalDialog(
        title = "Important!",
        "You need to select a valid file.",
        easyClose = TRUE
      ))
      return()
    }
    cat(file = stderr(), length(input$differential_callers))
    cat(file = stderr(), input$differential_callers)
    if (length(input$differential_callers) == 0) {
      showModal(modalDialog(
        title = "Important!",
        "You need to select a differential expression caller",
        easyClose = TRUE
      ))
      return()
    }
    # check if edgeR is in the list if so then run edgeR
    if (length(input$differential_callers) > 4) {
      showModal(modalDialog(
        title = "Important!",
        "You selected more than 4 differential expression callers.  Please revise your list.  In future versions of
        shinyDE, you will be able to run any number of DE callers because the program will be parallelized for speed.
        For now, we apologize for the inconvenience!",
        easyClose = TRUE
      ))
      return()
    }
    if ('edgeR' %in% input$differential_callers) {
    withProgress(message = "Running EdgeR", value = NULL, {
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
      ymerged_final <- merge(z, y$pseudo.counts, by.x = 0, by.y = 0, all.x = T, all.y = T)
      as.data.frame(ymerged_final)
      colnames(ymerged_final)[1] <- "temp"
      ymerged2 = merge(ymerged_final, results, by.x = 'temp', by.y = 0, all.x = T, all.y = T)
      ymerged3 = merge(annotationH, ymerged2, by.x = 'Ensembl', by.y = 'temp', all.x = T, all.y = T)
      edgerOverall <<- ymerged3
      edgerDE <<- subset(ymerged3, FDR < 0.05)
      gene_list <- as.vector(edgerDE$Ensembl)
      m <- y$pseudo.counts
      edgerFiltered <<- m[gene_list,]
      gtex_table <<- ymerged3
    })}

    if ('DESeq2' %in% input$differential_callers) {
    withProgress(message = "Running DESeq2", value = NULL, {
      df <- x[-1,]
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
      #gtex_table <<- yexport
    })}

    if ('EBSeq' %in% input$differential_callers) {
    withProgress(message = "Running EBSeq", value = NULL, {
      Sizes=MedianNorm(x)
      all <- x[-1,]
      g2 <- head(x, 1)
      g2 <- g2[-1]
      g2_mat <- as.matrix(g2)
      colnames(g2_mat) <- NULL
      grp <- as.vector(g2_mat)
      grp1 = replace(grp, grp == 1, 2)
      grp2 = replace(grp1, grp1 == 0, 1)
      myfactor=as.factor(grp2)
      EBOut=EBTest(Data=as.matrix(x), Conditions=myfactor, sizeFactors=Sizes, maxround=5)
      eb_matrix<-EBOut$PPMat
      ymerged=merge(annotation, eb_matrix, by.x="Ensembl", by.y=0, all.x=T, all.y=T)
    })}

    if ('baySeq' %in% input$differential_callers) {
    withProgress(message = "Running baySeq", value = NULL, {
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
    })}

    if ('PoissonSeq' %in% input$differential_callers) {
    withProgress(message = "Running PoissonSeq", value = NULL, {
      grp1 <- head(x, 1)
      grp2 <- grp1[,-1]
      z <- x[-1,]
      gname <- z[,1]
      n <- as.matrix(z[,-1])
      type = "twoclass"
      pair <- TRUE
      y <- grp2 + 1
      y2 <- unlist(y)
      dat <- list(n = n, y = y2, type = type, pair = pair, gname = gname)
      res <- PS.Main(dat = dat)
    })}

    if ('DEGseq' %in% input$differential_callers) {
    withProgress(message = "Running DEGseq", value = NULL, {
      count_matrix2 <- x[-1,]
      g1 <- head(x, 1)
      j = 1
      i = length(g1)
      c1 <- c()
      c2 <- c()
      while (j <= i) {
        if (g1[j] == 0) { c1 <- c(c1, j) }
        if (g1[j] == 1) { c2 <- c(c2, j) }
        j = j + 1
      }
      DEGexp(geneExpMatrix1 = count_matrix2, geneCol1 = 1, expCol1 = c1, groupLabel1 = "group1", geneExpMatrix2 = count_matrix2, geneCol2 = 1, expCol2 = c2, groupLabel2 = "group2", method = "MARS", outputDir = "./")
    })}


  })

  GT <- reactive({
    if (input$do) {
      as.data.frame(gtex_table)
    }
    else {
      gtex_table <<- c(2)
    }
  })

  BT <- reactive({

    if (input$do) {
      as.data.frame(bs_table)
    }

  })

  AT <- reactive({
    if (input$do) {
      as.data.frame(annot_tbl)
    }
  })
  output$mytable1 <- renderDataTable({
    if (input$do) {
      as.data.frame(GT())
    }

    #if (globali == 0) {
    #if (input$do) {
    #  cat("upon rendering\n")
    #  if(input$caption == "") {
    #    infile <- input$file1
    #    x <- read.table(infile$datapath, header = F, col.names = c("GeneSymbol"))
    #  }
    #  else
    #  {
    #    x <- t(data.frame(do.call(rbind, strsplit(input$caption, "\n", fixed = TRUE))))
    #    colnames(x) <- "GeneSymbol"
    #  }
  #
  #  z <- read.table("www/GTEx.txt", header = T, sep = "\t")
  #  y <- z
  #  if (input$gtex_filter == "norm") {
  #    y <- z[1:15]
  #  }
  #  final_tbl <- merge(x, y, by.x = "GeneSymbol", by.y = "Description")
  #
  #  gtex_table <<- final_tbl
  #  as.data.frame(gtex_table)
#    }}

  })

  output$GTEx_Heatmap <- renderPlot({

    if (globali == 0)
    {
      return(NULL)
    }

    gtex2 <- GT()
    rnames <- gtex2[,1]
    z <- data.matrix(gtex2[,3:ncol(gtex2)])
    z[z==0]<-0.00001
    zlog <- log(z,10)
    heatmap.2(zlog, labRow = rnames, tracecol = NA, col = greenred(10), margins = c(16,8), cexCol = 0.9)

    #infile <- input$file1
    #if (is.null(infile))
    #    return(NULL)
    #rnames <- gtex_table[,1]
    #z <- data.matrix(gtex_table[,3:ncol(gtex_table)])
    #z[z==0] <- 0.00001
    #zlog <- log(z,10)
    #heatmap.2(zlog, labRow = rnames, tracecol = NA, col = greenred(25), margins = c(16,8), cexCol = 0.9)


  }, height = 1000)

  output$mytable3 <- renderDataTable({
    as.data.frame(BT())
    #infile <- input$file1
    #if (is.null(infile))
    #  return(NULL)
    #withProgress(message = "Running BrainSpan", value = NULL, {
      #cat(input$bs_filter)
      #cat(input$bs_gender)

      #type1 <- switch(input$bs_filter, opt1 = "opt1", opt2 = "opt2", opt3 = "thirdoption")
      #file_to_open = paste(type1, input$bs_gender, sep = ".")
      #openfile = paste("www/", file_to_open, sep = "")

      #x <- read.table(infile$datapath, header = F, col.names = c("GeneSymbol"))
      #y <- read.table(openfile, header = T, sep = "\t")
      #final_tbl <- merge(x, y, by.x = "GeneSymbol", by.y = "GeneSymbol")
      #as.data.frame(final_tbl)
      #bs_table <<- final_tbl

    #})

  })


  output$BrainSpan_Heatmap <- renderPlot({
    bstable2 <- BT()

    if (globali == 0)
    {
      return(NULL)
    }
    rnames <- bstable2[,1]
    z <- data.matrix(bstable2[,3:ncol(bstable2)])
    z[z==0] <- 0.00001
    zlog <- log(z,10)
    heatmap.2(zlog, labRow = rnames, tracecol = NA, col = greenred(10), dendrogram = "row", Colv = F, margins = c(16,8), cexCol = 0.9)

  }, height = 1000)
  output$mytable5 <- renderDataTable({

    as.data.frame(AT())
    #infile <- input$file1
    #if (is.null(infile))
    #  return(NULL)
    #x <- read.table(infile$datapath, header = F, col.names = c("GeneSymbol"))
    #y <- read.table("www/exac_intolerance.txt", header = T, sep = "\t")
    #final_tbl <- merge(x, y, by.x = "GeneSymbol", by.y = "gene")
    #as.data.frame(final_tbl)
    #annot_tbl <<- final_tbl
  }, escape = FALSE)

  output$downloadData <- downloadHandler(
    filename = function() { paste("complete_results", ".csv", sep = "") },
    content = function(file) {
      a <- merge(gtex_table, bs_table, by.x = "GeneSymbol", by.y = "GeneSymbol")
      b <- merge(annot_tbl2, a, by.x = "GeneSymbol", by.y = "GeneSymbol")
      write.csv(b, file)
    }
  )

  output$mytable6 <- renderDataTable({

  })

  output$instr <- renderText({
    "
    The purpose of this site is to calculate differential expression (DE) across multiple R packages and return an overlap list between multiple DE callers.

    Just select from the list of differential expression callers on the left panel and upload a matrix file of counts.  The format of the count file is very specific.  The first line should include the sample IDs, the second line should include the group designation and the first column should be a list of Ensembl gene identifiers.

    For example a typical input would look something like this :

    Symbol	Sample1	Sample2	Sample3	Sample4	Sample5	Sample6
    Group	0	0	0	1	1	1
    ENSG000000001	2	4	3	10	32	15
    .
    .
    .

    Once the file is uploaded then click on 'Run Query' to run the differential expression.  Once the analysis has completed you will see a tabset panel with the results.

    Some notes on the results panel :
    - The heatmap is only generated from the list of genes that overlapped between all callers.  If you only ran 1 DE caller, then it's all the differentially expressed genes of that particular caller.  If you ran 4 callers then the heatmap is the intersection of the 4 callers.
    - Ontologies work the same way as the heatmap.  Only the intersection of genes is taken.  Also, note that ontology is only calculated for human, mouse, and rat genomes at this time.
    "


  })

})
