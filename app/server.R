library(shiny)
library(gplots)
library(dplyr)
library(edgeR)

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

  #HEATMAP
  #output here
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
