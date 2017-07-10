library(shiny)

shinyUI(pageWithSidebar(

  headerPanel("shinyDE"),

  sidebarPanel(
    tags$div(),


    #radioButtons("dist", "DE Method : ", c("edgeR"="edgeR", "DESeq"="DESeq", "baySeq"="baySeq", "ALL"="ALL")),
    checkboxGroupInput("differential_callers", label=h4("DE Callers:"), choices = list("edgeR" = "edgeR", "DESeq2" = "DESeq2", "baySeq" = "baySeq", "NOISeq" = "NOISeq", "SAMSeq" = "SAMSeq", "DEGseq" = "DEGSeq", "EBSeq" = "EBSeq", "PoissonSeq" = "PoissonSeq")),

    fileInput('file1', 'Choose matrix count file to upload:',
                accept = c(
                  'text/csv',
                  'text/comma-separated-values',
                  'text/tab-separated-values',
                  'text/plain',
                  '.csv',
                  '.tsv'
                )
      ),


    actionButton("do", "Run Query"),
    tags$br(),
    downloadButton("downloadData", "Download Results"),
    tags$br(),
    tags$br(),
    img(src = "UM.png", width="400px", height="75px")

  ),

  mainPanel(
  	tabsetPanel(type = "tabs", tabPanel("Instructions", verbatimTextOutput("instr")), tabPanel("Differential Expression Table", value="mytable1", dataTableOutput("mytable1")), tabPanel("Heatmap", plotOutput("GTEx_Heatmap"), height="auto", width="auto"), tabPanel("VennDiagram", plotOutput("venn"), height="auto", width="auto"), tabPanel("Ontologies", dataTableOutput("mytable5")))
    #tabsetPanel(type = "tabs", tabPanel("GTEx_Data Table", dataTableOutput("mytable1")), tabPanel("DESeq", dataTableOutput("mytable2")), tabPanel("baySeq", dataTableOutput("mytable3")), tabPanel("Intersecting Table", dataTableOutput("mytable4")), tabPanel("Venn", plotOutput("venn_it")), tabPanel("Heatmap",plotOutput("GTEx_Heatmap")))
  )

))
