
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#


#more plots
#post it online

library(shiny)
library(MetKMR)
library(GenomicRanges)
library(rtracklayer)
library(Biobase)
library(dplyr)

# This is for allowing to upload large files

options(shiny.maxRequestSize=1000*1024^2)

# Define UI for application 

demo_data <- load('data/analysis_fa.RData')

ui <- fluidPage(
  
  # Application title
  titlePanel(title=div(img(src="image.png"))),
  
  # Sidebar with the parameters for MetRKAT
  sidebarLayout(
    sidebarPanel(
      radioButtons('demo', 'Do you want to run a Demo or to analyze data',
                   choices = list('Demo Mode' = 'demo','Analysis Mode' = 'analysis')),
      
      helpText('The outcome variable must be a table with one column'),
      
      fileInput("output", 'Insert your outcome variable table' ,multiple = FALSE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")),
      checkboxInput("outcometype", "Is your outcome variable dichotomous?", TRUE),
      
      fileInput("file1", 'Insert your table of normalized beta values in csv format' ,multiple = FALSE,
                accept = c("text/csv",
                "text/comma-separated-values,text/plain",
                ".csv")),
      # Input: Checkbox if file has header ----
      checkboxInput("header", "Header", TRUE),
      
      # Input: Select separator ----
      radioButtons("sep", "Separator",
                   choices = c(Comma = ",",
                               Semicolon = ";",
                               Tab = "\t"),
                   selected = ","),
      
      fileInput("file2", 'Insert your annotation table in csv format' ,multiple = FALSE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")),
      # Input: Checkbox if file has header ----
      checkboxInput("header2", "Header", TRUE),
      
      # Input: Select separator ----
      radioButtons("sep2", "Separator",
                   choices = c(Comma = ",",
                               Semicolon = ";",
                               Tab = "\t"),
                   selected = ","),
      
      numericInput('wsize', 'Select the window size. This number must be 
                   positive', min = 0, max = 5000, value = 1000),
      
      numericInput('gap', 'Select the gap size between sliding windows.',
                   min = 0, max = 4999, value = 0),
      
      radioButtons('distmethod', 'Select the distance computation method',
                   choices = list('Euclidean' = 'euclidean',
                                  'Manhattan' = 'manhattan',
                                  'Both' =  'omnibus')),
      radioButtons('wmethod', 'Select the window generation method',
                   choices = list('Default' = 'default',
                                  'By gene' = 'genes',
                                  'By location' = 'location')),
      

      
      numericInput('max.na', 'Select the maximun NA proportion',
                   max=1, min = 0, value = 0.3, step = 0.01),
      

      actionButton("run", "Run Analysis")
      
      
      ),
    
    # Show the results
    mainPanel(h1('Results'),
              conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                               tags$div("Running... Analysis may take a while please be patient",id="loadmessage")),
              #p('Manhattan plot:'),
              plotOutput('manhattan'),
              p('Most differencially expressed CpGs:'),
              tableOutput("head"),


              helpText(" Click on the download button to download all the differentially methylated CpGs as .txt  format"),
              downloadButton('downloadData', 'Download'),
              textInput('which_ch', 'Which chromosome would you like to plot?',
                        value = 10),
              plotOutput('chromosome')
              

    )
      )
      )



# Define server logic to render the plots and outputs

server <- function(input, output,session) {
 
  
  data <- eventReactive(input$run, ({
    
    if(input$demo == 'demo'){
      data <- readRDS('data/analysis.RDS')
      return(data)
      
      
    }else{
      #load('data/annotation2.rda')
      betas <- read.csv(input$file1$datapath,
                     header = input$header,
                     sep = input$sep)
      annotation2<-read.csv(input$file2$datapath,
                          header = input$header2,
                          sep = input$sep2)
      phenotype<- read.csv(input$output$datapath,
                           header = F,
                           sep = " ")
     # betas<-read.table('data/betas.txt',sep="\t",header=T)
      if (input$distmethod=="manhattan"){
        distmethod="manhattan"
      }
      if(input$distmethod=="euclidean"){
        distmethod="euclidean"
      }
      if(input$distmethod=="omnibus"){
        distmethod=c("euclidean","manhattan")
      }
     
       data <- new("MetRKAT",
                  data = betas,
                  annotation = annotation2,
                  distmethod =  distmethod,
                  wsize = input$wsize, gap = input$gap,
                  max.na = input$max.na,wmethod = input$wmethod)
      data@intervals <- createIntervals(data)
          if(input$outcometype==TRUE){
                data@results <- applyRKAT(data, y = as.integer(phenotype[,1]))
                return(data) }
          else{
                data@results <- applyRKAT(data, y = as.integer(phenotype[,1]),out_type ='C')
        return(data)
        
      }

    }
  }))
  
  output$head <- renderTable({
    result_filtered<-as.data.frame(data()@results)[as.data.frame(data()@results)$pval<= 0.05, ]
    DMG_symbol<-data()@annotation[result_filtered$first_row, 'gene']
    result_filtered$gene<-DMG_symbol
    site_symbol<-data()@annotation[result_filtered$first_row, 'site']
    result_filtered$site<-site_symbol
    head(result_filtered[order(result_filtered$pval), ])
  })
  
  output$chromosome <- renderPlot({
    
    plotChromosome(data(),
                   paste('chr',input$which_ch, sep = ''),
                   pvals = input$distmethod)
    
    print(paste('chr',input$which_ch))
    
    
  })
  
  output$manhattan <- renderPlot({
    
    plotManhattan(data(), pvals =input$distmethod)
    
  })
  

  
  output$downloadData <- downloadHandler(
    
    # This function returns a string which tells the client
    # browser what name to use when saving the file.
    filename = function() {
      paste("result", "txt", sep = ".") # example : iris.Rdata
      
    },
    
    # This function should write data to a file given to it by
    # the argument 'file'.
    content = function(file) {
      result_filtered<-as.data.frame(data()@results)
      result_filtered<-result_filtered[result_filtered$pval<= 0.05, ]
      DMG_symbol<-data()@annotation[result_filtered$first_row, 'gene']
      result_filtered$gene<-DMG_symbol
      site_symbol<-data()@annotation[result_filtered$first_row, 'site']
      result_filtered$site<-site_symbol
      write.table(result_filtered[order(result_filtered$pval), ], file ,sep="\t",row.names = F)
    }
  )
}


# Run the application 
shinyApp(ui = ui, server = server)
