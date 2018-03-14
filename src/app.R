#cl_score_app/app.R
#needs new UI for the score of cl-tum

library(shiny)

library(randomForest)
library(reprtree)
library(ROCR)
library(knitr)
library(maftools)
source("C:/Users/Vladimir/Desktop/tcga/rf_functions_toimport.R")

load("C:/Users/Vladimir/Desktop/tcga/for_app.Rdata")


names_cl_endo <- as.vector(sapply(rownames(CN_endometrial_CCLE),
                                  function(x) substr(x, 1, nchar(x) - nchar("_endometrium"))))
rownames(CN_endometrial_CCLE) <- names_cl_endo


# Define UI for application that draws a histogram
ui <- fluidPage(
  titlePanel("Cell line scoring"),

  sidebarLayout(
    sidebarPanel("Select tumour type and cell line:",
                 selectInput("tum", label = "Tumour type:", 
                             list("Tumour cohorts" = c("UCEC", "BRCA"))),
                 actionButton("do_model", "Build model"),
                 hr(),
                 selectInput("cl", label = "Pick a CL of x type",
                             choices = c("AN3CA", "MFE-280", "JHUEM-3") ),
                 actionButton("do_score", "Score the cell line"),
                 hr(), hr(),
                 fileInput("fajla",
                           "Upload own TCGA subtypes for the selected type",
                           accept = c("text/csv", 
                                      "text/comma-separated-values,text/plain", ".csv")),
                 actionButton("barcode_update","Update the barcode classification")
                 ),
    mainPanel("Tumour classifier performance",
              plotOutput("mds"),
              hr(),
              textOutput({"rez"}),
              tableOutput("feele"))
  )
)





# server logic
server <- function(input, output, clientData, session) {
  
  #changing cell line selection based on tumour type selected
  observe({
      c_tum <- input$tum
      cls <- switch(c_tum,
             UCEC = names_cl_endo, #c("AN3CA", "MFE"),
             BRCA = c("MCF7", "CAL51"))
      updateSelectInput(session, "cl",
                        label = paste("Pick a cell line of", input$tum,"type:"),
                        choices = cls)
    })
  
  
  vals <- reactiveValues(rf = NULL, rfmds = hist(c(1,2,3)))
  
  #updating TCGA barcode classification
  observeEvent(input$barcode_update,
               {
                 req(input$fajla)
                 TCGA_histology <<- read.csv(input$fajla$datapath)
               })
  
  #building rf tumour model and generating the plot
  observeEvent(input$do_model,
                               { if(input$tum=="UCEC"){
                                 environment(rf_default) <- environment() #isuse R je odvratan, a i ja
                                 vals$rf <<- rf_default()
                                 dev.off()
                                 MDSplot(vals$rf, training_set$class, palette=c(3,2))
                                 vals$rfmds <<- recordPlot()
                               }else{vals$rf <<- NULL}
                                
                               }
                              )

  #scoring cell line
  rez_cl_tum <- eventReactive(input$do_score,
                           {
                            if(is.null(vals$rf) | input$tum=="BRCA"){ #because I don't have it for brca
                              paste(input$tum, input$cl)
                            }else{
                              CN_endometrial_CCLE <- rbind(training_set[1, 1:ncol(CN_endometrial_CCLE)],
                                                           CN_endometrial_CCLE)
                              CN_endometrial_CCLE <- CN_endometrial_CCLE[-1, ]
                              paste("Endometrioid subtype score of", input$cl, ":",
                                    predict(vals$rf, CN_endometrial_CCLE[input$cl,], type="prob")[1])
                            }
                                   
                          }
                        )

  
    
  output$mds <- renderPlot(vals$rfmds)
  output$rez <- renderText(rez_cl_tum())
  output$feele <- renderTable({
    req(input$fajla)
    df <- read.csv(input$fajla$datapath)
    return(head(df))
  })
  
}


# Running
shinyApp(ui = ui, server = server)
