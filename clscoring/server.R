# server logic
server <- function(input, output, session){ #clientData, 
  
  if(!exists("a")){x <- 100}
  
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
                                if(is.null(vals$rf) | input$tum=="BRCA"){ #because I don't have it for brca yet
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
  
  #cell line distance
  cl_dist <- eventReactive(input$do_score,
                           {
                             #paste("Testing CL distance text")
                           })
  
  
  
  output$mds <- renderPlot(vals$rfmds)
  output$pred <- renderText(cancer_pred)
  output$cl_dist <- renderText(cl_dist())
  output$rez <- renderText(rez_cl_tum())
  output$feele <- renderTable({
    req(input$fajla)
    df <- read.csv(input$fajla$datapath)
    return(head(df))
  })
  
}
