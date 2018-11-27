# server logic
server <- function(input, output, session){ #clientData, 
  
  #changes cell line selection based on tumour type selected
  observe({
    c_tum <- input$tum
    cls <- switch(c_tum,
                  UCEC = names_cl_endo, #c("AN3CA", "MFE"),
                  BRCA = c("MCF7", "CAL51"))
    updateSelectInput(session, "cl",
                      label = paste("Select a cell line of", input$tum,"type:"),
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
                 
                 print(input$datatypeSelected)
                 
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
  
  rez_all_cl_plot <- eventReactive(input$do_score_2,
                              {
                                if(is.null(vals$rf) | input$tum=="BRCA"){
                                  paste(input$tum, input$cl)
                                }else{
                                  #black magic
                                  CN_endometrial_CCLE <- rbind(training_set[1, 1:ncol(CN_endometrial_CCLE)],
                                                               CN_endometrial_CCLE)
                                  CN_endometrial_CCLE <- CN_endometrial_CCLE[-1, ]

                                  rez_list <- sapply(rownames(CN_endometrial_CCLE),
                                         function(cl) predict(vals$rf, CN_endometrial_CCLE[cl,], type="prob")) %>% t %>%
                                    as.data.frame %>% rownames_to_column(var = "Name")
                                  colnames(rez_list) <- c("Name","Class1","Class2")

                                  rez_list


                                }

                              }
  )
  
  
  
  
  output$mds <- renderPlot(vals$rfmds)
  output$pred <- renderText(cancer_pred)
  output$rez <- renderText(rez_cl_tum())

  output$plotscore <- renderPlot({
      if(is.null(vals$rf)) return(NULL)
      rez_all_cl_plot() %>% ggplot(aes(x=reorder(Name, Class1),y=Class1, fill=Class2)) +
               geom_bar(stat = "identity", colour="black") +
               coord_flip() + scale_fill_viridis_c() +
               guides(fill=F) +
              labs(x="Name",y="Classifier score")
    })

  output$feele <- renderTable({
    req(input$fajla)
    df <- read.csv(input$fajla$datapath)
    return(head(df))
  })
  
}
