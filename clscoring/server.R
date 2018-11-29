# server logic
server <- function(input, output, session){ #clientData, 
  
  #changes cell line selection based on tumour type selected
  observe({
    c_tum <- input$tum
    #data_type <- input$datatypeSelected
    #ucec_cl_x <- switch(data_type,)
    #brca_cl_x <- switch(data_type,)
    cls <- switch(c_tum,
                  UCEC = rownames(cl_cn_ucec), #c("AN3CA", "MFE"),
                  BRCA = rownames(brcacl_cn))
    updateSelectInput(session, "cl",
                      label = paste("Select a cell line of", input$tum,"type:"),
                      choices = cls)
  })
  
  #input$datatypeSelected
  
  
  vals <- reactiveValues(rf = NULL, rfmds = hist(c(1,2,3)), dtype=NULL,
                         pdat=NULL, cldat=NULL, histdat=NULL, histname="Default", i=0)
  
  #updating TCGA barcode classification
  observeEvent(input$barcode_update,
               {
                 req(input$fajla)
                 vals$histdat <<- read.csv(input$fajla$datapath)
                 vals$i<-vals$i+1
                 vals$histname <<- paste("Custom classification",vals$i)
                 #output$tcga_table <- head(vals$hisdat)
               })
  
  #building rf tumour model and generating the plot
  observeEvent(input$do_model,
               {
                 if(input$tum=="UCEC"){
                     vals$pdat <- switch(input$datatypeSelected,
                            cn=ucec_cn,
                            mut=ucec_mut,
                            exp=ucec_rsem)
                     vals$cldat <- switch(input$datatypeSelected,
                                          cn=cl_cn_ucec,
                                          mut=cl_maf_ucec,
                                          exp=cl_rsem_ucec)
                     vals$histdat<-UCEC_histology
                 }
                 
                 if(input$tum=="BRCA"){
                   vals$pdat <- switch(input$datatypeSelected,
                                       cn=brca_cn,
                                       mut=brca_mut,
                                       exp=brca_rsem)
                   vals$cldat <- switch(input$datatypeSelected,
                                       cn=brcacl_cn,
                                       mut=brcacl_mut,
                                       exp=brcalcl_rsem)
                   vals$histdat <- BRCA_histology
                 }
                 

                 
                environment(rf_default) <- environment() #isuse R je odvratan, a i ja
                vals$rf <<- rf_default(histology=vals$histdat,
                                        patient_data=vals$pdat)
                dev.off()
                MDSplot(vals$rf, training_set$class, palette=c(3,2))
                vals$rfmds <<- recordPlot()

                 
               }
  )
  
  #scoring cell line
  rez_cl_tum <- eventReactive(input$do_score,
                              {
                                if(is.null(vals$rf)){
                                  paste(input$tum, input$cl)
                                }
                                else{
                                  
                                  clpreds <- predict(vals$rf, vals$cldat, type="prob")
                                  paste("A subtype score of", input$cl, ":",
                                        clpreds[input$cl,1])
                                }
                                
                              }
  )
  
  rez_all_cl_plot <- eventReactive(input$do_score_2,
                              {
                                if(is.null(vals$rf)){
                                  paste(input$tum, input$cl)
                                }else{
                                  all_cl_barplot(cldf = vals$cldat, pat = vals$pdat , rf = vals$rf)
                                }

                              }
  )
  
  cl_dists <- eventReactive(input$do_dists,
                            {
                              
                              get_dists(cldata = vals$cldat, patdata = vals$pdat) %>%
                                cldist_min_per_type(vals$histdat) %>%
                                plot_shiny_distance()
                              
                            }
                            
                            )
  
  
  
  
  output$mds <- renderPlot(vals$rfmds)
  output$pred <- renderText(cancer_pred)
  output$rez <- renderText(rez_cl_tum())
  
  output$histname <- renderText(vals$histname)
  
  output$plotscore <- renderPlot({
      if(is.null(vals$rf)) return(NULL)
      rez_all_cl_plot() %>% ggplot(aes(x=reorder(Name, Class1),y=Class1, fill=Class2)) +
               geom_bar(stat = "identity", colour="black") +
               coord_flip() + scale_fill_viridis_c() +
               guides(fill=F) +
              labs(x="Name",y="Classifier score")
    })
  
  output$plotcldist <- renderPlot({
    cl_dists()
  })

  output$tcga_table <- renderTable({
    head(vals$histdat, 15)
  })
  
  output$tum <- renderText(input$tum)
  output$datatypeSelected <- renderText(
    switch(input$datatypeSelected,
           
           exp="Expression",
           cn="Copy number",
           mut="Mutation",
           combo="Combination")
    )
  
  
}
