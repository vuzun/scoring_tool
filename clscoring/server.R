# server logic
server <- function(input, output, session){ #clientData, 
  
  #changes cell line selection based on tumour type selected
  observe({
    c_tum <- input$tum
    cls <- switch(c_tum,
                  UCEC = rownames(cl_cn_ucec),
                  BRCA = rownames(brcacl_cn))
    updateSelectInput(session, "cl",
                      label = paste("Select a cell line of", input$tum,"type:"),
                      choices = cls)
  })
  
  #input$datatypeSelected
  
  
  #everything
  vals <- reactiveValues(rf = NULL, rfmds = NULL, dtype=NULL,
                         pdat=NULL, cldat=NULL, rez_a_cl=NULL, histdat=UCEC_histology,
                         histucec=UCEC_histology, histbrca=BRCA_histology,
                         histname="Default classification scheme", i=0,
                         building_flag=0,
                         ucecdefault=UCEC_histology, brcadefault=BRCA_histology,
                         error_samplesize=0)
  
  
  #updating TCGA barcode classification
  observeEvent(input$barcode_update,
               {
                 req(input$fajla)
                 userdef_hist <- read.csv(input$fajla$datapath)
                 datalist_for_ctype <- switch(
                   input$tum,
                   UCEC=list(ucec_rsem, ucec_cn, ucec_mut),
                   BRCA=list(brca_rsem, brca_cn, brca_mut)
                  )
                   
                 vals$histdat <<- do_new_hist_intersect(userdef_hist, datalist_for_ctype)
            
                 #vals$histdat <<- read.csv(input$fajla$datapath)
                 
                 vals$i<<-vals$i+1
                 vals$histname <<- paste("Custom classification scheme",vals$i)
                 output$tcga_table <- renderTable({
                   head(vals$histdat, 15)
                 })
               })
  
  
  #setting back the default classification scheme for the selected cancer type
  observeEvent(input$set_default_scheme,
               {
                 vals$histname <<- "Default classification scheme"
                 vals$histdat <<- ifelse(input$tum=="UCEC",vals$ucecdefault, vals$brcadefault)
               })
  
  #building rf tumour model and generating the plot
  observeEvent(input$do_model,
               {
                 progress <- Progress$new(min=0,max=2)
                 on.exit(progress$close())
                 progress$set(message="Processing...", value=1.5)
                 
                 vals$error_samplesize <<- 0
                 
                 if(input$tum=="UCEC"){
                     vals$pdat <- switch(input$datatypeSelected,
                            cn=ucec_cn,
                            mut=ucec_mut,
                            exp=ucec_rsem,
                            combo=data_combine(type="patient",ucec_cn, ucec_mut, ucec_rsem))
                     vals$cldat <- switch(input$datatypeSelected,
                                          cn=as.data.frame(cl_cn_ucec),
                                          mut=cl_maf_ucec,
                                          exp=cl_rsem_ucec,
                                          combo=data_combine(type="CL",
                                                             as.data.frame(cl_cn_ucec), cl_maf_ucec, cl_rsem_ucec,
                                                             patient_colnames = colnames(vals$pdat)))
                     vals$histdat<-UCEC_histology
                 }

                 if(input$tum=="BRCA"){
                   vals$pdat <- switch(input$datatypeSelected,
                                       cn=brca_cn,
                                       mut=brca_mut,
                                       exp=brca_rsem,
                                       combo=data_combine(type="patient",
                                                          brca_cn, brca_mut, brca_rsem))

                   vals$cldat <- switch(input$datatypeSelected,
                                       cn=brcacl_cn,
                                       mut=brcacl_mut,
                                       exp=brcalcl_rsem,
                                       combo=data_combine(type="CL",
                                                          brcacl_cn, brcacl_mut, brcalcl_rsem,
                                                          patient_colnames = colnames(vals$pdat)))
                   vals$histdat <- BRCA_histology
                 }
                 
                 vals$histdat <- vals$histdat[na.omit(match(rownames(vals$pdat), vals$histdat[[1]] )),]


                vals$building_flag <<- 1

                environment(rf_default) <- environment() #this is ugly
                rf_run <- tryCatch(rf_default(histology=vals$histdat,
                                        patient_data=vals$pdat),
                                     
                                     error=function(e){
                                       if(e$message=="Not enough samples for one of the classes."){
                                         vals$error_samplesize <<- 1
                                         vals$rfmds <<- NULL
                                       }else{
                                         stop(e$message)
                                       }
                                         }
                )
                
                vals$rf <<- rf_run$model
                vals$auc <- rf_run$auc
                
                if(vals$error_samplesize==0){

                  dev.new()
                  dev.control("enable")
                  plot(stats::cmdscale(1 - vals$rf$proximity, k = 2, eig = TRUE)$points,
                     col=ifelse(training_set$class==1, "orange", "blue"), pch=19,
                     xlab="Dim 1", ylab="Dim 2")
                  vals$rfmds <<- recordPlot()
                  dev.off()
                }


               }
  )
  
  #scoring cell line
  #rez_cl_tum <- eventReactive
  observeEvent(input$do_score,
                              {
                                if(is.null(vals$rf) | vals$error_samplesize>0){
                                  vals$rez_a_cl <<- NULL
                                  #paste(input$tum, input$cl)
                                }
                                else{
                                  
                                  clproximity <- predict(vals$rf,
                                                         rbind(vals$cldat[input$cl,],
                                                               training_set[,-ncol(training_set)]),
                                                         type="prob", proximity=TRUE)
                                  
                                  clpreds <- predict(vals$rf, vals$cldat, type="prob")
                                  dev.new()
                                  dev.control("enable")
                                  plot(stats::cmdscale(1 - clproximity$proximity, k = 3, eig = TRUE)$points,
                                       col=c("black",
                                             ifelse(training_set$class==1, "orange", "blue")),
                                       pch=c(4, rep(19, length(training_set$class))),
                                       xlab="Dim 1", ylab="Dim 2")
                                  vals$rfmds <<- recordPlot()
                                  dev.off()
                                  vals$rez_a_cl <<- paste("Subtype score of", input$cl, ":",
                                        clpreds[input$cl,1])
                                  

                                }
                                
                              }
  )
  
  #scoring all cell lines
  rez_all_cl_plot <- eventReactive(input$do_score_2,
                              {
                                if(is.null(vals$rf)){
                                  paste(input$tum, input$cl)
                                }else{
                                  all_cl_barplot(cldf = vals$cldat, pat = vals$pdat , rf = vals$rf)
                                }

                              }
  )
  
  #making the distances plot
  cl_dists <- eventReactive(input$do_dists,
                            {
                              
                              if(is.null(vals$pdat)) return(NULL)
                              
                              progress <- Progress$new(min=0,max=2)
                              on.exit(progress$close())
                              progress$set(message="Processing...", value=1.5)
                          
                              
                              get_dists(cldata = vals$cldat, patdata = vals$pdat) %>%
                                cldist_min_per_type(vals$histdat) %>%
                                plot_shiny_distance()
                              
                            }
                            
                            )
  
  
  
  

  output$error_sample <- reactive({
    vals$error_samplesize==1
  })
  outputOptions(output, 'error_sample', suspendWhenHidden=FALSE)

  output$auc <- renderText({
    if(is.null(vals$auc)) return(NULL)
    as.character(round(vals$auc,2))
  })
  
  output$mds <- renderPlot({
    if(is.null(vals$rfmds)) return(NULL)
    replayPlot(vals$rfmds)})
  
  output$rez <- renderText({
    if(is.null(vals$rf)) return(NULL)
    vals$rez_a_cl
    })
  #output$rez <- renderText(rez_cl_tum())
  
  output$histname <- renderText(vals$histname)
  output$histname2 <- renderText(vals$histname)
  
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
  
  output$classification_scheme <- downloadHandler(
    filename = "classification_scheme.csv", 
    content = function(file){
      write.csv(vals$histdat, file, row.names = F)###
    }
  )
  
  output$tum <- renderText(input$tum)
  output$datatypeSelected <- renderText(
    switch(input$datatypeSelected,
           
           exp="Expression",
           cn="Copy number",
           mut="Mutation",
           combo="Combination")
    )
  
  
}
