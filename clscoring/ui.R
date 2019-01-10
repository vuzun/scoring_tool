# Define UI for application that draws a histogram
ui <- navbarPage("Cell line scoring",
                 
                 tabPanel("Start",
                          sidebarLayout(
                            sidebarPanel("Select tumour and data type for the machine learning model:",br(),
                                         selectInput("tum", label = "Tumour type:", 
                                                     list("Tumour cohorts" = c("UCEC", "BRCA"))),
                                         br(),
                                         radioButtons("datatypeSelected", "Data types",
                                                      choices = list("Expression" = "exp",
                                                                     "Copy number" = "cn",
                                                                     "Mutation" = "mut",
                                                                     "Combination" = "combo"),
                                                      selected = "cn"
                                         ),
                                         actionButton("do_model", "Build a model"),
                                         hr(),
                                         selectInput("cl", label = "Select a cell line to score:",
                                                     choices = c("AN3CA", "MFE-280", "JHUEM-3") ),
                                         actionButton("do_score", "Score the cell line"),
                                         hr(),
                                         br(),
                                         tags$b("Calculate and plot all the cell line scores of this cancer type:"),br(),
                                         actionButton("do_score_2",
                                                      "Score and plot all")#,
                                         #hr(), hr(),
                                         
                            ),
                            mainPanel(textOutput("histname2"),
                                      br(),
                                      "Tumour classifier performance - AUC= 0.87",
                                      hr(),
                                      plotOutput("mds"),
                                      textOutput("rez"),
                                      textOutput("cl_dist"),
                                      plotOutput("plotscore"))
                          )
                 ),
                 tabPanel("Classification schemes",
                          sidebarLayout(
                            sidebarPanel(
                              fileInput("fajla",
                                        "Upload own TCGA subtypes for the selected cancer type",
                                        accept = c("text/csv", 
                                                   "text/comma-separated-values,text/plain", ".csv")),
                              actionButton("barcode_update","Update the classification scheme"),
                              hr(),
                              downloadButton("classification_scheme","Download the current scheme")
                              
                            ),
                            mainPanel("You can provide a different classification scheme or download the current one.",br(),
                                      "The classsification scheme table needs to have 2 columns:",br(),
                                      tags$ul(
                                        tags$li("Valid TCGA barcodes"),
                                        tags$li("Column containing only 2 different values")
                                      ),
                                      
                                      "Once the custom table is uploaded and classification updated, rebuild the model.", br(),
                                      hr(),br(),
                                      tags$b("Top of the current classification table:"),
                                      br(),
                                      tableOutput("tcga_table"))
                          )
                 ),
                 

                 
                 tabPanel("Cell line distances",
                          sidebarLayout(sidebarPanel(
                            h3("Current settings"),br(),
                            tags$b("Tumour type:"), textOutput("tum"),
                            br(),
                            tags$b("Classification scheme:"),textOutput("histname"),
                            hr(),
                            tags$b("Data type:"), textOutput("datatypeSelected")
                            
                            
                          ),mainPanel("Aside from the classification model, cell line similarity to patient samples can also be helpful.",
                                      "Here, number of patient samples that have a certain cell line as their closest cell line can be calculated for selected cancer type, data and classification scheme.",
                          br(),
                          actionButton("do_dists", "Count closest cell lines"),
                          br(),br(),br(),
                          plotOutput("plotcldist")
                          ))
                 ),
                 tabPanel("About",
                          sidebarLayout(sidebarPanel("This application was developed as a part of the PhD project: 'Application of computational
                                                     methods to the assessment of clinical relevance of preclinical cancer models'",
                                                     br(),
                                                     "Author: Vladimir Uzun, The University of Sheffield."

                          ),mainPanel(
                          ))
                 )
                 
)