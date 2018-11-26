# Define UI for application that draws a histogram
ui <- navbarPage("Cell line scoring",
                 
                 tabPanel("Start",
                          sidebarLayout(
                            sidebarPanel("Select tumour type and cell line:",
                                         selectInput("tum", label = "Tumour type:", 
                                                     list("Tumour cohorts" = c("UCEC", "BRCA"))),
                                         actionButton("do_model", "Build a model"),
                                         hr(),
                                         selectInput("cl", label = "Pick a CL of x type",
                                                     choices = c("AN3CA", "MFE-280", "JHUEM-3") ),
                                         actionButton("do_score", "Score the cell line")#,
                                         #hr(), hr(),
                                         
                            ),
                            mainPanel("Tumour classifier performance - F1= ",
                                      plotOutput("mds"),
                                      hr(),
                                      textOutput("rez"),
                                      textOutput("cl_dist"),
                                      tableOutput("feele"))
                          )
                 ),
                 
                 tabPanel("Update classifications",
                          sidebarLayout(
                            sidebarPanel(
                              fileInput("fajla",
                                        "Upload own TCGA subtypes for the selected type",
                                        accept = c("text/csv", 
                                                   "text/comma-separated-values,text/plain", ".csv")),
                              actionButton("barcode_update","Update the barcode classification")
                            ),
                            mainPanel("Hi!",br(),"Here you can provide own set of classifications in the format of a table of \"TCGA barcodes\"~\"class\" pairs.",br(),
                                      "TCGA barcodes need to be in proper format. Click update, then build mdoel again?")
                          )
                 ),
                 
                 tabPanel("Data type integration",
                          checkboxGroupInput("typeSelected",
                                             h3("Data types used"),
                                             choices = list("Expression" = 1,
                                                            "Copy number" = 2,
                                                            "Mutation" = 3),
                                             selected = c(1,2,3)
                                             )
                   
                 )
                 
)