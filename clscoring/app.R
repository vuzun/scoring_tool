#cl_score_app/app.R
#needs new UI for the score of cl-tum


# TODO:
#   * add all datas
#   * specifying mut/cn/exp profile (default top x from rf)
#   * ?best cl for a tumour

 

#have a preproc file


library(shiny)

# Running
options(shiny.trace=TRUE)
shinyApp(ui = ui, server = server)

# cl in do dists is char???

# replacement has 524 rows, data has 513 - cldist_min_per_type ---> histology mismatch