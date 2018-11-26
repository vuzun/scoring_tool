#cl_score_app/app.R
#needs new UI for the score of cl-tum


# TODO:
#   * add mut/exp
#   * specifying mut/cn/exp profile (default top x from rf)
#   * ?best cl for a tumour
#   *





#have a preproc file


library(shiny)

# Running
shinyApp(ui = ui, server = server)
