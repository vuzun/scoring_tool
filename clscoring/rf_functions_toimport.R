library(randomForest)
#library(reprtree)
library(ROCR)
library(knitr)
#library(maftools)

#separated from .Rmd for easier importing
#and because I'm a crummy coder

silly_two <- function(n=5){
  Sys.sleep(4)
  return(dim(TCGA_histology)[[1]]*n)
  
}

rf_default <- function(n=800,mt=70, viz=F){

  sample_hist <- TCGA_histology[[2]]
  
  classes <- unique(sample_hist)
  

  primary_class <- names(which.max(table(sample_hist)))

  
  histology_vector<<-ifelse(sample_hist == primary_class,0,1)
  histology_vector[is.na(histology_vector)]<<-0
  
  size_trainclass <<- min(table(histology_vector))*0.8
  
  ser_inds <<- which(histology_vector==1)
  oid_inds <<- which(histology_vector==0)
  ser_train <<- sample(ser_inds, size_trainclass)
  oid_train <<- sample(oid_inds, size_trainclass)
  
  training_inds<<- c(ser_train, oid_train)
  
  training_set<<-as.data.frame(CN_endometrial_TCGA[training_inds,])
  testing_set<<-as.data.frame(CN_endometrial_TCGA[-training_inds,])
  
  training_set$class<<-as.factor(histology_vector[training_inds])
  testing_set$class<<-as.factor(histology_vector[-training_inds])
  
  
  
  RFmodel <<- randomForest(class ~ ., data=training_set, importance=TRUE, ntree=n, keep.forest=TRUE, proximity=TRUE, mtry=mt)
  
  #how to interpret this better?
  if(viz==T)
    MDSplot(RFmodel, training_set$class, palette=c(3,2))
  
  cancer_pred <<- predict(RFmodel, newdata=testing_set)
  return(RFmodel)
  
  
}

#mess
#?wee
#?wee_cl
cl_cont_on_rf <- function(){
  RF_cellline_prediction_continuous <- predict(RFmodel, newdata=CN_endometrial_CCLE, type="prob")
  kable(RF_cellline_prediction_continuous[order(RF_cellline_prediction_continuous[,1]),], format = "markdown")
  
#brcasamples_kal <- as.data.frame(t(brcasamples_kal)) #wtf?!
  
  rf_kal<-randomForest(class ~ ., data=wee[c(inds_train_nontn, inds_train_tn),], importance=TRUE,
                       ntree=100, keep.forest=TRUE, proximity=TRUE)
  
  
  
  wee_cl <- log(as.data.frame(t(brcacl_kal))+0.0001)
  
  cl_predict_brcakal_rf <-predict(rf_kal, newdata = wee_cl, type="prob")
  
  
}