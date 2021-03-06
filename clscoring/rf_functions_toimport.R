
#' random forest wrapper for training/testiing samples with 'randomForest'
#'
#' @param histology Nx2 dataframe with class in the 2nd column.
#' @param patient_data NxM dataframe with already aligned rows to histology.
#' @param n Random forest parameter: number of trees. Defaults to 800.
#' @param mt Random forest parameter: number of different variable tries per node.
#' Defaults to sqrt(number of variables).
#' @param viz Flag for plotting or no plotting of MDS.
#'
#' @return R object of class 'randomForest'
#' @export
#'
#' @examples rf_default(UCEC_histology, UCEC_copy_number_TCGA)
rf_default <- function(histology=CN_histology,
                       patient_data=CN_endometrial_TCGA,
                       n=500, newmt=(-1), viz=F){

  
  ##plot(stats::cmdscale(1 - rfUC$proximity, eig = TRUE, k = 2)$points, col=ifelse(training_set$class==1, "yellow", "purple"), pch=19)
  
  histology <- histology[na.omit(match(rownames(patient_data), histology[[1]] )),]
  patient_data <- patient_data[na.omit(match(histology[[1]],rownames(patient_data))),]
  
  sample_hist <- histology[[2]]
  
  classes <- unique(sample_hist)

  
  mt <- sqrt(ncol(patient_data))
  if(newmt>0) mt <- newmt

  primary_class <- names(which.max(table(sample_hist)))

  histology_vector<<-ifelse(sample_hist == primary_class,0,1)
  histology_vector[is.na(histology_vector)]<<-0
  
  size_trainclass <<- min(table(histology_vector))*0.8
  
  # throws error, but this specific one should be caught higher up
  if(size_trainclass < 10){
    stop("Not enough samples for one of the classes.")
  }
  
  ser_inds <<- which(histology_vector==1)
  oid_inds <<- which(histology_vector==0)
  ser_train <<- sample(ser_inds, size_trainclass)
  oid_train <<- sample(oid_inds, size_trainclass)
  
  training_inds<<- c(ser_train, oid_train)

  patient_data$class <- as.factor(histology_vector)

  training_set<<-as.data.frame(patient_data[training_inds,])
  testing_set<<-as.data.frame(patient_data[-training_inds,])
  
  t_premodel <- Sys.time()
  
  RFmodel <- randomForest(class ~ ., data=training_set,
                           importance=TRUE, ntree=n, keep.forest=TRUE, proximity=TRUE, mtry=mt)
  paste0("--", Sys.time()-t_premodel) %>% print
  #how to interpret this better?
  if(viz==T)
    MDSplot(RFmodel, training_set$class, palette=c("blue","orange"))
  
  
  auc_score <- prediction(predict(RFmodel, newdata=testing_set,type="prob")[,2], testing_set$class) %>% 
    performance("auc")
  return(list("model"=RFmodel, "auc"=as.numeric(auc_score@y.values)))
  
  
}


data_combine <- function(type, cn, mut, exp, histology=NULL, subtypesize_min=100, patient_colnames=NULL){
  
  colnames(cn) <- colnames(cn) %>% paste0("_CN")
  colnames(mut) <- colnames(mut) %>% paste0("_MUT")
  colnames(exp) <- colnames(exp) %>% paste0("_EXP")
  
  if(type=="patient"){

    int_cnmut<- intersect(rownames(cn),rownames(mut)) 
    int_cnexp <- intersect(rownames(cn),rownames(exp))
    int_expmut <- intersect(rownames(exp),rownames(mut))
    int_all <- intersect(intersect(rownames(cn),rownames(mut)), rownames(exp))
    
    if(length(int_all) > subtypesize_min){
      cn <- cn[which(rownames(cn) %in% int_all),]
      mut <- mut[which(rownames(mut) %in% int_all),]
      exp <- exp[which(rownames(exp) %in% int_all),]
      
      #rescaling
      
      final_merge <- do.call(cbind, list(cn, mut, exp))
      
    }else if(max(int_expmut %>% length, int_cnexp %>% length, int_cnmut %>% length)>subtypesize_min){
      int_lens <- c(int_expmut %>% length, int_cnexp %>% length(), int_cnmut %>% length())
      datas <- switch(which.max(int_lens),
                      list(exp,mut),
                      list(cn, exp),
                      list(cn, mut)
             )
      final_merge <- do.call(cbind, datas)
      
      
    }else{
      print("Too small.") #throw an error?
      final_merge <- cn
    }
  }
  else{

    int_all <- intersect(intersect(rownames(cn), rownames(exp)), rownames(mut))
    
    cn <- cn[match(int_all, rownames(cn)),which(colnames(cn) %in% patient_colnames)]
    exp <- exp[match(int_all, rownames(exp)),which(colnames(exp) %in% patient_colnames)]
    mut <- mut[match(int_all, rownames(mut)),which(colnames(mut) %in% patient_colnames)]
    
  
        
  
    final_merge <- do.call(cbind, list(cn, exp, mut))
    final_merge <- final_merge[,match(patient_colnames ,colnames(final_merge))]
    final_merge %>% dim %>% paste("- dimension of the final merge")
    final_merge
  
  }
}

all_cl_barplot <- function(cldf, pat, rf){
  cldf <- rbind(pat[1,],cldf)
  cldf <- cldf[-1,]
  
  rez_list <-  predict(rf, cldf, type="prob") %>%
    as.data.frame %>% rownames_to_column(var = "Name")
  colnames(rez_list) <- c("Name","Class1","Class2")
  
  rez_list
}

do_distances_to_patients <- function(cl, pat_data){
  #needs rescaling when data combined
  #pat_data <- sapply(pat_data, function(co) as.character(co) %>% as.numeric) %>% as.data.frame
  #cl <- as.numeric(cl)
  #browser()
  apply(pat_data, 1, function(p_row) mean(abs(p_row-cl)) )
}


get_dists <- function(cldata, patdata){
 pat_rownames <- rownames(patdata)
 patdata <- data.frame(patdata %>% apply(2, function(x) as.numeric(as.character(x))))
 rownames(patdata) <- pat_rownames
 #browser()
 as.data.frame(apply(cldata, 1,
                     function(cl) do_distances_to_patients(cl %>% as.character %>% as.numeric, patdata)))
  
}

cldist_min_per_type <- function(cldist, histdata){

  minmax <- cldist %>%
    apply(1, function(ro) colnames(cldist)[c(which.min(ro), which.max(ro))]) %>%
    t %>% as.data.frame
  colnames(minmax) <- c("Min","Max")
  minmax$class <- histdata[[2]]
  
  classes <- as.character(unique(minmax$class))
  
  #browser()
  
  minmax <- minmax %>% group_by(Min) %>%
    summarise(c1=sum(as.double(class==classes[1])), c2=sum(as.double(class==classes[2])))
  
  colnames(minmax)[2:3] <- classes
  
  adding_missing_cls <- function(minmax, cldist){
    missing <- setdiff(cldist %>% colnames, minmax[[1]])
    
    addon_minmax <- cbind(missing, rep(0, length(missing)),rep(0, length(missing)))
    colnames(addon_minmax) <- colnames(minmax)
    
    minmax <- rbind(minmax, addon_minmax)
    minmax
    
  }
  
  minmax %>% adding_missing_cls(cldist)
}

plot_shiny_distance <- function(mind){
  mind %<>% gather(subtype, counts, -Min)
  mind$counts  %<>% as.numeric
  mind$subtype  %<>% as.factor
  colnames(mind) <- c("Name", "Subtype", "counts")
  
  mind %>% ggplot(aes(x=Name, y=counts, fill=Subtype)) +
    geom_bar(stat="identity", position = position_dodge(), colour="black") +
    theme(axis.text.x = element_text(angle = 45), axis.text.y = element_text(hjust = 0.5, vjust = 0.5)) +
    scale_fill_viridis_d() +
    labs(y="Samples with this cell line as closest")+
    coord_flip()
}


do_new_hist_intersect <- function(user_hist, datalist){
  bcodes_user <- user_hist[[1]]
  bcodes_exp <- rownames(datalist[[1]])
  bcodes_cn <- rownames(datalist[[2]])
  bcodes_mut <- rownames(datalist[[3]])
  
  paste(length(bcodes_user), length(bcodes_cn), length(bcodes_exp), length(bcodes_mut)) %>% 
    print
  
  bcodes_matched <- intersect(bcodes_user, bcodes_cn) %>% 
    union(intersect(bcodes_user, bcodes_exp)) %>% 
    union(intersect(bcodes_user, bcodes_mut))
  
  print(paste("Matched bcodes:", length(bcodes_matched)))
  
  user_hist[which(bcodes_matched %in% bcodes_user),]
}


#' Title
#'
#' @return
#' @export
#'
#' @examples
#cl_cont_on_rf <- function(){
#  RF_cellline_prediction_continuous <- predict(RFmodel, newdata=CN_endometrial_CCLE, type="prob")
#  kable(RF_cellline_prediction_continuous[order(RF_cellline_prediction_continuous[,1]),], format = "markdown")
  
#brcasamples_kal <- as.data.frame(t(brcasamples_kal)) #wtf?!
  
#  rf_kal<-randomForest(class ~ ., data=wee[c(inds_train_nontn, inds_train_tn),], importance=TRUE,
#                       ntree=100, keep.forest=TRUE, proximity=TRUE)
  
  
  
#  wee_cl <- log(as.data.frame(t(brcacl_kal))+0.0001)
  
#  cl_predict_brcakal_rf <-predict(rf_kal, newdata = wee_cl, type="prob")
  
  
#}