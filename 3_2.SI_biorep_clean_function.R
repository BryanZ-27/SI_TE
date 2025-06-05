####Categorical_variable ###
Categorical_variable <- c("strains_group","rep_location_name","rep_location_num",
                          "tube_type","tube_type_YPD","tube_type_SC",
                          "glu_concentration")


###function:Clean data after SI calculation  #This function is copy from SI_Data_for_QTL_and_other_analysis.R
###   purpose: to remove high Sd strains or not good biorep 
###input:SIData_Raw is the SI data that contain strains SI data , aggregate the biorep and have the sd_SI (biorep)
###      SIData_contain_bio is the SI data that contain each strain biorep SI data without summarize calculation
###      cutoff_type is the cutoff calculation method , now just have the boxplot
###      Categorical_variable is the vector for group , melt and dcast function
find_SI_biorep_outlier <- function(SIData_contain_bio,cutoff_type="CV",CV_cutoff=0.5,Categorical_variable=NULL){
  if (cutoff_type == "boxplot" ) {
    SIData_Raw <- SIData_contain_bio%>%
      group_by_at(c("strains",Categorical_variable[!Categorical_variable %in% 
                                                     c("rep_location_name","rep_location_num","tube_type")])) %>%
      dplyr::summarize(sd_SI = sd(SI),SI_mean_bio = mean(SI))
    #use box plot to find sd outlier strains
    outlier_object <- boxplot(SIData_Raw$sd_SI)
    outlier_cutoff <- outlier_object$out %>% sort() %>%.[1]
    
    SI_biorep_out_of_SD <- SIData_contain_bio%>%
      group_by_at(c("strains",Categorical_variable[!Categorical_variable %in% 
                                                     c("rep_location_name","rep_location_num","tube_type")])) %>%
      mutate(Group_Size = n()) %>%
      mutate(sd_SI = sd(SI),SI_mean_bio = mean(SI)) %>%
      filter(sd_SI >= outlier_cutoff)
    
    if(nrow(SI_biorep_out_of_SD) == 0 ){
      return(SI_biorep_out_of_SD %>% select(-Group_Size))
    }else{
      SI_biorep_outlier <- SI_biorep_out_of_SD%>%
        do({
          myDF <- .
          sum_of_Manhattan_Distance_vector <- c()
          for (i in 1:length(myDF[["SI"]])) {
            value <- myDF[["SI"]][i]
            other_value <- myDF[["SI"]][-i]
            sum_of_Manhattan_Distance <- sum(abs(value - other_value))
            sum_of_Manhattan_Distance_vector <- c(sum_of_Manhattan_Distance_vector,sum_of_Manhattan_Distance)
          }
          outputDF <- cbind(as.data.frame(myDF),sum_of_Manhattan_Distance_vector)
          if(outputDF$Group_Size[1]>3){
            outputDF <- outputDF %>%
              arrange(sum_of_Manhattan_Distance_vector) %>%
              filter(row_number() >3)
          }else if(outputDF$Group_Size[1]==3){
            outputDF <- outputDF %>%
              arrange(sum_of_Manhattan_Distance_vector) %>%
              filter(row_number() >2)
          }else{
            outputDF <- outputDF
          }
          outputDF
        }) %>% select(-Group_Size,-sum_of_Manhattan_Distance_vector)
      }
    
  }else if(cutoff_type == "CV" ){
    
    SI_biorep_out_of_CV <- SIData_contain_bio %>%
      group_by_at(c("strains",Categorical_variable[!Categorical_variable %in% 
                                                     c("rep_location_name","rep_location_num","tube_type")])) %>%
      mutate(Group_Size = n()) %>%
      mutate(sd_SI = sd(SI),SI_mean_bio = mean(SI)) %>%
      mutate(CV = sd_SI/SI_mean_bio)%>%
      filter(CV >= CV_cutoff)
    
    if(nrow(SI_biorep_out_of_CV) == 0 ){
      return(SI_biorep_out_of_CV %>% select(-Group_Size))
    }else{
      SI_biorep_outlier <- SI_biorep_out_of_CV%>%
        do({
          myDF <- .
          sum_of_Manhattan_Distance_vector <- c()
          for (i in 1:length(myDF[["SI"]])) {
            value <- myDF[["SI"]][i]
            other_value <- myDF[["SI"]][-i]
            sum_of_Manhattan_Distance <- sum(abs(value - other_value))
            sum_of_Manhattan_Distance_vector <- c(sum_of_Manhattan_Distance_vector,sum_of_Manhattan_Distance)
          }
          outputDF <- cbind(as.data.frame(myDF),sum_of_Manhattan_Distance_vector)
          if(outputDF$Group_Size[1]>3){
            outputDF <- outputDF %>%
              arrange(sum_of_Manhattan_Distance_vector) %>%
              filter(row_number() >3)
          }else if(outputDF$Group_Size[1]==3){
            outputDF <- outputDF %>%
              arrange(sum_of_Manhattan_Distance_vector) %>%
              filter(row_number() >2)
          }else{
            outputDF <- outputDF
          }
          outputDF
        }) %>% select(-Group_Size,-sum_of_Manhattan_Distance_vector)}
  }else{
    stop("not normal cutoff type")
  }
  return(SI_biorep_outlier)
}

###function:Clean data after SI calculation  #This function is copy from find_SI_biorep_outlier function and change for SR
###   purpose: to remove high Sd strains at day X "Survival Rate" or not good biorep 
###   input:
###      SRData_contain_bio is the Survival Rate data that contain each strain biorep Survival Rate data 
###      at DIFFERENT DAY TIME without summarize calculation
###      cutoff_type is the cutoff calculation method , now just have the CV method.
###      Categorical_variable is the vector for group , melt and dcast function

find_SR_biorep_outlier <- function(SRData_contain_bio,cutoff_type="CV",CV_cutoff=0.5,Categorical_variable=NULL){
  if (cutoff_type == "boxplot" ) {
    #### This argument function is not finish in SR function. NOW IS JUST CAN USE THE cutoff_type=="CV" ARGUMENT function.
  }else if(cutoff_type == "CV" ){
    SR_biorep_out_of_CV <- SRData_contain_bio %>%
      group_by_at(c("strains","time",Categorical_variable[!Categorical_variable %in% 
                                                            c("rep_location_name","rep_location_num","tube_type")])) %>%
      mutate(Group_Size = n()) %>%
      mutate(sd_SR = sd(Survival_n),SR_mean_bio = mean(Survival_n)) %>%
      mutate(CV = sd_SR/SR_mean_bio)%>%
      filter(CV >= CV_cutoff)
    
    if(nrow(SR_biorep_out_of_CV) == 0 ){
      return(SR_biorep_out_of_CV %>% select(-Group_Size))
    }else{
      SR_biorep_outlier <- SR_biorep_out_of_CV%>%
        do({
          myDF <- .
          sum_of_Manhattan_Distance_vector <- c()
          for (i in 1:length(myDF[["Survival_n"]])) {
            value <- myDF[["Survival_n"]][i]
            other_value <- myDF[["Survival_n"]][-i]
            sum_of_Manhattan_Distance <- sum(abs(value - other_value))
            sum_of_Manhattan_Distance_vector <- c(sum_of_Manhattan_Distance_vector,sum_of_Manhattan_Distance)
          }
          outputDF <- cbind(as.data.frame(myDF),sum_of_Manhattan_Distance_vector)
          if(outputDF$Group_Size[1]>3){
            outputDF <- outputDF %>%
              arrange(sum_of_Manhattan_Distance_vector) %>%
              filter(row_number() >3)
          }else if(outputDF$Group_Size[1]==3){
            outputDF <- outputDF %>%
              arrange(sum_of_Manhattan_Distance_vector) %>%
              filter(row_number() >2)
          }else{
            outputDF <- outputDF
          }
          outputDF
        }) %>% select(-Group_Size,-sum_of_Manhattan_Distance_vector)}
  }else{
    stop("not normal cutoff type")
  }
  return(SR_biorep_outlier)
}

