#description : calculate 11 type strain F_R_ratio , summary to TE_11_biorep_data , 
#and filter TE_11_biorep_data to TE_11_biorep_clean_bio_data , summary to TE_11_data
#save TE_11_biorep_data(biorep type) , TE_11_biorep_outlier_data(biorep type) , 
#TE_11_biorep_clean_bio_data(biorep type) , TE_11_data(strain type, for TE calculate)

#this function also suit for TE 13 type data 

#function script
######TE 11 biorep filter by CV ########
find_TE_11_biorep_outlier <- function(TE_11_Data_contain_bio,cutoff_type="CV",CV_cutoff=0.5,Categorical_variable=NULL){
   if(cutoff_type == "CV" ){
    
    TE_11_biorep_out_of_CV <- TE_11_Data_contain_bio %>%
      group_by_at(c("strain","mutant_type",
                    Categorical_variable[!Categorical_variable %in% c("rep_location_name","rep_location_num","tube_type")])) %>%
      mutate(Group_Size = n()) %>%
      mutate(sd_F_R_ratio = sd(F_R_ratio),F_R_ratio_mean_bio = mean(F_R_ratio)) %>%
      mutate(CV = sd_F_R_ratio/F_R_ratio_mean_bio)%>%
      filter(CV >= CV_cutoff)
    
    if(nrow(TE_11_biorep_out_of_CV) == 0 ){
      return(TE_11_biorep_out_of_CV %>% select(-Group_Size))
    }else{
      TE_11_biorep_outlier <- TE_11_biorep_out_of_CV%>%
        do({
          #create index for debug
          the_strain = unique(.$strain)
          cat("Working on...", the_strain,"\n")
          
          myDF <- .
          sum_of_Manhattan_Distance_vector <- c()
          for (i in 1:length(myDF[["F_R_ratio"]])) {
            value <- myDF[["F_R_ratio"]][i]
            other_value <- myDF[["F_R_ratio"]][-i]
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
  return(TE_11_biorep_outlier)
}

