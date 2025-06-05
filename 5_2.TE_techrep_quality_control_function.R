#function script
####find techrep outliers in 11/13 F_R_ratio value by CV #######
#FUNCTION copy from SRver5-1.R and modification in TE part
##function:find outlier in technically repeat##
##input: Blank600 is a data frame used to find outlier in technically repeat, 
##       default_rm_vector is a vector used to remove some strains you did not want
##       Categorical_variable is the vector for group , melt and dcast function
##       CV_cutoff is the cutoff of CV value
##       OD_time is the time to be choosed for techrep OD filter(find outlier) 
##output: ODoutlier_observation contain the strain outlier technically repeat location
##reconstruct######
##reconstruct in SIver5-1 part####
find_TE_techrep_outlier <- function(TE_techrep_data,default_rm_vector = NA,CV_cutoff=0.5,Categorical_variable= NULL){
  
  F_R_ratio_data <-  TE_techrep_data %>% filter(!strain %in% default_rm_vector) %>% 
    ## remove default strain vector  e.g. BY and blank
    mutate(F_R_ratio = `F`/`R`)
  F_R_ratio_Techrep_CV <- F_R_ratio_data %>% 
    group_by_at(c("strain","mutant_type","biorep","panel",
                  Categorical_variable[!Categorical_variable %in% c("rep_location_name","rep_location_num","tube_type")])) %>% 
    mutate(Group_Size = n()) %>%
    mutate(mean_F_R_ratio = mean(F_R_ratio),sd_F_R_ratio = sd(F_R_ratio)) %>%
    mutate(CV = sd_F_R_ratio/mean_F_R_ratio) %>%
    filter(CV >= CV_cutoff)
  
  if(nrow(F_R_ratio_Techrep_CV) == 0 ){
    return(F_R_ratio_Techrep_CV %>% select(-Group_Size))
  }else{
    F_R_outlier_observation <- F_R_ratio_Techrep_CV%>%
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
      }) %>% select(-Group_Size,-sum_of_Manhattan_Distance_vector)
    
    return(F_R_outlier_observation)
  }
}

#####tehcrep clean by anti_join ######
#Blank600 clean
TE_techrep_clean <- function(data,TEoutlier_data,remove_col=c(9:12)){
  Techrep_clean_data <- anti_join(data,TEoutlier_data,by = names(TEoutlier_data)[-remove_col])
  return(Techrep_clean_data)
}

