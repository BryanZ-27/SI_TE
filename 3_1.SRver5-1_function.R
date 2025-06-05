#library
library(reshape2)
library(stringr)
library(dplyr)

####Categorical_variable ###
Categorical_variable <- c("strains_group","rep_location_name","rep_location_num",
                          "tube_type","tube_type_YPD","tube_type_SC",
                          "glu_concentration")

##function:find outlier in technically repeat##
##input: Blank600 is a data frame used to find outlier in technically repeat, 
##       default_rm_vector is a vector used to remove some strains you did not want
##       Categorical_variable is the vector for group , melt and dcast function
##       CV_cutoff is the cutoff of CV value
##       OD_time is the time to be choosed for techrep OD filter(find outlier) 
##output: ODoutlier_observation contain the strain outlier technically repeat location
find_techrep_outlier <- function(Blank600_data,default_rm_vector = NA,CV_cutoff=0.5,OD_time="t50",Categorical_variable= NULL){
  
  OD_at_t_time <-  Blank600_data %>% filter(!strains %in% default_rm_vector) %>% 
    ## remove default strain vector  e.g. BY and blank
    melt(id.var = c("strains","biorep","techrep","time","panel",Categorical_variable)) %>%
    filter(variable == OD_time)
  OD_at_t_time_Techrep_CV <- OD_at_t_time %>% 
    group_by_at(c("strains","biorep","time","panel",
                  Categorical_variable[!Categorical_variable %in% c("rep_location_name","rep_location_num","tube_type")])) %>% 
    mutate(Group_Size = n()) %>%
    mutate(mean_value = mean(value),sd_value = sd(value)) %>%
    mutate(CV = sd_value/mean_value) %>%
    filter(CV >= CV_cutoff)
  
  if(nrow(OD_at_t_time_Techrep_CV) == 0 ){
    return(OD_at_t_time_Techrep_CV %>% select(-Group_Size))
  }else{
    ODoutlier_observation <- OD_at_t_time_Techrep_CV%>%
      do({
        myDF <- .
        sum_of_Manhattan_Distance_vector <- c()
        for (i in 1:length(myDF[["value"]])) {
          value <- myDF[["value"]][i]
          other_value <- myDF[["value"]][-i]
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
    
    return(ODoutlier_observation)
  }
}

#Blank600 clean
Blank_clean <- function(Blank600_data,ODoutlier_data){
  Blank_clean_data <- anti_join(Blank600_data,ODoutlier_data,by = names(ODoutlier_data)[-c(13,14,15,16,17)])
  return(Blank_clean_data)
}

###dbt
####clean Lnblank data by ODoutlier_data
LnBlank_clean <- function(LnBlank_data,ODoutlier_data){
  LnBlank_clean_data <- anti_join(LnBlank_data,ODoutlier_data,by = names(ODoutlier_data)[-c(13,14,15,16,17)])
  return(LnBlank_clean_data)
}

####doubling time calculation###
####function:doubling_time calculation , the outype base on the "well","retain_biorep","is.mean","sd" argument#######
####output type: well (2types,method:mean or row_number)
###############  retain_biorep(4types,method:mean or row_number,sd:have or not have techrep)
###############  strain(4types,method:mean or row_number,sd:have or not have biorep)
###output:
doubling_time_calculation <- function(LnBlankData,default_rm_vector = NA, Categorical_variable = NULL,
                                      OD_interval = c(0.2,0.5),well = FALSE, retain_biorep = FALSE, 
                                      is.mean = FALSE, daytime = c("day2") , sd = FALSE , is.all = F){
  dfDbt_ALL <- LnBlankData %>% filter(!strains %in% default_rm_vector) %>%
    melt(id.var = c("strains","biorep","techrep","time","panel",Categorical_variable)) %>%
    filter(strains != "blank") %>%
    group_by_at(c("strains","biorep","techrep","time","panel",Categorical_variable)) %>%
    mutate(doubling_time = 10*log(2)/ c(NA,diff(value)) ) %>%
    filter(value >= log(OD_interval[1]) & value <= log(OD_interval[2])) %>%
    filter(!variable %in% c(paste("t",seq(0,30,10),sep = ""))) %>% 
    filter(!doubling_time <= 0) %>%
    arrange(doubling_time)
  if(is.all == T){
    return(dfDbt_ALL)
  }
  #####calculate type: mean or second median##
  if (is.mean == T) {
    mean_dfDbt <- dfDbt_ALL %>% 
      dplyr::summarize(doubling_time = mean(doubling_time)) %>%
      ungroup()
  }else{
    dfDbt <- dfDbt_ALL %>% 
      filter(row_number() == 2) %>% 
      ungroup() %>%
      select(c(-variable,-value)) 
  }
  ####output type: well,biorep,strain
  if (well == T) {
    if (is.mean == T) {   
      return(mean_dfDbt)
    }else{
      return(dfDbt)
    }
  }else if(retain_biorep == T){
    if(is.mean == T){
      if (sd == FALSE) {
        mean_dfDbt_contain_biorep <-  mean_dfDbt %>% 
          filter(time %in% daytime) %>%   #### in mean method ,use day2 dbt as strain doubling_time 
          group_by_at(c("strains","biorep","time","panel",
                        Categorical_variable[!Categorical_variable %in% 
                                               c("rep_location_name","rep_location_num","tube_type")])) %>%
          dplyr::summarize(doubling_time = mean(doubling_time))
        return(mean_dfDbt_contain_biorep)
      }else{
        mean_dfDbt_contain_biorep <- mean_dfDbt %>% 
          filter(time %in% daytime) %>% 
          group_by_at(c("strains","biorep","time","panel",
                        Categorical_variable[!Categorical_variable %in% 
                                               c("rep_location_name","rep_location_num","tube_type")])) %>%
          dplyr::summarize(sd_tech = sd(doubling_time),doubling_time = mean(doubling_time))   ####this sd is techrep sd##
        return(mean_dfDbt_contain_biorep)
      }
      
    }else{
      if (sd == FALSE) {
        dfDbt_contain_biorep <- dfDbt_ALL %>% 
          filter(row_number() == 2) %>% 
          ungroup()%>% 
          select(c(-variable,-value)) %>% 
          group_by_at(c("strains","biorep","time","panel",
                        Categorical_variable[!Categorical_variable %in% 
                                               c("rep_location_name","rep_location_num","tube_type")])) %>%
          dplyr::summarize(doubling_time = mean(doubling_time))  
        return(dfDbt_contain_biorep)
      }else{
        dfDbt_contain_biorep <- dfDbt_ALL %>% 
          filter(row_number() == 2) %>% 
          ungroup()%>% 
          select(c(-variable,-value)) %>% 
          group_by(strains,biorep,time,strains_group,panel) %>% 
          dplyr::summarize(sd_tech = sd(doubling_time),doubling_time = mean(doubling_time))  
        ####contain other days data ,this sd is techrep sd##
        return(dfDbt_contain_biorep)
      }
      
    }
  }else{
    if(is.mean == T){
      if (sd == FALSE){
        mean_Dbt <- mean_dfDbt %>% 
          filter(time %in% daytime) %>%   #### in mean method ,use day2 dbt as strain doubling_time 
          group_by_at(c("strains","biorep","time","panel",
                        Categorical_variable[!Categorical_variable %in% 
                                               c("rep_location_name","rep_location_num","tube_type")])) %>%
          dplyr::summarize(doubling_time = mean(doubling_time)) %>%
          group_by_at(c("strains","time","panel",
                        Categorical_variable[!Categorical_variable %in% 
                                               c("rep_location_name","rep_location_num","tube_type")])) %>%
          dplyr::summarize(doubling_time = mean(doubling_time))
        return(mean_Dbt)
      }else{
        mean_Dbt <- mean_dfDbt %>% 
          group_by_at(c("strains","biorep","time","panel",
                        Categorical_variable[!Categorical_variable %in% 
                                               c("rep_location_name","rep_location_num","tube_type")])) %>%
          dplyr::summarize(doubling_time = mean(doubling_time)) %>%
          group_by_at(c("strains","time","panel",
                        Categorical_variable[!Categorical_variable %in% 
                                               c("rep_location_name","rep_location_num","tube_type")]))  %>%
          dplyr::summarize(sd_bio = sd(doubling_time),doubling_time = mean(doubling_time))
        ####contain other days data ,this sd is biorep sd
        ###be caution: the order of new colum must be sd_bio first and the follow is doubling_time
        return(mean_Dbt)
      }
    }else{
      if( sd == FALSE){
        median_Dbt <- dfDbt_ALL %>%
          filter(row_number() == 2) %>% 
          ungroup()%>% 
          select(c(-variable,-value)) %>% 
          group_by_at(c("strains","biorep","time","panel",
                        Categorical_variable[!Categorical_variable %in% 
                                               c("rep_location_name","rep_location_num","tube_type")])) %>%
          dplyr::summarize(doubling_time = mean(doubling_time)) %>% 
          group_by_at(c("strains","time","panel",
                        Categorical_variable[!Categorical_variable %in% 
                                               c("rep_location_name","rep_location_num","tube_type")])) %>%
          dplyr::summarize(doubling_time = mean(doubling_time)) %>%
          group_by_at(c("strains",Categorical_variable[!Categorical_variable %in% 
                                                         c("rep_location_name","rep_location_num","tube_type")])) %>% 
          dplyr::summarize(doubling_time = median(doubling_time)) 
        return(median_Dbt)
      }else{
        median_Dbt <- dfDbt_ALL  %>%
          filter(row_number() == 2) %>% 
          ungroup()%>% 
          select(c(-variable,-value)) %>% 
          group_by_at(c("strains","biorep","time","panel",
                        Categorical_variable[!Categorical_variable %in% 
                                               c("rep_location_name","rep_location_num","tube_type")])) %>%
          dplyr::summarize(doubling_time = mean(doubling_time)) %>% 
          group_by_at(c("strains","time","panel",Categorical_variable[!Categorical_variable %in% c("rep_location_name","rep_location_num","tube_type")])) %>%
          dplyr::summarize(sd_bio = sd(doubling_time),doubling_time = mean(doubling_time))
        ####contain other days data ,this sd is biorep sd##
        return(median_Dbt)
      }
    }
  }
  
}

####clean Delta data by ODoutlier_data
Delta_clean <- function(Delta_data,ODoutlier_data){
  Delta_clean_data <- anti_join(Delta_data,ODoutlier_data,by = names(ODoutlier_data)[-c(13,14,15,16,17)])
  return(Delta_clean_data)
}

###calculate timeshift##
###function: calculate timeshift that shiftOD value is not NA
### input: Delta600Data for calculate input
###        default_rm_vector is a vector used to remove some strains you did not want
###        Categorical_variable is the vector for group , melt and dcast function
###        argument: shiftOD that that choose the OD cutoff
### output: dfTimeshift
timeshift_calculation <- function(Delta600Data,default_rm_vector=NA,Categorical_variable=NULL,shiftOD = 0.3){
  dfTimeshift <- Delta600Data %>% filter(!strains %in% default_rm_vector) %>%
    melt(id.var=c("strains","biorep","techrep","time","panel",Categorical_variable) )%>%
    filter(strains != "blank") %>%
    group_by_at(c("strains","biorep","techrep","time","panel",Categorical_variable)) %>%
    dplyr::mutate( t_offset = 10 - 10 * (log(shiftOD) - log(c(NA,value[-length(value)]))) /  c(NA,diff(log(value))) ) %>%
    filter(value >= shiftOD) %>%
    filter(row_number() == 1) %>%
    mutate(t_shiftOD = as.numeric(substring(variable,2)) - t_offset ) 
  
  return(dfTimeshift)
}

###function: calculate timeshift that shiftOD is NA value and bind all the timeshift data
###input: Delta600Data for calculate input and search the shiftOD NA value
###       dfTimeshift that contain the calculated data and use for setdiff(anti_join)
###       doubling_timeData is the data from function "doubling_time_calculation", 
###       contain the strains doubling time for estimate the timeshift to the shiftOD
###       default_rm_vector is a vector used to remove some strains you did not want
###       Categorical_variable is the vector for group , melt and dcast function
###       argument: shiftOD that choose the OD cutoff
timeshift_NA_calculate_and_bind_data <- function(Delta600Data,dfTimeshift,doubling_timeData,
                                                 default_rm_vector = NA,Categorical_variable=NULL,shiftOD = 0.3){
  dfAllTimeshift_element <-  Delta600Data %>% filter(!strains %in% default_rm_vector)
  
  dfNATimeshift_data <- anti_join(dfAllTimeshift_element,dfTimeshift,
                                  by=c("strains","biorep","techrep","time","panel",Categorical_variable))
  
  dfNATimeshift_OD <- dfNATimeshift_data %>% 
    filter(!strains %in% default_rm_vector) %>%
    melt(id.var = c("strains","biorep","techrep","time","panel",Categorical_variable)) %>%
    mutate(variable = as.numeric(substr(variable,2,9999))) %>%
    group_by_at(c("strains","biorep","techrep","time","panel",
                  Categorical_variable[!Categorical_variable %in% 
                                         c("rep_location_name","rep_location_num","tube_type")])) %>%
    arrange(-variable) %>%
    filter(!is.na(value)) %>% 
    filter(row_number()==1)%>%
    mutate(value = ifelse(value <= 0 , 0.0001,value)) %>% 
    ungroup()
  
  dfNATimeshift <- dfNATimeshift_OD %>%
    merge(doubling_timeData,by=c("strains",
                                 Categorical_variable[!Categorical_variable %in% 
                                                        c("rep_location_name","rep_location_num","tube_type")])) %>%
    mutate(t_shiftOD = (log(shiftOD)-log(value))*doubling_time/log(2)+variable) %>%
    mutate(variable = paste("t",variable,sep = "")) %>%
    select(-doubling_time) %>%
    select(strains,biorep,techrep,time,panel,one_of(Categorical_variable),everything())
  
  dfAllTimeshift <- rbind(dfTimeshift %>% select(-t_offset) %>% ungroup(),dfNATimeshift)
  
  return(dfAllTimeshift)
}

###survival rate calculation
###function: survival rate calculation
###input: doubling_timeData is the strain's doubling time data for survival rate calculate,
###       from function "doubling_time_calculation"
###       AllTimeshiftData is the strain's time shift in different day for survival rate calculate,
###       from function "timeshift_calculation" and "timeshift_NA_calculate_and_bind_data"
###       default_rm_vector is a vector used to remove some strains you did not want
###       Categorical_variable is the vector for group , melt and dcast function
###       Relative_day is the day choosed to calculate survival

survival_rate_calculate <- function(doubling_timeData,AllTimeshiftData,default_rm_vector = NA,Relative_day = "day2",
                                    Categorical_variable=NULL){
  AllTimeShiftData_tech_mean <- AllTimeshiftData %>%
    group_by_at(c("strains","biorep","time","panel",
                  Categorical_variable[!Categorical_variable %in% 
                                         c("rep_location_name","rep_location_num","tube_type")])) %>%
    dplyr::summarize(sd_timeshift_tech = sd(t_shiftOD),mean_t_shiftOD = mean(t_shiftOD)) %>% 
    group_by_at(c("strains","biorep","time","panel",
                  Categorical_variable[!Categorical_variable %in% c("rep_location_name","rep_location_num","tube_type")]))
  ###survival rate calculate
  SurvivalRate_retain_bio <- AllTimeShiftData_tech_mean %>%
    left_join(AllTimeShiftData_tech_mean %>% filter(time == Relative_day),
              by = group_vars(AllTimeShiftData_tech_mean)[-c(3,4)]) %>%
    filter(!is.na(time.y)) %>% 
    ###SRver5-1 reconstruction function : filter Relative_day survival rate IS NA data , prevent from SI calculate BUG! 
    mutate(delta_t = as.numeric(mean_t_shiftOD.x - mean_t_shiftOD.y)) %>%
    select(strains,biorep,time = time.x ,panel = panel.x,
           Categorical_variable[!Categorical_variable %in% c("rep_location_name","rep_location_num","tube_type")],
           sd_timeshift_tech = sd_timeshift_tech.x,mean_t_shiftOD = mean_t_shiftOD.x,delta_t) %>%
    left_join(doubling_timeData,by = c("strains",
                                       Categorical_variable[!Categorical_variable %in% 
                                                              c("rep_location_name","rep_location_num","tube_type")])) %>%
    mutate(Survival_n = 1 / 2^(delta_t / doubling_time)) 
  
  SurvivalRate <- SurvivalRate_retain_bio %>%
    group_by_at(c("strains","time","panel",Categorical_variable[!Categorical_variable %in% c("rep_location_name","rep_location_num","tube_type")])) %>%
    dplyr::summarize(Survival_n_mean_bio = mean(Survival_n),sd_Survival_n_bio = sd(Survival_n))
  
  return(SurvivalRate_retain_bio)
}

###SI calculate 
###function: SI calculation
###input:SurvivalRateData is the strain's survival data from function "survival_rate_calculate"
###       default_rm_vector is a vector used to remove some strains you did not want
###       Categorical_variable is the vector for group , melt and dcast function
###       Relative_day is the day choosed to calculate survival (the same as survival_rate_calculate function)
###       Remove_not_enough_datapoint is the argument for control the biorep 
###       which do not have enough(6 data point) 20220922 debug,
###       datapoint_number is the argument for control data point number 
###       when the subset SR data not contain the enough data point 20220923 debug,
SurvivalIntegral_calculate <- function(SurvivalRateData,default_rm_vector = NA,Relative_day="day2",
                                       Categorical_variable = NULL,retain_biorep = F,
                                       Remove_not_enough_datapoint=T,datapoint_number=6){
  #RM data point not contain 6 data point
  if(Remove_not_enough_datapoint==T){
    SurvivalRateData <- SurvivalRateData %>%
      group_by_at(c("strains","biorep",
                    Categorical_variable[!Categorical_variable %in% 
                                           c("rep_location_name","rep_location_num","tube_type")])) %>%
      filter(n()==datapoint_number)
  }
  if(retain_biorep == T){
    SurvivalIntegral_biorep <- SurvivalRateData %>% 
      filter(!strains %in% default_rm_vector) %>%
      group_by_at(c("strains","biorep",
                    Categorical_variable[!Categorical_variable %in% 
                                           c("rep_location_name","rep_location_num","tube_type")])) %>%
      mutate(num_time = as.numeric(gsub("day(\\d+)","\\1",time))) %>%
      filter(num_time >= as.numeric(gsub("day(\\d+)","\\1",Relative_day))) %>% 
      ###SRver5-1 reconstruction function : filter >= Relative_day survival rate data
      dplyr::arrange(num_time) %>%
      mutate(SI = c(NA,(Survival_n[-length(Survival_n)]+Survival_n[-1])/2 * diff(num_time))) %>%
      dplyr::summarize(SI = sum(SI,na.rm = T))
    
    return(SurvivalIntegral_biorep)
  }else{
    SurvivalIntegral <-  SurvivalRateData %>% 
      filter(!strains %in% default_rm_vector) %>%
      group_by_at(c("strains","biorep",
                    Categorical_variable[!Categorical_variable %in% 
                                           c("rep_location_name","rep_location_num","tube_type")])) %>%
      mutate(num_time = as.numeric(gsub("day(\\d+)","\\1",time))) %>%
      filter(num_time >= as.numeric(gsub("day(\\d+)","\\1",Relative_day))) %>% 
      ###SRver5-1 reconstruction function : filter >= Relative_day survival rate data
      dplyr::arrange(num_time) %>%
      mutate(SI = c(NA,(Survival_n[-length(Survival_n)]+Survival_n[-1])/2 * diff(num_time))) %>%
      dplyr::summarize(SI = sum(SI,na.rm = T)) %>%
      group_by_at(c("strains",Categorical_variable[!Categorical_variable %in% 
                                                     c("rep_location_name","rep_location_num","tube_type")])) %>%
      dplyr::summarize(sd_SI = sd(SI),SI_mean_bio = mean(SI))
    
    return(SurvivalIntegral)
  }
}

