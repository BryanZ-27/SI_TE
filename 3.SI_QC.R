###quality control of survival integral
#library
library(reshape2)
library(stringr)
library(dplyr)
library(doParallel)

#source function
source("./3_1.SRver5-1_function.R")
source("./3_2.SI_biorep_clean_function.R")
source("./3_3.SI_remove_wildtype_and_outliers_function.R")

#main
####load data####
load("./3.raw_survival_data.Rdata")

####find techrep outlier, first time####
ODoutlier_1_data <- find_techrep_outlier(survivalData$Blank600,default_rm_vector =NA,CV_cutoff = 0.2,OD_time = "t50",
                                         Categorical_variable=Categorical_variable)
Blank600_clean_data <- Blank_clean(Blank600_data = survivalData$Blank600,ODoutlier_data = ODoutlier_1_data)

####find techrep outlier, second time####
ODoutlier_2_data <- find_techrep_outlier(Blank600_clean_data,default_rm_vector =NA,CV_cutoff = 0.2,OD_time = "t50",
                                         Categorical_variable=Categorical_variable)

#bind_rows ODoutlier data to get techrep outlier in first time and second time #####
ODoutlier_data <- bind_rows(ODoutlier_1_data,ODoutlier_2_data)

##clean techrep outlier before doubling time calculation
LnBlank_clean_data <- LnBlank_clean(LnBlank_data = survivalData$LnBlank,ODoutlier_data = ODoutlier_data)

##clean techrep outlier before timeshift calculation
Delta_clean_data <- Delta_clean(Delta_data = survivalData$Delta600,ODoutlier_data = ODoutlier_data)

# add doubling_time calculation new function.
doubling_time_calculation_day2_min <- function(LnBlankData,default_rm_vector = NA, Categorical_variable = NULL,
                                               OD_interval = c(0.2,0.5),well = FALSE, retain_biorep = FALSE, 
                                               is.mean = FALSE,is.min=FALSE,daytime = c("day2") , sd = FALSE , 
                                               is.all = F){
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
  if (is.mean == T && is.min == F) {
    mean_dfDbt <- dfDbt_ALL %>% 
      dplyr::summarize(doubling_time = mean(doubling_time)) %>%
      ungroup()
  }else if(is.min == T){
    min_dfDbt <- dfDbt_ALL %>% 
      filter(row_number() == 1) %>% 
      ungroup() %>%
      select(c(-variable,-value)) 
  }else{
    dfDbt <- dfDbt_ALL %>% 
      filter(row_number() == 2) %>% 
      ungroup() %>%
      select(c(-variable,-value)) 
  }
  ####output type: well,biorep,strain
  if (well == T) {
    if (is.mean == T && is.min == F) {   
      return(mean_dfDbt)
    }else if(is.min == T){
      return(min_dfDbt)
    }else{
      return(dfDbt)
    }
  }else if(retain_biorep == T){
    if(is.mean == T && is.min == F){
      if (sd == FALSE) {
        mean_dfDbt_contain_biorep <-  mean_dfDbt %>% 
          filter(time %in% daytime) %>%   #### in mean method ,use day2 dbt as strain doubling_time 
          group_by_at(c("strains","biorep","time","panel",Categorical_variable[!Categorical_variable %in% c("rep_location_name","rep_location_num","tube_type")])) %>%
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
      
    }else if(is.min == T){
      if (sd == FALSE) {
        min_dfDbt_contain_biorep <-  dfDbt_ALL %>% 
          filter(time %in% daytime) %>%
          filter(row_number() == 1) %>% 
          ungroup()%>% 
          select(c(-variable,-value)) %>% 
          group_by_at(c("strains","biorep","time","panel",
                        Categorical_variable[!Categorical_variable %in% 
                                               c("rep_location_name","rep_location_num","tube_type")])) %>%
          dplyr::summarize(doubling_time = mean(doubling_time))  
        return(min_dfDbt_contain_biorep)
      }else{
        min_dfDbt_contain_biorep <- dfDbt_ALL %>% 
          # filter(time %in% daytime) %>% 
          filter(row_number() == 1) %>% 
          ungroup()%>% 
          select(c(-variable,-value)) %>% 
          group_by_at(c("strains","biorep","time","panel",
                        Categorical_variable[!Categorical_variable %in% 
                                               c("rep_location_name","rep_location_num","tube_type")])) %>% 
          dplyr::summarize(sd_tech = sd(doubling_time),
                           doubling_time = mean(doubling_time))  ####contain other days data ,this sd is techrep sd##
        return(min_dfDbt_contain_biorep)
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
          group_by_at(c("strains","biorep","time","panel",
                        Categorical_variable[!Categorical_variable %in% 
                                               c("rep_location_name","rep_location_num","tube_type")])) %>% 
          dplyr::summarize(sd_tech = sd(doubling_time),
                           doubling_time = mean(doubling_time))  ####contain other days data ,this sd is techrep sd##
        return(dfDbt_contain_biorep)
      }
    }
  }else{
    if(is.mean == T && is.min == F){
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
          dplyr::summarize(sd_bio = sd(doubling_time),
                           doubling_time = mean(doubling_time))  ####contain other days data ,this sd is biorep sd
        ###be caution: the order of new colum must be sd_bio first and the follow is doubling_time
        return(mean_Dbt)
      }
    }else if(is.min == T){
      if( sd == FALSE){
        day2_min_Dbt <- dfDbt_ALL %>%
          filter(time %in% daytime) %>%
          filter(row_number() == 1) %>% 
          ungroup()%>% 
          select(c(-variable,-value)) %>% 
          group_by_at(c("strains","biorep","time","panel",
                        Categorical_variable[!Categorical_variable %in% 
                                               c("rep_location_name","rep_location_num","tube_type")])) %>%
          dplyr::summarize(doubling_time = mean(doubling_time)) %>% 
          group_by_at(c("strains","time","panel",
                        Categorical_variable[!Categorical_variable %in% 
                                               c("rep_location_name","rep_location_num","tube_type")])) %>%
          dplyr::summarize(doubling_time = mean(doubling_time))
        return(day2_min_Dbt)
      }else{
        min_Dbt <- dfDbt_ALL  %>%
          # filter(time %in% daytime) %>%
          filter(row_number() == 1) %>% 
          ungroup()%>% 
          select(c(-variable,-value)) %>% 
          group_by_at(c("strains","biorep","time","panel",
                        Categorical_variable[!Categorical_variable %in% 
                                               c("rep_location_name","rep_location_num","tube_type")])) %>%
          dplyr::summarize(doubling_time = mean(doubling_time)) %>% 
          group_by_at(c("strains","time","panel",
                        Categorical_variable[!Categorical_variable %in% 
                                               c("rep_location_name","rep_location_num","tube_type")])) %>%
          dplyr::summarize(sd_bio = sd(doubling_time),
                           doubling_time = mean(doubling_time))   ####contain other days data ,this sd is biorep sd##
        return(min_Dbt)
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
          group_by_at(c("strains","time","panel",
                        Categorical_variable[!Categorical_variable %in% 
                                               c("rep_location_name","rep_location_num","tube_type")])) %>%
          dplyr::summarize(sd_bio = sd(doubling_time),
                           doubling_time = mean(doubling_time))   ####contain other days data ,this sd is biorep sd##
        return(median_Dbt)
      }
    }
  } 
}

#main
#doubling time calculation,median type & day2 mean type & day2 min type, 
#summary to strain for timeshift and survival rate calculate,contatin strains_group,tube_type,glu....Categorical_variable
#1.Argument vector
phe_name="doubling_time"
calculate_type_vector=c("day2_min","day2_mean","median")
OD_interval_vector=c("0.2_0.5","0.25_0.75","0.3_0.8")

#2.Create doubling_time_argument_data
doubling_time_argument_data <- expand.grid(phe_name=phe_name,calculate_type=calculate_type_vector,
                                           OD_interval=OD_interval_vector,stringsAsFactors = F) %>% 
  mutate(rowname=paste(phe_name,calculate_type,OD_interval,sep = "_")) %>%
  tibble::column_to_rownames(var = "rowname")
#3.According to doubling_time_argument_data, split as list (one row one list)
#and use mclapply to calculate different arugment doubling time.
doubling_time_output_list <- mclapply(split(doubling_time_argument_data,f = rownames(doubling_time_argument_data)),
                                      FUN=function(argument_data){
  phe_name=argument_data$phe_name
  calculate_type=argument_data$calculate_type
  OD_interval=as.numeric(str_split(argument_data$OD_interval,pattern = "_",simplify = T))
  if(calculate_type == "day2_min"){
    cat(calculate_type)
    output <- doubling_time_calculation_day2_min(LnBlankData = LnBlank_clean_data,default_rm_vector = "empty",
                                                 Categorical_variable = Categorical_variable,
                                                 well = F,retain_biorep = F,is.mean = F,is.min = T,daytime = "day2",
                                                 sd = F,is.all = F,OD_interval =  OD_interval)%>%
      ungroup %>% select(c(-time,-panel))
    
  }else if(calculate_type == "day2_mean"){
    cat(calculate_type)
    output <- doubling_time_calculation_day2_min(LnBlankData = LnBlank_clean_data,default_rm_vector = "empty",
                                                 Categorical_variable = Categorical_variable,
                                                 well = F,retain_biorep = F,is.mean = T,is.min = F,daytime = "day2",
                                                 sd = F,is.all = F,OD_interval =  OD_interval)%>%
      ungroup %>% select(c(-time,-panel))
  }else if(calculate_type == "median"){
    cat(calculate_type)
    output <- doubling_time_calculation_day2_min(LnBlankData = LnBlank_clean_data,default_rm_vector = "empty",
                                                 Categorical_variable = Categorical_variable,
                                                 well = F,retain_biorep = F,is.mean = F,is.min = F,daytime = "day2",
                                                 sd = F,is.all = F,OD_interval = OD_interval )
  }else{stop("ERROR WITH NOT CONTAIN FUNCTION TYPE")}
  return(output)
},mc.cores = 10)

#TIMESHIFT calculation PART######################
##Time shift calculation in the range(shift_OD = 0.3,0.4,0.5,0.6)
#1.Argument vector
phe_name="dfTimeshift"
shift_OD_vector=c(0.3,0.4,0.5,0.6)

#2.Create dfTimeshift_argument_data
dfTimeshift_argument_data <- expand.grid(phe_name=phe_name,shift_OD=shift_OD_vector,stringsAsFactors = F) %>% 
  mutate(rowname=paste(phe_name,shift_OD,sep = "_")) %>%
  tibble::column_to_rownames(var = "rowname")

#3.According to dfTimeshift_argument_data, split as list (one row one list)and use mclapply to calculate different argument Timeshift.
dfTimeshift_output_list <- mclapply(split(dfTimeshift_argument_data,f = rownames(dfTimeshift_argument_data)),
                                    FUN = function(argument_data){
  shift_OD = argument_data$shift_OD
  output <- timeshift_calculation(Delta600Data = Delta_clean_data,default_rm_vector = "empty",
                                  Categorical_variable = Categorical_variable,shiftOD = shift_OD)
  return(output)
},mc.cores = 10)

##Time shift calculation that use for doubling time to predict and bind all time shift data
#1.Create AllTimeshift_argument_data
AllTimeshift_argument_data <- expand.grid(doubling_time_argument=names(doubling_time_output_list),
                                          dfTimeshift_argument=names(dfTimeshift_output_list),stringsAsFactors = F) %>%
  mutate(rowname=paste("AllTimeshift",gsub("doubling_time_(.*)","\\1",doubling_time_argument),
                       gsub("dfTime(.*)","\\1",dfTimeshift_argument),sep = "_")) %>%
  tibble::column_to_rownames(var = "rowname")

#2.According to AllTimeshift_argument_data, split as list (one row one list)
#and use mclapply to calculate different argument AllTimeshift.
AllTimeshift_output_list <- mclapply(split(AllTimeshift_argument_data,f = rownames(AllTimeshift_argument_data)) ,
                                     FUN = function(argument_data){
  doubling_time_data_name <- argument_data$doubling_time_argument
  dfTimeshift_data_name <-  argument_data$dfTimeshift_argument
  shift_OD <- as.numeric(str_extract(dfTimeshift_data_name,pattern = "\\d+\\.\\d+"))
  output <- timeshift_NA_calculate_and_bind_data(Delta600Data = Delta_clean_data,
                                                 dfTimeshift = dfTimeshift_output_list[[dfTimeshift_data_name]],
                                                 doubling_timeData = doubling_time_output_list[[doubling_time_data_name]],
                                                 default_rm_vector = "empty",Categorical_variable = Categorical_variable,
                                                 shiftOD = shift_OD)
  return(output)
},mc.cores = 20)

##Survival Rate calculate PART e.g day2 Relative_day##############
#1.Create SurvivalRate_argument_data
SurvivalRate_argument_data <- AllTimeshift_argument_data %>% 
  tibble::rownames_to_column() %>%
  dplyr::rename(AllTimeshift_argument=rowname) %>%
  left_join(expand.grid(AllTimeshift_argument=rownames(AllTimeshift_argument_data),
                        Relative_day=c("day2","day4","day6","day9","day12"),stringsAsFactors = F)) %>%
  mutate(rowname=paste("SurvivalRate",gsub("AllTimeshift_(.*)","\\1",AllTimeshift_argument),Relative_day,sep = "_")) %>%
  tibble::column_to_rownames(var = "rowname")

#2.According to SurvivalRate_argument_data,split as list (one row one list)
#and use mclapply to calculate different argument SurvivalRate.
SurvivalRate_output_list <- mclapply(split(SurvivalRate_argument_data,f = rownames(SurvivalRate_argument_data)),
                                     FUN = function(argument_data){
  doubling_time_data_name <- argument_data$doubling_time_argument
  AllTimeshift_data_name <-  argument_data$AllTimeshift_argument
  Relative_day=argument_data$Relative_day
  output <-  survival_rate_calculate(doubling_timeData = doubling_time_output_list[[doubling_time_data_name]],
                                     AllTimeshiftData = AllTimeshift_output_list[[AllTimeshift_data_name]],
                                     default_rm_vector = "empty",Categorical_variable = Categorical_variable,
                                     Relative_day = Relative_day)
  return(output)
},mc.cores = 30 )

##Survival Integral calculate e.g retain biorep######
#1.create SurvivalIntegral_argument_data
SurvivalIntegral_argument_data <- SurvivalRate_argument_data %>%
  tibble::rownames_to_column() %>%
  dplyr::rename(SurvivalRate_argument=rowname) %>%
  mutate(rowname1=paste("SurvivalIntegral_for_biorep",gsub("SurvivalRate_(.*)","\\1",SurvivalRate_argument),sep = "_")) %>%
  mutate(rowname2=paste("SurvivalIntegral_for_strain",gsub("SurvivalRate_(.*)","\\1",SurvivalRate_argument),sep = "_")) 

SurvivalIntegral_for_biorep_argument_data <- SurvivalIntegral_argument_data %>% 
  tibble::column_to_rownames(var = "rowname1")
SurvivalIntegral_for_strain_argument_data <- SurvivalIntegral_argument_data %>% 
  tibble::column_to_rownames(var = "rowname2")

#2.According to SurvivalIntegral_argument_data,split as list (one row one list)
SurvivalIntegral_for_biorep_output_list <- mclapply(split(SurvivalIntegral_for_biorep_argument_data,
                                                          f = rownames(SurvivalIntegral_for_biorep_argument_data)),
                                                    FUN = function(argument_data){
  SurvivalRate_data_name <- argument_data$SurvivalRate_argument
  Relative_day=argument_data$Relative_day
  output <-  SurvivalIntegral_calculate(SurvivalRateData = SurvivalRate_output_list[[SurvivalRate_data_name]],
                                        default_rm_vector = "empty",Categorical_variable = Categorical_variable,
                                        retain_biorep = T,Relative_day = Relative_day,Remove_not_enough_datapoint = F)
  return(output)
},mc.cores = 30 )

SurvivalIntegral_for_strain_output_list <- mclapply(split(SurvivalIntegral_for_strain_argument_data,
                                                          f = rownames(SurvivalIntegral_for_strain_argument_data)),
                                                    FUN = function(argument_data){
  SurvivalRate_data_name <- argument_data$SurvivalRate_argument
  Relative_day=argument_data$Relative_day
  output <-  SurvivalIntegral_calculate(SurvivalRateData = SurvivalRate_output_list[[SurvivalRate_data_name]],
                                        default_rm_vector = "empty",
                                        Categorical_variable = Categorical_variable,retain_biorep = F,
                                        Relative_day = Relative_day,Remove_not_enough_datapoint = F)
  return(output)
},mc.cores = 30 )

# MAIN#######
#1.According to argument_data, pull the data which will be calculated######
#In this script , will choose different doubling_time, different_shift_OD and Relative_day=day2, 
#SR data for following calculate
Argument_data_subset <- SurvivalDataAnalysised$Argument_data %>% filter(Relative_day=="day2")

#2.Calculate SI_for_biorep_ori for following biorep control.####
# Use mclapply to calculate each SI_for_biorep_ori data
SI_for_biorep_ori_Argument_data <- Argument_data_subset %>%
  mutate(SI_for_biorep_ori_name=paste("SI_for_biorep_ori",
                                      gsub("SurvivalRate_(.*)","\\1",SurvivalRate_argument),sep = "_")) %>%
  tibble::column_to_rownames(var = "SI_for_biorep_ori_name")

SI_for_biorep_ori_output_list <- mclapply(split(SI_for_biorep_ori_Argument_data,
                                                f = rownames(SI_for_biorep_ori_Argument_data)),
                                          FUN = function(argument_data){
  SurvivalRate_data_name <- argument_data$SurvivalRate_argument
  Relative_day = argument_data$Relative_day
  output <- SurvivalIntegral_calculate(
    SurvivalRateData = SurvivalDataAnalysised$SurvivalRate_data_list[[SurvivalRate_data_name]],default_rm_vector = NA,
    Relative_day = Relative_day,Categorical_variable = Categorical_variable,retain_biorep = T,
    Remove_not_enough_datapoint = T,datapoint_number = 6)
  return(output)
},mc.cores = 30)

#3.Clean biorep by sd SI , CV < 0.2 cutoff ,first time#####
SI_for_biorep_ori_outlier_1_output_list <- mclapply(SI_for_biorep_ori_output_list,FUN = function(SI_for_biorep_ori_data){
  SI_for_biorep_ori_outlier_1_data <- find_SI_biorep_outlier(SIData_contain_bio = SI_for_biorep_ori_data,CV_cutoff = 0.2,
                                                             Categorical_variable = Categorical_variable)
  return(SI_for_biorep_ori_outlier_1_data)
},mc.cores = 30)

#4.According to step 3 , remove the biorep from SI_for_biorep_ori_output_list#######
SI_for_biorep_ori_clean_1_output_list <- mclapply(names(SI_for_biorep_ori_output_list),
                                                  FUN = function(SI_for_biorep_ori_data_name){
  SI_for_biorep_ori_data <- SI_for_biorep_ori_output_list[[SI_for_biorep_ori_data_name]]
  SI_for_biorep_ori_outlier_1_data <- SI_for_biorep_ori_outlier_1_output_list[[SI_for_biorep_ori_data_name]]
  SI_for_biorep_ori_clean_1_data <- SI_for_biorep_ori_data %>% anti_join(SI_for_biorep_ori_outlier_1_data)
  return(SI_for_biorep_ori_clean_1_data)
},mc.cores = 30)
names(SI_for_biorep_ori_clean_1_output_list) <- names(SI_for_biorep_ori_output_list)

#5.Clean biorep by sd SI , CV < 0.2 cutoff ,second time#####
SI_for_biorep_ori_outlier_2_output_list <- mclapply(SI_for_biorep_ori_clean_1_output_list,
                                                    FUN = function(SI_for_biorep_ori_clean_1_data){
  SI_for_biorep_ori_outlier_2_data <- find_SI_biorep_outlier(SIData_contain_bio = SI_for_biorep_ori_clean_1_data,
                                                             CV_cutoff = 0.2,Categorical_variable = Categorical_variable)
  return(SI_for_biorep_ori_outlier_2_data)
},mc.cores = 30)

# bind the outlier data of the first & second time 
SI_for_biorep_ori_outlier_Argument_data <- SI_for_biorep_ori_Argument_data %>% 
  tibble::rownames_to_column() %>%
  dplyr::rename(SI_for_biorep_ori_name = rowname) %>%
  mutate(SI_for_biorep_ori_outlier_name=paste("SI_for_biorep_ori_outlier",
                                              gsub("SI_for_biorep_ori_(.*)","\\1",SI_for_biorep_ori_name),sep = "_")) %>%
  tibble::column_to_rownames(var = "SI_for_biorep_ori_outlier_name")

SI_for_biorep_ori_outlier_output_list <- mclapply(split(SI_for_biorep_ori_outlier_Argument_data,
                                                        f = rownames(SI_for_biorep_ori_outlier_Argument_data)),
                                                  FUN = function(argument_data){
  SI_for_biorep_ori_outlier_list_name <- argument_data$SI_for_biorep_ori_name
  output <-  bind_rows(SI_for_biorep_ori_outlier_1_output_list[[SI_for_biorep_ori_outlier_list_name]],
                       SI_for_biorep_ori_outlier_2_output_list[[SI_for_biorep_ori_outlier_list_name]])
  return(output)
},mc.cores = 30 )

#6.According to SI_for_biorep_ori_outlier_output_list data , remove the correspond SR data ####
SurvivalRate_clean_SI_biorep_ori_Argument_data <- SI_for_biorep_ori_outlier_Argument_data %>%
  tibble::rownames_to_column() %>%
  dplyr::rename(SI_for_biorep_ori_outlier_name = rowname) %>%
  mutate(SurvivalRate_clean_SI_biorep_ori_name = paste("SurvivalRate_clean_SI_biorep_ori",
                                                       gsub("SI_for_biorep_ori_(.*)","\\1",SI_for_biorep_ori_name),
                                                       sep = "_")) %>%
  tibble::column_to_rownames(var = "SurvivalRate_clean_SI_biorep_ori_name")

SurvivalRate_clean_SI_biorep_ori_output_list <- mclapply(
  split(SurvivalRate_clean_SI_biorep_ori_Argument_data,f = rownames(SurvivalRate_clean_SI_biorep_ori_Argument_data)),
  FUN = function(argument_data){
  SurvivalRate_data_name = argument_data$SurvivalRate_argument
  SI_for_biorep_ori_outlier_name = argument_data$SI_for_biorep_ori_outlier_name
  output <- SurvivalDataAnalysised$SurvivalRate_data_list[[SurvivalRate_data_name]] %>% 
    anti_join(SI_for_biorep_ori_outlier_output_list[[SI_for_biorep_ori_outlier_name]])
  return(output)
} ,mc.cores = 30)

#7.Use subset daytime SR data ,calculate pre_subset_SI_for_biorep_data,
#for rm similar SI value but different SR curve biorep, here use day2-day9 subset data for calculate SI.####
Subset_SurvivalRate_clean_SI_biorep_ori_day2_day9_output_list <- mclapply(
  SurvivalRate_clean_SI_biorep_ori_output_list,FUN = function(SurvivalRate_clean_SI_biorep_ori_data){
  Subset_SR_data <- SurvivalRate_clean_SI_biorep_ori_data %>%
    filter(!time %in% c("day12","day15"))
  return(Subset_SR_data)
},mc.cores = 30 )

#8.Calculate pre_subset_SI_data####
Subset_SI_for_biorep_clean_SI_biorep_ori_day2_day9_output_list <- mclapply(
  Subset_SurvivalRate_clean_SI_biorep_ori_day2_day9_output_list,
  FUN = function(Subset_SurvivalRate_clean_SI_biorep_ori_data){
  Subset_SI_for_biorep_data <- SurvivalIntegral_calculate(
    SurvivalRateData = Subset_SurvivalRate_clean_SI_biorep_ori_data ,default_rm_vector = NA,
    Relative_day = "day2",Categorical_variable = Categorical_variable,retain_biorep = T,Remove_not_enough_datapoint = T,
    datapoint_number = 4)
  return(Subset_SI_for_biorep_data)
},mc.cores = 30)

#9.find SR curve outlier#######
SI_for_biorep_curve_outlier_Argument_data <- SurvivalRate_clean_SI_biorep_ori_Argument_data %>%
  tibble::rownames_to_column() %>%
  dplyr::rename(SurvivalRate_clean_SI_biorep_ori_name=rowname) %>%
  mutate(SI_for_biorep_curve_outlier_name = paste(
    "SI_for_biorep_curve_outlier",gsub("SI_for_biorep_ori_(.*)","\\1",SI_for_biorep_ori_name),sep = "_")) %>%
  tibble::column_to_rownames(var = "SI_for_biorep_curve_outlier_name")

SI_for_biorep_curve_outlier_output_list <- mclapply(split(SI_for_biorep_curve_outlier_Argument_data,
                                                          f = rownames(SI_for_biorep_curve_outlier_Argument_data)),
                                                    FUN = function(argument_data){
  Subset_SI_for_biorep_clean_SI_biorep_ori_data_name <- argument_data$SurvivalRate_clean_SI_biorep_ori_name 
  output <- find_SI_biorep_outlier(
    SIData_contain_bio = 
      Subset_SI_for_biorep_clean_SI_biorep_ori_day2_day9_output_list[[Subset_SI_for_biorep_clean_SI_biorep_ori_data_name]],
    CV_cutoff = 0.25,Categorical_variable = Categorical_variable)
  return(output)
},mc.cores = 30 )

#10.According to SI_biorep_outlier & SI_biorep_outlier_for_curve , clean SR data for final SI calculation######
SurvivalRate_clean_SI_biorep_ori_and_curve_Argument_data <- SI_for_biorep_curve_outlier_Argument_data %>%
  tibble::rownames_to_column() %>%
  dplyr::rename(SI_for_biorep_curve_outlier_name=rowname) %>%
  mutate(SurvivalRate_clean_SI_biorep_ori_and_curve_data_name = paste("SurvivalRate_clean_SI_biorep_ori_and_curve",
                                                                      gsub("SI_for_biorep_ori_(.*)","\\1",
                                                                           SI_for_biorep_ori_name),sep = "_")) %>%
  tibble::column_to_rownames(var = "SurvivalRate_clean_SI_biorep_ori_and_curve_data_name")

SurvivalRate_clean_SI_biorep_ori_and_curve_output_list <- mclapply(split(
  SurvivalRate_clean_SI_biorep_ori_and_curve_Argument_data,f = rownames(
    SurvivalRate_clean_SI_biorep_ori_and_curve_Argument_data)),FUN = function(argument_data){
  SurvivalRate_data_name = argument_data$SurvivalRate_argument
  SI_for_biorep_ori_outlier_name = argument_data$SI_for_biorep_ori_outlier_name
  SI_for_biorep_curve_outlier_name=argument_data$SI_for_biorep_curve_outlier_name
  
  output <- SurvivalDataAnalysised$SurvivalRate_data_list[[SurvivalRate_data_name]] %>%
    anti_join(SI_for_biorep_ori_outlier_output_list[[SI_for_biorep_ori_outlier_name]]) %>%
    anti_join(SI_for_biorep_curve_outlier_output_list[[SI_for_biorep_curve_outlier_name]])
  return(output)
},mc.cores = 30)

#11.Calculate the clean SI for biorep data , following remove the strain which just has only one biorep#####
SI_for_biorep_clean_SI_biorep_ori_and_curve_Argument_data <- SurvivalRate_clean_SI_biorep_ori_and_curve_Argument_data %>%
  tibble::rownames_to_column() %>%
  dplyr::rename(SurvivalRate_clean_SI_biorep_ori_and_curve_data_name=rowname) %>%
  mutate(SI_for_biorep_clean_SI_biorep_ori_and_curve_name = paste("SI_for_biorep_clean_SI_biorep_ori_and_curve",
                                                                  gsub("SI_for_biorep_ori_(.*)","\\1",
                                                                       SI_for_biorep_ori_name),sep = "_") ) %>%
  tibble::column_to_rownames(var = "SI_for_biorep_clean_SI_biorep_ori_and_curve_name")

SI_for_biorep_clean_SI_biorep_ori_and_curve_output_list <- mclapply(
  split(SI_for_biorep_clean_SI_biorep_ori_and_curve_Argument_data,f = rownames(
    SI_for_biorep_clean_SI_biorep_ori_and_curve_Argument_data)),FUN = function(argument_data){
  SurvivalRate_clean_SI_biorep_ori_and_curve_data_name <- 
    argument_data$SurvivalRate_clean_SI_biorep_ori_and_curve_data_name
  output <- SurvivalIntegral_calculate(
    SurvivalRateData = 
      SurvivalRate_clean_SI_biorep_ori_and_curve_output_list[[SurvivalRate_clean_SI_biorep_ori_and_curve_data_name]],
    default_rm_vector = NA,Relative_day = "day2",Categorical_variable = Categorical_variable,retain_biorep = T,
    Remove_not_enough_datapoint = T,datapoint_number = 6)
  
  return(output)
},mc.cores = 30)

#12.Aggregate SI data to strain , remove the strain which just has only one biorep (sd == NA)######
SI_for_strain_clean_SI_biorep_ori_and_curve_rm_na_sd_Argument_data <- 
  SI_for_biorep_clean_SI_biorep_ori_and_curve_Argument_data %>%
  tibble::rownames_to_column() %>%
  dplyr::rename(SI_for_biorep_clean_SI_biorep_ori_and_curve_name=rowname) %>%
  mutate(SI_for_strain_clean_SI_biorep_ori_and_curve_rm_na_sd_name=paste(
    "SI_for_strain_clean_SI_biorep_ori_and_curve_rm_na_sd",
    gsub("SI_for_biorep_ori_(.*)","\\1",SI_for_biorep_ori_name),sep = "_") ) %>%
  tibble::column_to_rownames(var = "SI_for_strain_clean_SI_biorep_ori_and_curve_rm_na_sd_name")

SI_for_strain_clean_SI_biorep_ori_and_curve_rm_na_sd_output_list <- 
  mclapply(split(SI_for_strain_clean_SI_biorep_ori_and_curve_rm_na_sd_Argument_data,
                 f = rownames(SI_for_strain_clean_SI_biorep_ori_and_curve_rm_na_sd_Argument_data)),
           FUN = function(argument_data){
  SI_for_biorep_clean_SI_biorep_ori_and_curve_name <- argument_data$SI_for_biorep_clean_SI_biorep_ori_and_curve_name
  output <- SI_for_biorep_clean_SI_biorep_ori_and_curve_output_list[[SI_for_biorep_clean_SI_biorep_ori_and_curve_name]] %>%
    group_by(strains,strains_group,tube_type_YPD,tube_type_SC,glu_concentration) %>%
    dplyr::summarise(sd_SI=sd(SI),SI_mean_bio = mean(SI)) %>%
    filter(!is.na(sd_SI))
  return(output)
},mc.cores = 30)

# MAIN#######
#####label the BY4741 panel manually#####
####12 panel(strains group) in total
# 1.0126BR2A1A8 (not contain in Res data)
# 2.4A1A10 
# 3.4A11B8(BY4741 not grow in SC,no doubling time)
# 4.4B9C6(BY4741 not grow in SC,no doubling time)
# 5.11A1A10
# 6.11A11B8
# 7.11B9C6
# 8.11C7D4
# 9.BYRM11D5D12 (not contain in Res data , because day4 data missing)
# 10.11E1E10
# 11.11E11F8
# 12.11F9G6
BY4741_strains_group <- c("0126BR2A1A8","4A1A10","4A11B8","4B9C6","11A1A10",
                          "11A11B8","11B9C6","11C7D4","BYRM11D5D12","11E1E10",
                          "11E11F8","11F9G6")

#0.Argument_data####
Argument_data <- SI_biorep_clean_list$Argument_data

#########
#IN THIS PART , list name will not change#####
#########

#1.get WT(BY/RM) data and rm WT data######
#BY######
BY_for_biorep_data_output_list <- mclapply(split(Argument_data,
                                                 f = Argument_data$SI_for_biorep_clean_SI_biorep_ori_and_curve_name),
                                           FUN = function(argument_data){
  SI_for_biorep_clean_SI_biorep_ori_and_curve_name=argument_data$SI_for_biorep_clean_SI_biorep_ori_and_curve_name
  output <- Filter_wild_type_strain_and_rename(
    SI_for_strain_or_biorep_data = SI_biorep_clean_list$SI_for_biorep_clean_SI_biorep_ori_and_curve_output_list[[SI_for_biorep_clean_SI_biorep_ori_and_curve_name]],
    wild_type_name = "B",BY4741_strains_group = BY4741_strains_group) 
  return(output)
},mc.cores = 30)
#for strain data , BY strain data just have the biorep strains_group data
BY_for_strain_data_output_list <-  mclapply(split(Argument_data,f = rownames(Argument_data)),FUN = function(argument_data){
  SI_for_strain_clean_SI_biorep_ori_and_curve_name=rownames(argument_data)
  output <- Filter_wild_type_strain_and_rename(
    SI_for_strain_or_biorep_data = SI_biorep_clean_list$SI_for_strain_clean_SI_biorep_ori_and_curve_rm_na_sd_output_list[[SI_for_strain_clean_SI_biorep_ori_and_curve_name]],
    wild_type_name = "B",BY4741_strains_group = BY4741_strains_group) 
  return(output)
},mc.cores = 30)

#RM########
RM_for_biorep_data_output_list <- mclapply(split(Argument_data,f = Argument_data$SI_for_biorep_clean_SI_biorep_ori_and_curve_name),FUN = function(argument_data){
  SI_for_biorep_clean_SI_biorep_ori_and_curve_name=argument_data$SI_for_biorep_clean_SI_biorep_ori_and_curve_name
  output <- Filter_wild_type_strain_and_rename(
    SI_for_strain_or_biorep_data = SI_biorep_clean_list$SI_for_biorep_clean_SI_biorep_ori_and_curve_output_list[[SI_for_biorep_clean_SI_biorep_ori_and_curve_name]],
    wild_type_name = "R",BY4741_strains_group = BY4741_strains_group) 
  return(output)
},mc.cores = 30)
RM_for_strain_data_output_list <- mclapply(split(Argument_data,f = rownames(Argument_data)),FUN = function(argument_data){
  SI_for_strain_clean_SI_biorep_ori_and_curve_name=rownames(argument_data)
  output <- Filter_wild_type_strain_and_rename(
    SI_for_strain_or_biorep_data = SI_biorep_clean_list$SI_for_strain_clean_SI_biorep_ori_and_curve_rm_na_sd_output_list[[SI_for_strain_clean_SI_biorep_ori_and_curve_name]],
    wild_type_name = "R",BY4741_strains_group = BY4741_strains_group) 
  return(output)
},mc.cores = 30)

#remove WT data#######
SI_for_biorep_clean_SI_biorep_ori_and_curve_rm_wt_output_list <- mclapply(
  split(Argument_data,f = Argument_data$SI_for_biorep_clean_SI_biorep_ori_and_curve_name),FUN = function(argument_data){
  SI_for_biorep_clean_SI_biorep_ori_and_curve_name=argument_data$SI_for_biorep_clean_SI_biorep_ori_and_curve_name
  output <- SI_wild_type_remove(
    SI_for_strain_or_biorep_data = SI_biorep_clean_list$SI_for_biorep_clean_SI_biorep_ori_and_curve_output_list[[SI_for_biorep_clean_SI_biorep_ori_and_curve_name]],
    String = "^B|R")
  return(output)
},mc.cores = 30)

SI_for_strain_clean_SI_biorep_ori_and_curve_rm_wt_output_list <- mclapply(
  split(Argument_data,f = rownames(Argument_data)),FUN = function(argument_data){
  SI_for_strain_clean_SI_biorep_ori_and_curve_name=rownames(argument_data)
  output <- SI_wild_type_remove(
    SI_for_strain_or_biorep_data = SI_biorep_clean_list$SI_for_strain_clean_SI_biorep_ori_and_curve_rm_na_sd_output_list[[SI_for_strain_clean_SI_biorep_ori_and_curve_name]],
    String = "^B|R")
  return(output)
},mc.cores = 30)

#2.get glu_concentration==0.5 data#######
Glu_concentration_for_biorep_data_output_list <- mclapply(
  split(Argument_data,f = Argument_data$SI_for_biorep_clean_SI_biorep_ori_and_curve_name),FUN = function(argument_data){
  SI_for_biorep_clean_SI_biorep_ori_and_curve_name=argument_data$SI_for_biorep_clean_SI_biorep_ori_and_curve_name
  output <- Filter_glu_concentration(
    SI_for_strain_or_biorep_data = SI_for_biorep_clean_SI_biorep_ori_and_curve_rm_wt_output_list[[SI_for_biorep_clean_SI_biorep_ori_and_curve_name]],
    is.filter.out = F)
  return(output)
},mc.cores = 30)

# Glu_concentration_for_strain_data_output_list is NULL because Glu_concentration_for_strain biorep is NA, 
#filter by !is.na(sd(SI))
Glu_concentration_for_strain_data_output_list <- mclapply(split(Argument_data,f = rownames(Argument_data)),
                                                          FUN = function(argument_data){
  SI_for_strain_clean_SI_biorep_ori_and_curve_name=rownames(argument_data)
  output <- Filter_glu_concentration(
    SI_for_strain_or_biorep_data = SI_for_strain_clean_SI_biorep_ori_and_curve_rm_wt_output_list[[SI_for_strain_clean_SI_biorep_ori_and_curve_name]],
    is.filter.out = F)
  return(output)
},mc.cores = 30) 

#remove glu_concentration data
SI_for_biorep_clean_SI_biorep_ori_and_curve_rm_wt_gc_output_list <- mclapply(
  split(Argument_data,f = Argument_data$SI_for_biorep_clean_SI_biorep_ori_and_curve_name),FUN = function(argument_data){
  SI_for_biorep_clean_SI_biorep_ori_and_curve_name=argument_data$SI_for_biorep_clean_SI_biorep_ori_and_curve_name
  output <- Filter_glu_concentration(
    SI_for_strain_or_biorep_data = SI_for_biorep_clean_SI_biorep_ori_and_curve_rm_wt_output_list[[SI_for_biorep_clean_SI_biorep_ori_and_curve_name]],
    is.filter.out = T)
  return(output)
},mc.cores = 30)
SI_for_strain_clean_SI_biorep_ori_and_curve_rm_wt_gc_output_list <- mclapply(
  split(Argument_data,f = rownames(Argument_data)),FUN = function(argument_data){
  SI_for_strain_clean_SI_biorep_ori_and_curve_name=rownames(argument_data)
  output <- Filter_glu_concentration(
    SI_for_strain_or_biorep_data = SI_for_strain_clean_SI_biorep_ori_and_curve_rm_wt_output_list[[SI_for_strain_clean_SI_biorep_ori_and_curve_name]],
    is.filter.out = T)
  return(output)
},mc.cores = 30)

#3.remove duplicate (rm tube_type_SC != "C" OR "J" ) , contain the same strain but different strains_group data (in biorep format) ,#######
#  and remove SI theoretical outlier(>15) , and summarize to SI_for_strain_rm_wt_gc_dp######
#3.1 get duplicate data (SI value with tube_type_SC == "G","Y" OR tube_type_SC == "J" when contain strain with tube_type_SC == "C" data)#####
#retrun SI_for_biorep data 
Tube_type_duplicate_for_biorep_data_output_list <- mclapply(
  split(Argument_data,f = Argument_data$SI_for_biorep_clean_SI_biorep_ori_and_curve_name),
  FUN = function(argument_data){
  SI_for_biorep_clean_SI_biorep_ori_and_curve_name=argument_data$SI_for_biorep_clean_SI_biorep_ori_and_curve_name
  output <- remove_duplicate_reconstruct(
    SI_for_biorep_or_strain_data = SI_for_biorep_clean_SI_biorep_ori_and_curve_rm_wt_gc_output_list[[SI_for_biorep_clean_SI_biorep_ori_and_curve_name]],
    is.biorep.data = T,is.filter.out = F,Group_var = "strains_group")
  return(output)
},mc.cores = 30)

Tube_type_duplicate_for_strain_data_output_list <- mclapply(
  split(Argument_data,f = rownames(Argument_data)),FUN = function(argument_data){
  SI_for_strain_clean_SI_biorep_ori_and_curve_name=rownames(argument_data)
  output <- remove_duplicate_reconstruct(
    SI_for_biorep_or_strain_data = SI_for_strain_clean_SI_biorep_ori_and_curve_rm_wt_gc_output_list[[SI_for_strain_clean_SI_biorep_ori_and_curve_name]],
    is.biorep.data = F,is.filter.out = F,Group_var = "strains_group")
  return(output)
},mc.cores = 30)

#3.2 get remove duplicate tube_type SI_for_biorep_rm_wt_gc_dp_data###########
SI_for_biorep_clean_SI_biorep_ori_and_curve_rm_wt_gc_dp_data_output_list <- mclapply(
  split(Argument_data,f = Argument_data$SI_for_biorep_clean_SI_biorep_ori_and_curve_name),
  FUN = function(argument_data){
  SI_for_biorep_clean_SI_biorep_ori_and_curve_name=argument_data$SI_for_biorep_clean_SI_biorep_ori_and_curve_name
  output <- remove_duplicate_reconstruct(
    SI_for_biorep_or_strain_data = SI_for_biorep_clean_SI_biorep_ori_and_curve_rm_wt_gc_output_list[[SI_for_biorep_clean_SI_biorep_ori_and_curve_name]],
    is.biorep.data = T,is.filter.out = T,Group_var = "strains_group")
  return(output)
},mc.cores = 30)

SI_for_strain_clean_SI_biorep_ori_and_curve_rm_wt_gc_dp_data_output_list <- mclapply(
  split(Argument_data,f = rownames(Argument_data)),FUN = function(argument_data){
  SI_for_strain_clean_SI_biorep_ori_and_curve_name=rownames(argument_data)
  output <- remove_duplicate_reconstruct(
    SI_for_biorep_or_strain_data = SI_for_strain_clean_SI_biorep_ori_and_curve_rm_wt_gc_output_list[[SI_for_strain_clean_SI_biorep_ori_and_curve_name]],
    is.biorep.data = F,is.filter.out = T,Group_var = "strains_group")
  return(output)
},mc.cores = 30)

#3.3 remove duplicate strain in different strains_group##############
# also use the function remove_duplicate_reconstruct
Strains_group_duplicate_for_biorep_data_list_output <-  mclapply(
  split(Argument_data,f = Argument_data$SI_for_biorep_clean_SI_biorep_ori_and_curve_name),
  FUN = function(argument_data){
  SI_for_biorep_clean_SI_biorep_ori_and_curve_name=argument_data$SI_for_biorep_clean_SI_biorep_ori_and_curve_name
  output <- remove_duplicate_reconstruct(
    SI_for_biorep_or_strain_data = SI_for_biorep_clean_SI_biorep_ori_and_curve_rm_wt_gc_dp_data_output_list[[SI_for_biorep_clean_SI_biorep_ori_and_curve_name]],
    is.biorep.data = T,is.filter.out = F,Group_var = NULL)
  return(output)
},mc.cores = 30)

Strains_group_duplicate_for_strain_data_list_output <- mclapply(
  split(Argument_data,f = rownames(Argument_data)),FUN = function(argument_data){
  SI_for_strain_clean_SI_biorep_ori_and_curve_name=rownames(argument_data)
  output <- remove_duplicate_reconstruct(
    SI_for_biorep_or_strain_data = SI_for_strain_clean_SI_biorep_ori_and_curve_rm_wt_gc_dp_data_output_list[[SI_for_strain_clean_SI_biorep_ori_and_curve_name]],
    is.biorep.data = F,is.filter.out = F,Group_var = NULL)
  return(output)
},mc.cores = 30)

SI_for_biorep_clean_SI_biorep_ori_and_curve_rm_wt_gc_dp_dpsg_data_output_list <- mclapply(
  split(Argument_data,f = Argument_data$SI_for_biorep_clean_SI_biorep_ori_and_curve_name),
  FUN = function(argument_data){
  SI_for_biorep_clean_SI_biorep_ori_and_curve_name=argument_data$SI_for_biorep_clean_SI_biorep_ori_and_curve_name
  output <- remove_duplicate_reconstruct(
    SI_for_biorep_or_strain_data = SI_for_biorep_clean_SI_biorep_ori_and_curve_rm_wt_gc_dp_data_output_list[[SI_for_biorep_clean_SI_biorep_ori_and_curve_name]],
    is.biorep.data = T,is.filter.out = T,Group_var = NULL)
  return(output)
},mc.cores = 30)

SI_for_strain_clean_SI_biorep_ori_and_curve_rm_wt_gc_dp_dpsg_data_output_list <- mclapply(
  split(Argument_data,f = rownames(Argument_data)),FUN = function(argument_data){
  SI_for_strain_clean_SI_biorep_ori_and_curve_name=rownames(argument_data)
  output <- remove_duplicate_reconstruct(
    SI_for_biorep_or_strain_data = SI_for_strain_clean_SI_biorep_ori_and_curve_rm_wt_gc_dp_data_output_list[[SI_for_strain_clean_SI_biorep_ori_and_curve_name]],
    is.biorep.data = F,is.filter.out = T,Group_var = NULL)
  return(output)
},mc.cores = 30)

#3.4 remove wrong marker data ################
SI_for_biorep_clean_SI_biorep_ori_and_curve_rm_wt_gc_dp_dpsg_wm_data_output_list <- mclapply(
  split(Argument_data,f = Argument_data$SI_for_biorep_clean_SI_biorep_ori_and_curve_name),
  FUN = function(argument_data){
  SI_for_biorep_clean_SI_biorep_ori_and_curve_name=argument_data$SI_for_biorep_clean_SI_biorep_ori_and_curve_name
  output <- SI_for_biorep_clean_SI_biorep_ori_and_curve_rm_wt_gc_dp_dpsg_data_output_list[[SI_for_biorep_clean_SI_biorep_ori_and_curve_name]] %>%
    filter(!strains %in% ErrorStrains_vector)
  return(output)
},mc.cores = 30)

SI_for_strain_clean_SI_biorep_ori_and_curve_rm_wt_gc_dp_dpsg_wm_data_output_list <- mclapply(
  split(Argument_data,f = rownames(Argument_data)),FUN = function(argument_data){
  SI_for_strain_clean_SI_biorep_ori_and_curve_name=rownames(argument_data)
  output <-SI_for_strain_clean_SI_biorep_ori_and_curve_rm_wt_gc_dp_dpsg_data_output_list[[SI_for_strain_clean_SI_biorep_ori_and_curve_name]] %>%
    filter(!strains %in% ErrorStrains_vector)
  return(output)
},mc.cores = 30)

#3.5 remove theoretical outlier >15 in biorep data (A0572)#############
Theoretial_outlier_for_biorep_data <- mclapply(
  split(Argument_data,f = Argument_data$SI_for_biorep_clean_SI_biorep_ori_and_curve_name),
  FUN = function(argument_data){
  SI_for_biorep_clean_SI_biorep_ori_and_curve_name=argument_data$SI_for_biorep_clean_SI_biorep_ori_and_curve_name
  output <- SI_for_biorep_clean_SI_biorep_ori_and_curve_rm_wt_gc_dp_dpsg_wm_data_output_list[[SI_for_biorep_clean_SI_biorep_ori_and_curve_name]] %>%
    filter(strains %in% c("A0572"))
  return(output)
},mc.cores = 30)

Theoretial_outlier_for_strain_data <-  mclapply(
  split(Argument_data,f = rownames(Argument_data)),FUN = function(argument_data){
  SI_for_strain_clean_SI_biorep_ori_and_curve_name=rownames(argument_data)
  output <-SI_for_strain_clean_SI_biorep_ori_and_curve_rm_wt_gc_dp_dpsg_wm_data_output_list[[SI_for_strain_clean_SI_biorep_ori_and_curve_name]] %>%
    filter(strains %in%  c("A0572"))
  return(output)
},mc.cores = 30)

SI_for_biorep_clean_SI_biorep_ori_and_curve_rm_wt_gc_dp_dpsg_wm_ol15_data_output_list <-  mclapply(
  split(Argument_data,f = Argument_data$SI_for_biorep_clean_SI_biorep_ori_and_curve_name),
  FUN = function(argument_data){
  SI_for_biorep_clean_SI_biorep_ori_and_curve_name=argument_data$SI_for_biorep_clean_SI_biorep_ori_and_curve_name
  output <- SI_for_biorep_clean_SI_biorep_ori_and_curve_rm_wt_gc_dp_dpsg_wm_data_output_list[[SI_for_biorep_clean_SI_biorep_ori_and_curve_name]] %>%
    filter(!strains %in% c("A0572"))
  return(output)
},mc.cores = 30)

SI_for_strain_clean_SI_biorep_ori_and_curve_rm_wt_gc_dp_dpsg_wm_ol15_data_output_list <-  mclapply(
  split(Argument_data,f = rownames(Argument_data)),FUN = function(argument_data){
  SI_for_strain_clean_SI_biorep_ori_and_curve_name=rownames(argument_data)
  output <-SI_for_strain_clean_SI_biorep_ori_and_curve_rm_wt_gc_dp_dpsg_wm_data_output_list[[SI_for_strain_clean_SI_biorep_ori_and_curve_name]] %>%
    filter(!strains %in%  c("A0572"))
  return(output)
},mc.cores = 30)

#SAVE OBJECT#####
SI_rm_wt_gc_dp_dpsg_wm_ol15_list <- list(Marker_wrong_strains=ErrorStrains_vector,
                                         BY_for_biorep_data_output_list=BY_for_biorep_data_output_list,
                                         BY_for_strain_data_output_list=BY_for_strain_data_output_list,
                                         RM_for_biorep_data_output_list=RM_for_biorep_data_output_list,
                                         RM_for_strain_data_output_list=RM_for_strain_data_output_list,
                                         SI_for_biorep_clean_SI_biorep_ori_and_curve_rm_wt_output_list=
                                           SI_for_biorep_clean_SI_biorep_ori_and_curve_rm_wt_output_list,
                                         SI_for_strain_clean_SI_biorep_ori_and_curve_rm_wt_output_list=
                                           SI_for_strain_clean_SI_biorep_ori_and_curve_rm_wt_output_list,
                                         Glu_concentration_for_biorep_data_output_list=
                                           Glu_concentration_for_biorep_data_output_list,
                                         Glu_concentration_for_strain_data_output_list=
                                           Glu_concentration_for_strain_data_output_list,
                                         SI_for_biorep_clean_SI_biorep_ori_and_curve_rm_wt_gc_output_list=
                                           SI_for_biorep_clean_SI_biorep_ori_and_curve_rm_wt_gc_output_list,
                                         SI_for_strain_clean_SI_biorep_ori_and_curve_rm_wt_gc_output_list=
                                           SI_for_strain_clean_SI_biorep_ori_and_curve_rm_wt_gc_output_list,
                                         Tube_type_duplicate_for_biorep_data_output_list=
                                           Tube_type_duplicate_for_biorep_data_output_list,
                                         Tube_type_duplicate_for_strain_data_output_list=
                                           Tube_type_duplicate_for_strain_data_output_list,
                                         SI_for_biorep_clean_SI_biorep_ori_and_curve_rm_wt_gc_dp_data_output_list=
                                           SI_for_biorep_clean_SI_biorep_ori_and_curve_rm_wt_gc_dp_data_output_list,
                                         SI_for_strain_clean_SI_biorep_ori_and_curve_rm_wt_gc_dp_data_output_list=
                                           SI_for_strain_clean_SI_biorep_ori_and_curve_rm_wt_gc_dp_data_output_list,
                                         Strains_group_duplicate_for_biorep_data_list_output=
                                           Strains_group_duplicate_for_biorep_data_list_output,
                                         Strains_group_duplicate_for_strain_data_list_output=
                                           Strains_group_duplicate_for_strain_data_list_output,
                                         SI_for_biorep_clean_SI_biorep_ori_and_curve_rm_wt_gc_dp_dpsg_data_output_list=
                                           SI_for_biorep_clean_SI_biorep_ori_and_curve_rm_wt_gc_dp_dpsg_data_output_list,
                                         SI_for_strain_clean_SI_biorep_ori_and_curve_rm_wt_gc_dp_dpsg_data_output_list=
                                           SI_for_strain_clean_SI_biorep_ori_and_curve_rm_wt_gc_dp_dpsg_data_output_list,
                                         Theoretial_outlier_for_biorep_data=Theoretial_outlier_for_biorep_data,
                                         Theoretial_outlier_for_strain_data=Theoretial_outlier_for_strain_data,
                                         SI_for_biorep_clean_SI_biorep_ori_and_curve_rm_wt_gc_dp_dpsg_wm_ol15_data_output_list=
                                           SI_for_biorep_clean_SI_biorep_ori_and_curve_rm_wt_gc_dp_dpsg_wm_ol15_data_output_list,
                                         SI_for_strain_clean_SI_biorep_ori_and_curve_rm_wt_gc_dp_dpsg_wm_ol15_data_output_list=
                                           SI_for_strain_clean_SI_biorep_ori_and_curve_rm_wt_gc_dp_dpsg_wm_ol15_data_output_list)

# save
set_dir <- "./"
save(SI_rm_wt_gc_dp_dpsg_wm_ol15_list,file = paste(set_dir,"3.SI_data",Sys.Date(),".Rdata",sep = ""))

