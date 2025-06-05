###quality control of translation error
#library
library(reshape2)
library(stringr)
library(dplyr)

# source function
source("./5_1.split_TE_data_and_F_R_value_control_function.R")
source("./5_2.TE_techrep_quality_control_function.R")
source("./5_3.TE_11_F_R_ratio_calculate_and_biorep_control_function.R")

#####load data
#this data contains machine fixed data,the first part data,and machine fixed supplment data
#and this Res data had reshape format to machine fixed data type , rename BY strain , and removed duplicate strain
load("./5.raw_TE_data.Rdata")

#main
#get outlier BY panel

######Quality control#########
##use normal method control data 
##(input TE_Res_data_all to run : 1.BY filter , 2.11 13 F R value filter , 3.11 biorep filter & 13 biorep filter)
#Filter NO.1 : find the not good panel which BY TE is uncorrected
BY_data <- TE_Res_data_all %>%
  filter(strain == "BY4716")

#outlier panel
outlierPanel <- BY_data %>% 
  mutate(F_R_ratio = `F`/`R`) %>%
  group_by(strain,mutant_type,biorep,panel) %>%
  dplyr::summarize(F_R_ratio = mean(F_R_ratio)) %>%
  dcast(strain + panel ~ mutant_type , value.var = "F_R_ratio") %>%
  mutate(TE = `13`/`11`) %>%
  filter(TE >= 1e-2) %>%
  pull(panel)

TE_Res_data_all_rm_panel_outlier <-  TE_Res_data_all %>% filter(!panel %in% outlierPanel)
TE_Res_data_all_outlier_panel_data <- TE_Res_data_all %>% filter(panel %in% outlierPanel)

# MAIN#######
#Filter NO.2: split strain data to 11 strains data & 13 strains data , do F R value filter respectively
TE_Res_data_split_list <- TE_Res_data_all_rm_panel_outlier %>% filter(!panel %in% outlierPanel) %>%
  group_by(mutant_type) %>%
  group_split()

TE_11_F_R_outlier <- TE_F_R_value_filter(data = TE_Res_data_split_list[[1]],mutant_type = 11,is.unique = F,
                                         return.outlier = T,F_11 = 1e+06,R_11_illumiate = 1e+06,R_11 = 1e+08)
TE_13_F_R_outlier <- TE_F_R_value_filter(data = TE_Res_data_split_list[[2]],mutant_type = 13,is.unique = F,
                                         return.outlier = T,F_13 = 500,R_13_illumiate = 1e+06,R_13_threshold = 1e+08)

TE_11_clean_F_R_data <- TE_F_R_value_filter(data = TE_Res_data_split_list[[1]],mutant_type = 11,is.unique = F,
                                            F_11 = 1e+06,R_11_illumiate = 1e+06,R_11 = 1e+08)
TE_13_clean_F_R_data <- TE_F_R_value_filter(data = TE_Res_data_split_list[[2]],mutant_type = 13,is.unique = F,
                                            return.outlier = F,F_13 = 500,R_13_illumiate = 1e+06,R_13_threshold = 1e+08)

#rm the biorep which just have only one techrep in 11/13
TE_11_clean_F_R_data_rm_o1_techrep <- remove_only_one_techrep(data = TE_11_clean_F_R_data,return.one.techrep = F)
TE_13_clean_F_R_data_rm_o1_techrep <- remove_only_one_techrep(data = TE_13_clean_F_R_data,return.one.techrep = F)

#get 11/13 only one techrep data 
TE_11_one_techrep <- remove_only_one_techrep(data = TE_11_clean_F_R_data,return.one.techrep = T)
TE_13_one_techrep <- remove_only_one_techrep(data = TE_13_clean_F_R_data,return.one.techrep = T)

#main
######11 strain techrep outlier#############
####find techrep outlier, first time####
F_R_11_outlier_1_data <- find_TE_techrep_outlier(TE_techrep_data = TE_11_clean_F_R_data_rm_o1_techrep,
                                                 default_rm_vector = NA,CV_cutoff = 0.2,Categorical_variable = NULL)

TE_11_clean_techrep_data_1 <- TE_techrep_clean(data = TE_11_clean_F_R_data_rm_o1_techrep,TEoutlier_data = F_R_11_outlier_1_data)

####find techrep outlier, second time####
F_R_11_outlier_2_data <- find_TE_techrep_outlier(TE_techrep_data = TE_11_clean_techrep_data_1,
                                                 default_rm_vector =NA,CV_cutoff = 0.2,Categorical_variable=NULL)

#bind_rows ODoutlier data to get techrep outlier in first time and second time #####
F_R_11_outlier_data <- bind_rows(F_R_11_outlier_1_data,F_R_11_outlier_2_data)

##clean techrep outlier before mean 11 techrep data
TE_11_techrep_clean_data <- TE_techrep_clean(data = TE_11_clean_F_R_data_rm_o1_techrep,TEoutlier_data = F_R_11_outlier_data)

#######13 strain techrep outlier############
F_R_13_outlier_1_data <- find_TE_techrep_outlier(TE_techrep_data = TE_13_clean_F_R_data_rm_o1_techrep,
                                                 default_rm_vector = NA,CV_cutoff = 0.2,Categorical_variable = NULL)
TE_13_clean_techrep_data_1 <- TE_techrep_clean(data = TE_13_clean_F_R_data_rm_o1_techrep,TEoutlier_data = F_R_13_outlier_1_data)

####find techrep outlier, second time####
F_R_13_outlier_2_data <- find_TE_techrep_outlier(TE_techrep_data = TE_13_clean_techrep_data_1,
                                                 default_rm_vector =NA,CV_cutoff = 0.2,Categorical_variable=NULL)

#bind_rows ODoutlier data to get techrep outlier in first time and second time #####
F_R_13_outlier_data <- bind_rows(F_R_13_outlier_1_data,F_R_13_outlier_2_data)

##clean techrep outlier before mean 13 techrep data
TE_13_techrep_clean_data <- TE_techrep_clean(data = TE_13_clean_F_R_data_rm_o1_techrep,TEoutlier_data = F_R_13_outlier_data)

#save
TE_11_techrep_clean_list <- list(F_R_11_outlier_data=F_R_11_outlier_data,TE_11_techrep_clean_data=TE_11_techrep_clean_data)
TE_13_techrep_clean_list <- list(F_R_13_outlier_data=F_R_13_outlier_data,TE_13_techrep_clean_data=TE_13_techrep_clean_data)
TE_techrep_clean_list <- list(TE_11_techrep_clean_list=TE_11_techrep_clean_list,
                              TE_13_techrep_clean_list=TE_13_techrep_clean_list)

set_dir <- "./"
save(TE_techrep_clean_list,file = paste(set_dir,"5.TE_techrep_data_",Sys.Date(),".Rdata",sep = ""))

#main
#####before TE 11 data biorep clean , should first calculate F_R_ratio value for each techrep , 
#and summarize TE_11_techrep_data to TE_11_biorep_data ###########
TE_11_biorep_data <- TE_techrep_clean_list$TE_11_techrep_clean_list$TE_11_techrep_clean_data %>%
  mutate(F_R_ratio = `F`/`R`) %>%
  group_by(strain,mutant_type,biorep,panel) %>%
  dplyr::summarise(sd_F_R_ratio = sd(F_R_ratio) , F_R_ratio = mean(F_R_ratio))

#####TE 11 biorep data clean part ##########
TE_11_biorep_outlier_1_data <- find_TE_11_biorep_outlier(TE_11_Data_contain_bio = TE_11_biorep_data,cutoff_type = "CV",
                                                         CV_cutoff = 0.2,Categorical_variable = NULL)

#anti_join to get first biorep clean
TE_11_biorep_clean_1 <- TE_11_biorep_data %>% anti_join(TE_11_biorep_outlier_1_data,
                                                        by=c("strain","mutant_type","biorep","panel"))

#find biorep outlier second time
TE_11_biorep_outlier_2 <- find_TE_11_biorep_outlier(TE_11_Data_contain_bio = TE_11_biorep_clean_1,cutoff_type = "CV",
                                                    CV_cutoff = 0.2,Categorical_variable = NULL)

#get outlier in the first and second time
TE_11_biorep_outlier_data <- bind_rows(TE_11_biorep_outlier_1_data,TE_11_biorep_outlier_2)

#anti_join outlier to get TE 11 biorep clean data for save and summary TE 11 biorep data to TE 11 strain data
TE_11_for_biorep_clean_bio_data <- TE_11_biorep_data  %>% anti_join(TE_11_biorep_outlier_data,
                                                                    by=c("strain","mutant_type","biorep","panel"))

# #summary to strain data
TE_11_for_strain_clean_bio_data <- TE_11_for_biorep_clean_bio_data %>%
  group_by_at(c("strain","mutant_type","biorep")) %>%
  dplyr::summarize(sd_F_R_ratio = sd(F_R_ratio),F_R_ratio = mean(F_R_ratio)) %>%
  group_by_at(c("strain","mutant_type")) %>%
  dplyr::summarize(sd_F_R_ratio = sd(F_R_ratio),F_R_ratio_mean_bio = mean(F_R_ratio))

TE_11_for_strain_clean_bio_have_biorep_data <- TE_11_for_strain_clean_bio_data %>%
  filter(!is.na(sd_F_R_ratio))

TE_11_for_strain_clean_bio_only_one_biorep_data <-  TE_11_for_strain_clean_bio_data %>%
  filter(is.na(sd_F_R_ratio))

#save
TE_11_biorep_clean_list <- list(TE_11_biorep_data=TE_11_biorep_data,
                                TE_11_biorep_outlier_data=TE_11_biorep_outlier_data,
                                TE_11_for_biorep_clean_bio_data=TE_11_for_biorep_clean_bio_data,
                                TE_11_for_strain_clean_bio_data=TE_11_for_strain_clean_bio_data,
                                TE_11_for_strain_clean_bio_have_biorep_data=TE_11_for_strain_clean_bio_have_biorep_data,
                                TE_11_for_strain_clean_bio_only_one_biorep_data=TE_11_for_strain_clean_bio_only_one_biorep_data)

set_dir <- "./"
save(TE_11_biorep_clean_list,file = paste(set_dir,"5.TE_11_biorep_data_",Sys.Date(),".Rdata",sep = ""))

#main
#####before TE 13 data biorep clean , should first calculate F_R_ratio value for each techrep , 
#and summarize TE_13_techrep_data to TE_13_biorep_data ###########
TE_13_biorep_data <- TE_techrep_clean_list$TE_13_techrep_clean_list$TE_13_techrep_clean_data %>%
  mutate(F_R_ratio = `F`/`R`) %>%
  group_by(strain,mutant_type,biorep,panel) %>%
  dplyr::summarise(sd_F_R_ratio = sd(F_R_ratio) , F_R_ratio = mean(F_R_ratio))

#####TE 13 biorep data clean part ##########
TE_13_biorep_outlier_1_data <- find_TE_11_biorep_outlier(TE_11_Data_contain_bio = TE_13_biorep_data,cutoff_type = "CV",
                                                         CV_cutoff = 0.2,Categorical_variable = NULL)

#anti_join to get first biorep clean
TE_13_biorep_clean_1 <- TE_13_biorep_data %>% anti_join(TE_13_biorep_outlier_1_data,
                                                        by=c("strain","mutant_type","biorep","panel"))

#find biorep outlier second time
TE_13_biorep_outlier_2 <- find_TE_11_biorep_outlier(TE_11_Data_contain_bio = TE_13_biorep_clean_1,cutoff_type = "CV",
                                                    CV_cutoff = 0.2,Categorical_variable = NULL)

#get outlier in the first and second time
TE_13_biorep_outlier_data <- bind_rows(TE_13_biorep_outlier_1_data,TE_13_biorep_outlier_2)

#anti_join outlier to get TE 13 biorep clean data for save and summary TE 13 biorep data to TE 13 strain data
TE_13_for_biorep_clean_bio_data <- TE_13_biorep_data  %>% anti_join(TE_13_biorep_outlier_data,
                                                                    by=c("strain","mutant_type","biorep","panel"))

# #summary to strain data
TE_13_for_strain_clean_bio_data <- TE_13_for_biorep_clean_bio_data %>%
  group_by_at(c("strain","mutant_type")) %>%
  dplyr::summarize(sd_F_R_ratio = sd(F_R_ratio),F_R_ratio_mean_bio = mean(F_R_ratio))

#save
TE_13_biorep_clean_list <- list(TE_13_biorep_data=TE_13_biorep_data,
                                TE_13_biorep_outlier_data=TE_13_biorep_outlier_data,
                                TE_13_for_biorep_clean_bio_data=TE_13_for_biorep_clean_bio_data,
                                TE_13_for_strain_clean_bio_data=TE_13_for_strain_clean_bio_data)

set_dir <- "./"
save(TE_13_biorep_clean_list,file = paste(set_dir,"5.TE_13_biorep_data_",Sys.Date(),".RData",sep = ""))

##main
#####TE_11_for_strain_clean_bio_data object ######
TE_11_for_strain_clean_bio_have_biorep_data <- TE_11_biorep_clean_list$TE_11_for_strain_clean_bio_have_biorep_data

#####TE_13_for_biorep_clean_bio_data object ######
TE_13_for_biorep_clean_bio_data <- TE_13_biorep_clean_list$TE_13_for_biorep_clean_bio_data

#calculate TE#########
#TE_for_biorep
TE_for_biorep_data <- TE_13_for_biorep_clean_bio_data %>% 
  left_join(TE_11_for_strain_clean_bio_have_biorep_data , by=c("strain")) %>%
  mutate(TE=F_R_ratio/F_R_ratio_mean_bio) %>%
  filter(!is.na(TE))

#TE_for_strain
TE_for_strain_data <- TE_for_biorep_data %>%
  group_by(strain,biorep) %>%
  dplyr::summarise(sd_TE = sd(TE),TE = mean(TE)) %>%
  group_by(strain) %>%
  dplyr::summarise(sd_TE = sd(TE),TE_mean_bio = mean(TE)) %>%
  filter(strain != "A0756") 

TE_for_strain_have_biorep_data <- TE_for_strain_data %>%
  filter(!is.na(sd_TE))

TE_for_strain_only_one_biorep_data <- TE_for_strain_data %>%
  filter(is.na(sd_TE))

#save
TE_data_list <- list(TE_for_biorep_data=TE_for_biorep_data,
                     TE_for_strain_data=TE_for_strain_data,
                     TE_for_strain_have_biorep_data=TE_for_strain_have_biorep_data,
                     TE_for_strain_only_one_biorep_data=TE_for_strain_only_one_biorep_data)

set_dir <- "./"
save(TE_data_list,file = paste(set_dir,"5.TE_data_",Sys.Date(),".Rdata",sep = ""))

