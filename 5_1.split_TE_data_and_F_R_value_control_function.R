#function script
#####TE_F_R_value_filter #########
###function: F R value filter in 11/13 strains####
###input: data is the Raw TE data splited by different mutant_type (11/13),
#it contains strains F R value of the correspond mutant_type
###argument:is.unique tell function whether unique techrep and rename strain name to normal type 
TE_F_R_value_filter <- function(data,mutant_type,return.outlier = F,is.unique=F,F_11 = 1e+05 , R_11_illumiate = 0 , 
                                R_11 = 1e+08 , F_13=50, F_13_threshold = 5e+04,R_13_illumiate = 1e+05 , R_13_threshold = 1e+08){
  if(mutant_type==11){
    F_R_outlier <- data %>%
      filter(`F` < F_11 | `R` < R_11_illumiate |`R` > R_11 | `F` >= `R`) %>%
      mutate(problem_of_F = ifelse(`F` < F_11,"not illuminated",NA),
             problem_of_R_1 = ifelse( `R` < R_11_illumiate,"illuminated protein maybe mutation",NA),
             problem_of_R_2 = ifelse( `R` > R_11,"over threshold",NA),
             problem_of_F_R_1 =  ifelse(  `F` >= `R`,"F larger than R",NA))  
    
  }else if(mutant_type==13){
    F_R_outlier <- data %>%
      filter(`F`<=F_13 | `F` >= F_13_threshold | `R` < R_13_illumiate | `R` > R_13_threshold | `F` >= `R`) %>%
      mutate(problem_of_F =  ifelse(`F` <= F_13,"so small F",NA),
             problem_of_F_2 = ifelse(`F` >= F_13_threshold,"so large F",NA),
             problem_of_R_1 = ifelse(`R` < R_13_illumiate,"not illuminated",NA),
             problem_of_R_2 = ifelse( `R` > R_13_threshold,"over threshold",NA),
             problem_of_F_R_1 =  ifelse(  `F` >= `R`,"F larger than R",NA))
  }else{
    stop("not normal mutant_type")
  }
  
  if(return.outlier == T){
    if(is.unique == T){
      source("~/swj/Project/SI_TE_project/Function_dir/SI_TE_project_FUNCTION.R")
      F_R_outlier_unique <- change_strain_name_to_normal_strain_name(F_R_outlier) %>%
        group_by(strain,mutant_type,biorep,panel) %>%
        filter(row_number() == 1)
      return(F_R_outlier_unique)
    }else{
      return(F_R_outlier)
    }
  }else{
    clean_data <- data %>% anti_join(F_R_outlier)
    return(clean_data)
  }
}

#filter the only one techrep after TE_F_R_value_filter######
remove_only_one_techrep <- function(data,return.one.techrep = F){
  if(return.one.techrep == T){
    Output_data <-  data %>% group_by(strain,mutant_type,biorep,panel) %>%
      mutate(Group_size = n()) %>%
      filter(Group_size ==1) %>%
      select(-Group_size) %>%
      ungroup()
    return(Output_data)
  }else{
    Output_data <-  data %>% group_by(strain,mutant_type,biorep,panel) %>%
      mutate(Group_size = n()) %>%
      filter(Group_size >1) %>%
      select(-Group_size) %>%
      ungroup()
    return(Output_data)
  }
}

