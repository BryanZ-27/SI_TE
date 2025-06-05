#function script

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

####marker wrong strain########
#function: remove wrong strains because of wrong markers
ErrorStrains <- c("11C10",
                  "11H06",
                  "11H03",
                  "10A04",
                  "10A05",
                  "10A07",
                  "10C04",
                  "10G06",
                  "08A05",
                  "08F08",
                  "08F10",
                  "08F12",
                  "08G03",
                  "08G05",
                  "08G07",
                  "08G09",
                  "08G10",
                  "08G11",
                  "08G12",
                  "08H04",
                  "08H10",
                  "08H11",
                  "08H12",
                  "04B09",
                  "04C01",
                  "04D01",
                  "04E05",
                  "04G01",
                  "04C07",
                  "04C10",
                  "04E02",
                  "04E05",
                  "04E12",
                  "04F04",
                  "04G07",
                  "04H05",
                  "04H09",
                  "01A06",
                  "01B03",
                  "01D12",
                  "01E08",
                  "01B06",
                  "01C04",
                  "01G08",
                  "01H05")

ErrorStrains_rename <- function(ErrorStrains){
  strainHash <- data.frame(locationLine=c("A","B","C","D","E","F","G","H"),Value=c(0:7))
  ErrorStrains_vector<- ErrorStrains %>% sapply(function(i){
    x <-i %>%
      str_match("(\\d+)(\\D)(\\d+)")
    c <- as.character(strainHash[strainHash$locationLine == x[3],]$Value*12 + as.numeric(x[4])) %>%
      str_pad(width=2,side = "left",pad = "0")
    
    tt <- i %>%
      str_replace("\\D\\d+",as.character(c)) %>%
      str_pad(width = 5,side = "left",pad = "A")
    return(tt)
  })
  return(ErrorStrains_vector)
}

ErrorStrains_vector <- ErrorStrains_rename(ErrorStrains = ErrorStrains)

####Categorical_variable ###
Categorical_variable <- c("strains_group","rep_location_name","rep_location_num",
                          "tube_type","tube_type_YPD","tube_type_SC",
                          "glu_concentration")

#filter wild type strain and rename############
Filter_wild_type_strain_and_rename <- function(SI_for_strain_or_biorep_data, wild_type_name,BY4741_strains_group=c()){
  if(wild_type_name == "B"){
    BY_data <- SI_for_strain_or_biorep_data %>% 
      filter(str_detect(strains,"^B")) %>%
      mutate(strains=ifelse(strains_group %in% BY4741_strains_group , "BY4741","BY4716"))
    return(BY_data)
  }else if(wild_type_name == "R"){
    RM_data <-  SI_for_strain_or_biorep_data %>% 
      filter(str_detect(strains,"^R")) 
    return(RM_data)
  }
}

#filter out wild type strain in SI biorep and strain data################
SI_wild_type_remove <- function(SI_for_strain_or_biorep_data , String){
  wild_type_remove_data <- SI_for_strain_or_biorep_data %>%
    filter(!str_detect(strains,String))
  return(wild_type_remove_data)
}

#Filter or filter out glu concentration == 0.5%####################
Filter_glu_concentration <- function(SI_for_strain_or_biorep_data,is.filter.out = T){
  if(is.filter.out==F){
    Outputdata <- SI_for_strain_or_biorep_data %>%
      filter(glu_concentration == 0.5)
  }else{
    Outputdata <- SI_for_strain_or_biorep_data %>%
      filter(glu_concentration != 0.5)
  }
 return(Outputdata)
}

#remove duplicate and retain the tube_type == "C" or "J"####################
remove_duplicate_reconstruct <- function(SI_for_biorep_or_strain_data,is.biorep.data = T ,is.filter.out = T,
                                         Group_var="strains_group"){
  if(is.biorep.data == T){
    SI_for_biorep_or_strain_data <- SI_for_biorep_or_strain_data %>% group_by_at(c("strains","biorep",Group_var))
  }else{
    SI_for_biorep_or_strain_data <- SI_for_biorep_or_strain_data %>% group_by_at(c("strains",Group_var))
  }
  SI_rm_dp <- SI_for_biorep_or_strain_data %>%
    filter(!tube_type_SC %in% c("G","Y")) %>%
    filter(row_number()==1)
  
  SI_duplicate <- SI_for_biorep_or_strain_data %>% anti_join(SI_rm_dp)
  
  if(is.filter.out == T){
    return(SI_rm_dp)
  }else{
    return(SI_duplicate)
  }
}

remove_duplicate_strain_group_reconstruct <- function(SI_for_biorep_or_strain_data,is.biorep.data = T ,
                                                      is.filter.out = T,Group_var="strains_group"){
  if(is.biorep.data == T){
    SI_for_biorep_or_strain_data <- SI_for_biorep_or_strain_data %>% group_by_at(c("strains","biorep",Group_var))
  }else{
    SI_for_biorep_or_strain_data <- SI_for_biorep_or_strain_data %>% group_by_at(c("strains",Group_var))
  }
  SI_rm_dp <- SI_for_biorep_or_strain_data %>%
    filter(!tube_type_SC %in% c("G","Y")) %>%
    filter(row_number()==1 & strains !="A0606") %>%
    bind_rows(SI_for_biorep_or_strain_data %>% filter(strains =="A0606" & strains_group == "6A1-A12(6A1-A3rep)"))
  
  SI_duplicate <- SI_for_biorep_or_strain_data %>% anti_join(SI_rm_dp)
  
  if(is.filter.out == T){
    return(SI_rm_dp)
  }else{
    return(SI_duplicate)
  }
}

remove_duplicate <- function(SI_for_biorep_data,SI_for_strain_data,is.filter.out = T,Categorical_variable){
  SI_for_strain_to_semi_join <- SI_for_strain_data%>% group_by(strains) %>%
    mutate(group_size = n()) %>%
    filter(tube_type_SC %in% c("C","J"))%>%  #get tube_type_SC == "C" or "J"
    mutate(group_size = n()) %>%
    mutate(rm_J_dp = ifelse(group_size >1 & tube_type_SC %in% c("J","G","Y"),1,0)) %>% 
    ##remove the J tube_type when group_size == 2 , contatin all C tube_type_SC
    filter(rm_J_dp != 1)

  SI_biorep_rm_dp <- SI_for_biorep_data %>% 
    semi_join(SI_for_strain_to_semi_join,by=c("strains",
                                              Categorical_variable[!Categorical_variable %in% 
                                                                     c("rep_location_name","rep_location_num","tube_type")]))

  SI_biorep_duplicate <- SI_for_biorep_data %>% 
    anti_join(SI_biorep_rm_dp,by=c("strains","biorep",
                                   Categorical_variable[!Categorical_variable %in% 
                                                          c("rep_location_name","rep_location_num","tube_type")]))
  if(is.filter.out == T){
    return(SI_biorep_rm_dp)
  }else{
    return(SI_biorep_duplicate)
  }
}

