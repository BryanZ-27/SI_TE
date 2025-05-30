###library
library(plyr)
library(dplyr)
library(tidyr)
library(purrr)
library(generics)
library(doParallel)


###function
Calculate_t <- function(E,Ed_cutoff,Alpha,e0){
  t <- log((Alpha*Ed_cutoff)/(Alpha*e0-e0+E)-Ed_cutoff/(Alpha*e0-e0+E)+E/(Alpha*e0-e0+E))/log(Alpha)
  return(t)
}


###main
##model simulation
Orgel_model_simulation_Argument_in_main <- expand_grid(
  Ed_cutoff = 0.06,
  Alpha     = 1.5,
  e0        = 0,
  El        = 5e-4,
  Em        = c(1e-3, 1.5e-3, 2e-3, 2.5e-3, 3e-3),
  P         = seq(from=0,to=0.8,by=0.2),
  E         = seq((2.5e-4),(3e-3),length.out = 30000),
) 

Ed_range <- c(0.03, 0.045, 0.06, 0.075, 0.09)
Alpha_range <- c(1.25, 1.375, 1.5, 1.625, 1.75)
El_range <- c(2.5e-4, 3.75e-4, 5e-4, 6.25e-4, 7.5e-4)

Orgel_model_simulation_Argument_in_Sup_Ed <- mclapply(Ed_range,function(Ed_value){
  Output_data <- Orgel_model_simulation_Argument_in_main %>% 
    mutate(Ed_cutoff = Ed_value)
  return(Output_data)
}) %>% bind_rows()
Orgel_model_simulation_Argument_in_Sup_Alpha <- mclapply(Alpha_range,function(Alpha_value){
  Output_data <- Orgel_model_simulation_Argument_in_main %>% 
    mutate(Alpha = Alpha_value)
  return(Output_data)
}) %>% bind_rows()
Orgel_model_simulation_Argument_in_Sup_El <- mclapply(El_range,function(El_value){
  Output_data <- Orgel_model_simulation_Argument_in_main %>% 
    mutate(El = El_value)
  return(Output_data)
}) %>% bind_rows()

Orgel_model_simulation_Argument <- bind_rows(Orgel_model_simulation_Argument_in_main,
                                             Orgel_model_simulation_Argument_in_Sup_Ed,
                                             Orgel_model_simulation_Argument_in_Sup_Alpha,
                                             Orgel_model_simulation_Argument_in_Sup_El) %>%
  group_by(Ed_cutoff,Alpha,e0,El,Em,P)

Orgel_model_simulation_result <-   
  mclapply(X = c(1:dim(Orgel_model_simulation_Argument)[1]),FUN = function(Argument_num){
    Argument_Data = Orgel_model_simulation_Argument[Argument_num,]
    Ed_cutoff = Argument_Data$Ed_cutoff
    Alpha     = Argument_Data$Alpha
    e0        = Argument_Data$e0
    E         = Argument_Data$E
    
    El        = Argument_Data$El
    Em        = Argument_Data$Em
    P         = Argument_Data$P
    
    t <- log((Alpha*Ed_cutoff)/(Alpha*e0-e0+E)-Ed_cutoff/(Alpha*e0-e0+E)+E/(Alpha*e0-e0+E))/log(Alpha)
    
    return(data.frame(E=E,t=t,
                      Ed_cutoff=Ed_cutoff,Alpha=Alpha,e0=e0,El=El,
                      Em=Em,P=P))
  },mc.cores = detectCores()-10) %>% rbind.fill()

#filter out model with NA or NaN
Orgel_model_simulation_result_filter_out <-  Orgel_model_simulation_result %>% 
  na.omit()

Orgel_model_simulation_result_filter_out_list <-  
  Orgel_model_simulation_result_filter_out %>%
  dplyr::group_by(Ed_cutoff,Alpha,e0,El,Em,P) %>%
  group_split()

Orgel_model_simulation_result_filter_out_under_El_Em_cutoff <- mclapply(Orgel_model_simulation_result_filter_out_list,FUN = function(Orgel_model_simulation_result_filter_out_data){
  Orgel_model_simulation_result_under_El_Em_cutoff <- Orgel_model_simulation_result_filter_out_data %>%
    filter(E<=Em,E>=El) %>% filter(P == 0)
  if (nrow(Orgel_model_simulation_result_under_El_Em_cutoff) != 0) {
    return(Orgel_model_simulation_result_under_El_Em_cutoff)
  }
},mc.cores = detectCores()-10)

Orgel_model_simulation_result_filter_out_under_El_Em_cutoff <- Orgel_model_simulation_result_filter_out_under_El_Em_cutoff[sapply(Orgel_model_simulation_result_filter_out_under_El_Em_cutoff, function(x) !is.null(x) && length(x) > 0)]
Orgel_model_simulation_result_filter_out_under_El_Em_cutoff <- c(Orgel_model_simulation_result_filter_out_under_El_Em_cutoff, 
                                                                 Orgel_model_simulation_result_filter_out_under_El_Em_cutoff, 
                                                                 Orgel_model_simulation_result_filter_out_under_El_Em_cutoff, 
                                                                 Orgel_model_simulation_result_filter_out_under_El_Em_cutoff, 
                                                                 Orgel_model_simulation_result_filter_out_under_El_Em_cutoff, 
                                                                 Orgel_model_simulation_result_filter_out_under_El_Em_cutoff, 
                                                                 Orgel_model_simulation_result_filter_out_under_El_Em_cutoff, 
                                                                 Orgel_model_simulation_result_filter_out_under_El_Em_cutoff, 
                                                                 Orgel_model_simulation_result_filter_out_under_El_Em_cutoff, 
                                                                 Orgel_model_simulation_result_filter_out_under_El_Em_cutoff)

Fake_datapoint_Data_unif <- mclapply(1:length(Orgel_model_simulation_result_filter_out_under_El_Em_cutoff),function(list_number){
  set.seed(list_number)
  Orgel_model_simulation_result_under_Ed_cutoff_Alpha_Em_cutoff <- Orgel_model_simulation_result_filter_out_under_El_Em_cutoff[[list_number]]
  
  Ed_cutoff = unique(Orgel_model_simulation_result_under_Ed_cutoff_Alpha_Em_cutoff$Ed_cutoff)
  Alpha     = unique(Orgel_model_simulation_result_under_Ed_cutoff_Alpha_Em_cutoff$Alpha)
  e0        = unique(Orgel_model_simulation_result_under_Ed_cutoff_Alpha_Em_cutoff$e0)
  El        = unique(Orgel_model_simulation_result_under_Ed_cutoff_Alpha_Em_cutoff$El)
  Em        = unique(Orgel_model_simulation_result_under_Ed_cutoff_Alpha_Em_cutoff$Em)
  P         = unique(Orgel_model_simulation_result_under_Ed_cutoff_Alpha_Em_cutoff$P)
  
  Fake_E <- approx(
    cumsum(Orgel_model_simulation_result_under_Ed_cutoff_Alpha_Em_cutoff$t)/sum(Orgel_model_simulation_result_under_Ed_cutoff_Alpha_Em_cutoff$t),
    Orgel_model_simulation_result_under_Ed_cutoff_Alpha_Em_cutoff$E,
    runif(400)
  )$y
  
  Fake_datapoint <- mclapply(1:length(Fake_E),FUN = function(rownumber){
    set.seed(rownumber)
    E_value <- Fake_E[rownumber]
    Fake_T_value <- runif(1,min = 0,max = Calculate_t(E = E_value,Ed_cutoff = Ed_cutoff,Alpha = Alpha,e0 = e0))
    return(data.frame(Fake_E=E_value,Fake_T=Fake_T_value+3))
  },mc.cores = 10) %>% bind_rows()
  
  T_threshold <- Calculate_t(E = Fake_E,
                             Ed_cutoff = Ed_cutoff,
                             Alpha = Alpha,
                             e0 = e0)
  
  Fake_datapoint_summary_data <- data.frame(Ed_cutoff=Ed_cutoff,Alpha=Alpha,e0=e0,
                                            El=El,Em=Em,P=P,
                                            Fake_datapoint,T_threshold=T_threshold, 
                                            seed_number = list_number)
  
  return(Fake_datapoint_summary_data)
  
},mc.cores = 10)

cor_filter_P_data_200 <- Fake_datapoint_Data_unif %>% bind_rows() %>%
  group_by(Ed_cutoff,Alpha,e0,El,Em,P,seed_number) %>%
  nest() %>%
  mutate(cor_res = purrr::map(data,.f = ~ generics::tidy(cor.test(~Fake_E + Fake_T,method = "spearm",data = .)))) %>%
  unnest(cor_res) 

#save
save(cor_filter_P_data_200, file = paste("../data/cor_filter_P_data_200_", Sys.Date(), ".Rdata", sep = ""))

