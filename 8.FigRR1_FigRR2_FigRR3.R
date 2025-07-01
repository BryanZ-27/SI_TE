#Figure RR1
#set
saveDir <- "./"
#library
library(dplyr)
library(ggplot2)
library(stringr)


#load
"./5.TE_techrep_data2.Rdata"


#main
TE_techrep_data <- TE_11_clean_F_R_data_rm_o1_techrep %>% 
  merge(TE_13_clean_F_R_data_rm_o1_techrep, by = c("strain", "biorep", "techrep")) %>% filter(strain != "BY4716") %>% 
  mutate(F_R_ratio_11 = F.x/R.x, F_R_ratio_13 = F.y/R.y, TE = F_R_ratio_13/F_R_ratio_11) %>% 
  group_by(strain, biorep) %>% summarise(CV = sd(TE)/mean(TE), mean_TE = mean(TE)) %>% filter(biorep != 4) %>% 
  ungroup() %>% group_by(strain) %>% mutate(n = n(), mean_CV = mean(CV)) %>% filter(n > 2) # %>% filter(strain %in% t$strain_rename)

TE_matrix_data <- acast(TE_techrep_data, strain ~ biorep, value.var = "CV") %>% as_tibble() %>% 
  mutate(na_count = rowSums(is.na(select(., everything())))) %>% filter(na_count < 1) %>% select(1:3) %>% as.matrix()

cv_plot <- ggplot(TE_techrep_data, aes(x = biorep, y = reorder(strain, mean_CV), fill = CV)) +
  geom_tile(linewidth = 0.5) +
  #geom_text(aes(label = round(CV, 2)), color = "black", size = 5) +
  scale_fill_gradient(low = "yellow", high = "red", 
                      limits = c(0, 1.5),
                      name = "CV among technical replicates") +
  labs(x = "Biological replicate", y = "Strain") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 10, face = "bold"),
    axis.line = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8), 
    legend.position = "top"
  )
cv_plot

TE_techrep_data_perBiorep <- TE_techrep_data %>% select(1:3) %>% 
  tidyr::pivot_wider(id_cols = "strain", names_from = biorep, values_from = CV, names_prefix = "CV_biorep")

cor.test(TE_techrep_data_perBiorep$CV_biorep1, TE_techrep_data_perBiorep$CV_biorep2) #R = 0.2982756, p-value = 3.858e-07
cor.test(TE_techrep_data_perBiorep$CV_biorep1, TE_techrep_data_perBiorep$CV_biorep3) #R = 0.2596794, p-value = 1.115e-05

cv_plot_1_2 <- ggplot(TE_techrep_data_perBiorep, aes(CV_biorep1, CV_biorep2)) +
  geom_point(size = 0.7) +
  geom_smooth(method = "lm", se = F) +
  annotate("text", x = 0.8, y = 1.55, parse = T, size = 3.5,
           label = expression(paste("Pearson's ", italic("R"), " = 0.30, ", italic("P"), " = 3.86×10"^"-7"))) +
  scale_x_continuous(limits = c(0, 1.5)) +
  scale_y_continuous(limits = c(0, 1.55)) +
  labs(x = "CV among technical replicates \nof biological replicate 1", 
       y = "CV among technical replicates \n of biological replicate 2") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 10, face = "bold"),
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8)
  )
cv_plot_1_2

cv_plot_1_3 <- ggplot(TE_techrep_data_perBiorep, aes(CV_biorep1, CV_biorep3)) +
  geom_point(size = 0.7) +
  geom_smooth(method = "lm", se = F) +
  annotate("text", x = 0.8, y = 1.55, parse = T, size = 3.5,
           label = expression(paste("Pearson's ", italic("R"), " = 0.26, ", italic("P"), " = 1.12×10"^"-5"))) +
  scale_x_continuous(limits = c(0, 1.5)) +
  scale_y_continuous(limits = c(0, 1.55)) +
  labs(x = "CV among technical replicates \nof biological replicate 1", 
       y = "CV among technical replicates \nof biological replicate 3") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 10, face = "bold"),
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8)
  )
cv_plot_1_3


#save
ggsave(paste0(saveDir, "cv_plot_", Sys.Date(), ".pdf"), cv_plot, units = "cm", width = 10, height = 16)
ggsave(paste0(saveDir, "cv_plot_1_2_", Sys.Date(), ".pdf"), cv_plot_1_2, units = "cm", width = 8, height = 8)
ggsave(paste0(saveDir, "cv_plot_1_3_", Sys.Date(), ".pdf"), cv_plot_1_3, units = "cm", width = 8, height = 8)


#Figure RR2
# library######
library(plyr)
library(dplyr)
library(stringr)
library(reshape2)
library(generics)
library(purrr)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(doParallel)
library(parallel)
# source function

# load data

#MAIN#######
Calculate_t <- function(E,Ed_cutoff,Alpha,e0){
  t <- log((Alpha*Ed_cutoff)/(Alpha*e0-e0+E)-Ed_cutoff/(Alpha*e0-e0+E)+E/(Alpha*e0-e0+E))/log(Alpha)
  return(t)
}

Orgel_model_simulation_Argument_in_main <- expand_grid(
  Ed_cutoff = 0.06,
  Alpha     = 1.5,
  e0        = 0,
  El        = 5*10^-4,
  Em        = c(1*10^-3,1.5*10^-3,2*10^-3,3*10^-3),
  P         = seq(from=0,to=0.8,by=0.2),
  E         = seq((1*10^-4),(3*10^-3),length.out = 30000),
) 

Ed_range <- c(0.01,0.02,0.04,0.08,0.10)
Alpha_range <- c(1.01,1.2,1.8,2.0)
El_range <- c(1*10^-4,3*10^-4,7*10^-4,9*10^-4)

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

Orgel_model_simulation_result_filter_out %>%
  filter(Alpha %in% c(1.01,1.5)) %>%
  filter(Ed_cutoff >0.0024 , Ed_cutoff < 0.0026 | Ed_cutoff >0.0084 , Ed_cutoff < 0.0086) %>%
  filter(t>1)%>%
  ggplot()+
  geom_line(aes(x=t,y=E,group=interaction(Ed_cutoff,Alpha,e0)))+
  facet_grid(Alpha~Ed_cutoff)

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

Fake_datapoint_Data_Rnorm <- mclapply(1:length(Orgel_model_simulation_result_filter_out_under_El_Em_cutoff),function(list_number){
  set.seed(list_number)
  Orgel_model_simulation_result_under_Ed_cutoff_Alpha_Em_cutoff <- Orgel_model_simulation_result_filter_out_under_El_Em_cutoff[[list_number]]
  
  
  Ed_cutoff = unique(Orgel_model_simulation_result_under_Ed_cutoff_Alpha_Em_cutoff$Ed_cutoff)
  Alpha     = unique(Orgel_model_simulation_result_under_Ed_cutoff_Alpha_Em_cutoff$Alpha)
  e0        = unique(Orgel_model_simulation_result_under_Ed_cutoff_Alpha_Em_cutoff$e0)
  
  El        = unique(Orgel_model_simulation_result_under_Ed_cutoff_Alpha_Em_cutoff$El)
  Em        = unique(Orgel_model_simulation_result_under_Ed_cutoff_Alpha_Em_cutoff$Em)
  P         = unique(Orgel_model_simulation_result_under_Ed_cutoff_Alpha_Em_cutoff$P)
  
  FakeApprox <- approx(
    cumsum(Orgel_model_simulation_result_under_Ed_cutoff_Alpha_Em_cutoff$t)/sum(Orgel_model_simulation_result_under_Ed_cutoff_Alpha_Em_cutoff$t),
    Orgel_model_simulation_result_under_Ed_cutoff_Alpha_Em_cutoff$E,
    runif(460))
  Fake_E <- FakeApprox$y
  FakeEsT <- FakeApprox$x
  FakeRnorm <- rnorm(1000, mean = mean(FakeEsT), sd = sd(FakeEsT))
  if (min(FakeRnorm) < 0) {FakeRnorm <- FakeRnorm + abs(min(FakeRnorm))}
  maxT <- Calculate_t(E = Fake_E,Ed_cutoff = Ed_cutoff,Alpha = Alpha,e0 = e0)+4
  
  Fake_datapoint <- mclapply(1:length(Fake_E),FUN = function(rownumber){
    set.seed(rownumber)
    E_value <- Fake_E[rownumber]
    if (is.na(E_value)) {Fake_T_value <- -1}else {
      FakeRnormMod <- FakeRnorm/mean(FakeRnorm)*5
      Fake_T_value <- rnorm(100, mean = mean(FakeRnormMod)+3, sd = sd(FakeRnormMod))
      if(Fake_T_value > maxT[rownumber]) {Fake_T_value <- -1}
    }
    return(data.frame(Fake_E=E_value,Fake_T=Fake_T_value))
  },mc.cores = 10) %>% bind_rows()
  
  Fake_datapoint_summary_data <- data.frame(Ed_cutoff=Ed_cutoff,Alpha=Alpha,e0=e0,
                                            El=El,Em=Em,P=P,
                                            Fake_datapoint,T_threshold=maxT,seed_number = list_number) %>% 
    filter(Fake_T > 0)
  
  return(Fake_datapoint_summary_data)
  
},mc.cores = 10)


cor_filter_P_data_200_1 <- Fake_datapoint_Data_Rnorm %>% bind_rows() %>%
  group_by(Ed_cutoff,Alpha,e0,El,Em,P,seed_number) %>%
  nest() %>%
  mutate(cor_res = purrr::map(data,.f = ~ generics::tidy(cor.test(~Fake_E + Fake_T,method = "spearm",data = .)))) %>%
  unnest(cor_res) 

#library
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(dplyr)
#source function

#load data#####

#main
fig1d_data_original <- cor_filter_P_data_200_1 %>% 
  ungroup() %>% 
  filter(Alpha==1.5,Ed_cutoff == 0.06,El==5*10^-4, Em == 2e-3) %>% 
  select(1:8)

fig1d_data_final <- tibble(Ed_cutoff = 0.06,Alpha=1.5,e0 = 0,El=5*10^-4, Em = 2e-3, P = 0, 
                           seed_number = NA, rho = NA)
myThreshold <- seq(3.78, 12.59, length.out = 6)
for(i in 1 : nrow(fig1d_data_original)){
  df <- fig1d_data_original$data[i][[1]]
  for (value_threshold in myThreshold[1:5]) {
    new_data <- df %>% filter(Fake_T >= value_threshold)
    new_ratio <- tibble(Ed_cutoff = 0.06,Alpha=1.5,e0 = 0,El=5*10^-4, Em = 2e-3, P = value_threshold, seed_number = i, 
                        rho = cor.test(new_data$Fake_T, new_data$Fake_E, method = "spearman")$estimate)
    fig1d_data_final <- rbind(fig1d_data_final, new_ratio) %>% filter(!is.na(rho))
  }
}
fig1e_plotData <- Fake_datapoint_Data_Rnorm %>%
  bind_rows() %>% filter(Ed_cutoff == 0.06 & Alpha == 1.5 & El == 5e-4 & Em == 2e-3 & P == 0)

fig1d_data_final <- fig1d_data_final %>% mutate(p.value = pvalue)

pvalue <- rep(c(0.4479, 0.5116, 0.7502, 0.003644, 0.0004993), 10)

fig1d_data <- fig1d_data_final %>% 
  mutate(adjust_rho = -rho*10) %>% 
  group_by(P) %>%
  dplyr::summarise(mean_rho = mean(rho), sd_rho = sd(rho), 
                   adjust_mean_rho = mean(adjust_rho)+2+0.5, adjust_sd_rho = sd(adjust_rho), 
                   mean_p.value = ifelse(mean(p.value)*5 > 1, 1, mean(p.value)*5), sd_p.value = sd(p.value), 
                   fisher_P = pchisq(-2 * sum(log(p.value*5)), df = 2 * length(p.value*5), lower.tail = FALSE)) %>%
  ungroup() %>% 
  mutate(adjust_mean_p.value = log10(mean_p.value)*2-2-0.5)

fig1d_color <- c("#99cc99", "#66cc66", "#669933", "#339933", "#336633")

#4 5.5 7 8.5 10

FIG.1.E <- fig1d_data %>% ggplot()+
  geom_bar(aes(x=P,y=adjust_mean_rho, fill = as.factor(P)), 
           stat = "identity",position =position_dodge(width = 0.15),width = 0.10)+
  geom_bar(aes(x=P,y=adjust_mean_p.value, fill = as.factor(P)), 
           stat = "identity",position =position_dodge(width = 0.15),width = 0.01)+
  geom_errorbar(aes(x = P, ymin = adjust_mean_rho - adjust_sd_rho, ymax = adjust_mean_rho + adjust_sd_rho), 
                stat = "identity",position = position_dodge(width = 0.15),width = 0.10, color = "black") + 
  geom_point(aes(x=P,y=adjust_mean_p.value, fill = as.factor(P)), 
             position = position_dodge(width = 0.15), size=30, color = fig1d_color, shape=21)+
  geom_abline(slope = 0,intercept=log10(0.05)*2-2-0.5,color="red",linewidth=1.2,linetype=2)+
  annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymax = 0.5, ymin = -0.5,fill = "white") +
  annotate(geom = "text", y = 0, x = fig1d_data$P %>% unique(), label = fig1d_data$P %>% unique(), size = 4.5) +
  annotate(geom = "text",label="alpha==1.5",x = 4.1,y = 9,size=4.5,parse = TRUE)+
  annotate(geom = "text",label="~italic(D)==0.06",x = 4.1,y = 8,size=4.5,parse = TRUE)+
  annotate(geom = "text",label="~~~~~italic(U)==2%*%10^-3",x = 4.1,y = 7,size=4.5,parse = TRUE)+
  annotate(geom = "text",label="~~~~~italic(L)==5%*%10^-4",x = 4.1,y = 6,size=4.5,parse = TRUE)+
  ylab("Spearman’s correlation between translational\nerror rate and lifespan\nP value(log10 scale)\trho")+
  ###plot modify###
  scale_fill_manual(values = fig1d_color)+
  scale_color_manual(values = fig1d_color)+
  scale_y_continuous(limits = c(-9, 9), 
                     breaks=c(-8.5, -6.5, -4.5, -2.5, -0.5, 0.5, 2.5, 4.5, 6.5, 8.5), 
                     labels = c(expression(10^-3), expression(10^-2), 
                                expression(10^-1), 0, expression(10^1), 
                                "0.2", 0, "-0.2", "-0.4", "-0.6"))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5,size = 18, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title = element_text(size=16,face = "bold"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank()
  )
FIG.1.E

ggsave(paste("./Fig.RR2_", Sys.Date(), ".pdf", sep = ""), FIG.1.E, width = 32, height = 24, units = "cm")

#zoom
zoom1_data <- fig1d_data_original %>% 
  filter(seed_number == 618) %>% 
  pull(data) %>% as.data.frame()

zoom1_plot <- zoom1_data %>% ggplot()+
  geom_point(aes(x=Fake_E,y=Fake_T), size = 1)+
  #stat_cor(aes(x=Fake_E,y=Fake_T),method = "spearman",label.x = 8e-4, label.y = 10,cor.coef.name = "rho")+
  geom_smooth(aes(x=Fake_E,y=Fake_T),method = "lm",se = F, color = fig1d_color[1])+
  scale_x_continuous(limits = c(0.0005, 0.0020), breaks = c(seq(0.0005, 0.0020, 0.0005)), labels = c(seq(5, 20, 5)))+
  scale_y_continuous(breaks = c(seq(0, 10, 5)), labels = c(seq(0, 10, 5)), limits = c(4.5, 12))+
  labs(x="",
       y="Lifespan")+
  theme_classic()+theme(plot.title = element_text(hjust = 0.5,size = 18, face = "bold"),
                        axis.text=element_text(size=14,face = "bold"),
                        axis.title = element_text(size=16,face = "bold"),
                        axis.text.x = element_text(vjust = 1,hjust = 1),
                        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
                        # legend.position = c(.1,.85),legend.background = element_rect(colour = "black",size = 1.2,fill = "lightyellow"),legend.title = element_text(size=12,face = "bold"))
  )
zoom1_plot

zoom2_data <- fig1d_data_original %>% 
  filter(seed_number == 618) %>% 
  pull(data) %>% as.data.frame() %>% filter(Fake_T >= 7.304)

zoom2_plot <- zoom2_data %>% ggplot()+
  geom_point(aes(x=Fake_E,y=Fake_T), size = 1)+
  #stat_cor(aes(x=Fake_E,y=Fake_T),method = "spearman",label.x = 8e-4, label.y = 10,cor.coef.name = "rho")+
  geom_smooth(aes(x=Fake_E,y=Fake_T),method = "lm",se = F, color = fig1d_color[3])+
  scale_x_continuous(limits = c(0.0005, 0.0020), breaks = c(seq(0.0005, 0.0020, 0.0005)), labels = c(seq(5, 20, 5)))+
  scale_y_continuous(breaks = c(seq(0, 10, 5)), labels = c(seq(0, 10, 5)), limits = c(4.5, 12))+
  labs(x="",
       y="Lifespan")+
  theme_classic()+theme(plot.title = element_text(hjust = 0.5,size = 18, face = "bold"),
                        axis.text=element_text(size=14,face = "bold"),
                        axis.title = element_text(size=16,face = "bold"),
                        axis.text.x = element_text(vjust = 1,hjust = 1),
                        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
                        # legend.position = c(.1,.85),legend.background = element_rect(colour = "black",size = 1.2,fill = "lightyellow"),legend.title = element_text(size=12,face = "bold"))
  )
zoom2_plot

zoom3_data <- fig1d_data_original %>% 
  filter(seed_number == 618) %>% 
  pull(data) %>% as.data.frame() %>% filter(Fake_T >= 10.828)

zoom3_plot <- zoom3_data %>% ggplot()+
  geom_point(aes(x=Fake_E,y=Fake_T), size = 1)+
  #stat_cor(aes(x=Fake_E,y=Fake_T),method = "spearman",label.x = 8e-4, label.y = 10,cor.coef.name = "rho")+
  geom_smooth(aes(x=Fake_E,y=Fake_T),method = "lm",se = F, color = fig1d_color[5])+
  scale_x_continuous(limits = c(0.0005, 0.0020), breaks = c(seq(0.0005, 0.0020, 0.0005)), labels = c(seq(5, 20, 5)))+
  scale_y_continuous(breaks = c(seq(0, 10, 5)), labels = c(seq(0, 10, 5)), limits = c(4.5, 12))+
  labs(x="",
       y="Lifespan")+
  theme_classic()+theme(plot.title = element_text(hjust = 0.5,size = 18, face = "bold"),
                        axis.text=element_text(size=14,face = "bold"),
                        axis.title = element_text(size=16,face = "bold"),
                        axis.text.x = element_text(vjust = 1,hjust = 1),
                        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
                        # legend.position = c(.1,.85),legend.background = element_rect(colour = "black",size = 1.2,fill = "lightyellow"),legend.title = element_text(size=12,face = "bold"))
  )
zoom3_plot

###save
ggsave(paste("./Fig.RR2.1_", Sys.Date(), ".pdf", sep = ""), zoom1_plot, width = 8, height = 6, units = "cm")
ggsave(paste("./Fig.RR2.2_", Sys.Date(), ".pdf", sep = ""), zoom2_plot, width = 8, height = 6, units = "cm")
ggsave(paste("./Fig.RR2.3_", Sys.Date(), ".pdf", sep = ""), zoom3_plot, width = 8, height = 6, units = "cm")


#Figure RR3
# library
library(tidyr)
library(purrr)
library(doParallel)
library(generics)
library(broom)

# source function
Transform_strain_name <- function(strainname,name_type){
  strainname_matrix <- matrix(data = c(1:96),nrow = 8,ncol = 12,byrow = TRUE,
                              dimnames = list(c("A","B","C","D","E","F","G","H"),c(1:12)))
  s <- str_split(strainname,pattern = "-",simplify = TRUE)
  plate <- as.numeric(s[,1])
  RN <- str_extract(s[,2],pattern = "\\D")
  CN <- as.numeric(str_extract(s[,2],pattern = "\\d+"))
  strain_number <- c()
  for (i in 1 : length(RN)) {
    strain_number <- c(strain_number, strainname_matrix[RN[i], CN[i]])
  }
  output_strainname <- paste0("A",str_pad(plate,width = 2,pad = "0"),str_pad(strain_number,width = 2,pad = "0"))
  return(output_strainname)
}

# load data#####
#1# SI data#####
load("./3.SI_data.Rdata")

SI_804_data <- SI_rm_wt_gc_dp_dpsg_wm_ol15_list$SI_for_strain_clean_SI_biorep_ori_and_curve_rm_wt_gc_dp_dpsg_wm_ol15_data_output_list$SI_for_strain_clean_SI_biorep_ori_and_curve_rm_na_sd_day2_min_0.2_0.5_shift_0.3_day2

#2# TE data#####
TE_env <- new.env()
load("./5.TE_data.Rdata",envir = TE_env)
load("./5.TE_11_biorep_data.Rdata",envir = TE_env)
load("./5.TE_13_biorep_data.Rdata",envir = TE_env)

TE_11_for_strain_data <- TE_env$TE_11_biorep_clean_list$TE_11_for_strain_clean_bio_have_biorep_data %>%
  dplyr::rename(sd_F_R_ratio_11=sd_F_R_ratio,F_R_ratio_mean_bio_11=F_R_ratio_mean_bio)

TE_13_for_strain_data <- TE_env$TE_13_biorep_clean_list$TE_13_for_strain_clean_bio_data %>%
  filter(!is.na(sd_F_R_ratio))%>%
  dplyr::rename(sd_F_R_ratio_13=sd_F_R_ratio,F_R_ratio_mean_bio_13=F_R_ratio_mean_bio)

#356 observation
TE_271_and_11_13_data <- TE_env$TE_data_list$TE_for_strain_have_biorep_data %>%
  full_join(TE_11_for_strain_data) %>%
  full_join(TE_13_for_strain_data, by = c("strain")) %>%
  filter(!strain %in% c("BY4716","10-G12","RM11")) %>% #remove the wrong strain in the first batch and RM11 strain .
  filter(!strain %in% c("1-A10","4-C7","4-G7","4-E8","10-G12","1-G6")) %>% #remove the wrong strain with wrong genoytpe
  rowwise() %>%
  mutate(strain_rename=Transform_strain_name(strain))

#353 observation
TE_271_and_11_13_RM3point <- TE_271_and_11_13_data %>% 
  filter(!strain_rename %in% c("A1052","A1085","A1044"))

#345 observation
TE_271_and_11_13_RM11point <- TE_271_and_11_13_RM3point %>%
  filter(!strain_rename %in% c("A0141","A1035","A1029","A1024")) %>%
  filter(!strain_rename %in% c("A0190","A0645","A0863","A0306")) 

# MAIN#######
#1.Construct TE_SI_data (use TE_data left join SI_data,remove NA)######
TE_SI_data <- TE_271_and_11_13_RM11point %>%
  left_join(SI_804_data,by = c("strain_rename"="strains")) %>%
  na.omit() %>%
  ungroup()

#2.Construct list with different remove SI fraction.#########
#@## 2.1.cor list#####
SI_TE_cor_diff_SI_cutoff_output_list <-  mclapply(mc.cores = 60, X = c(seq(from=1,to=0.1,by=-0.1)) ,FUN = function(cutoff){
  output <- TE_SI_data %>% 
    dplyr::mutate(N=n()) %>%
    arrange(SI_mean_bio) %>%
    filter(row_number() <=N*cutoff) %>%
    nest() %>%
    mutate(res = purrr::map(data, .f = ~ tidy(cor.test(~TE_mean_bio +SI_mean_bio, data = .,method="spearman")))) %>%
    unnest(res) %>%
    mutate(observation_number = purrr::map(data,.f = function(data){nrow(data)})) %>%
    mutate(cut_off = cutoff) %>%
    unnest(observation_number)
  return(output)
})

#@## 2.2 point list correspond to cor list###
SI_TE_point_plot_diff_SI_cutoff_output_list <- mclapply(X = c(seq(from=1,to=0.1,by=-0.1)) ,FUN = function(cutoff){
  output <- TE_SI_data %>% 
    dplyr::mutate(N=n()) %>%
    arrange(SI_mean_bio) %>%
    filter(row_number() <=N*cutoff) %>%
    mutate(cut_off = cutoff)
  return(output)
})

#3.Plot obj###########
FIG4.A_For_plot_ob <- SI_TE_point_plot_diff_SI_cutoff_output_list %>% bind_rows()
FIG4.B_For_plot_ob <- SI_TE_cor_diff_SI_cutoff_output_list %>% bind_rows()

#4.Plot #######
#@###4.1 Figure 4A#######
new_fig4a_data <- FIG4.A_For_plot_ob %>% 
  filter(cut_off != "0.2") %>% 
  filter(cut_off != "0.4") %>% 
  filter(cut_off != "0.6") %>% 
  filter(cut_off != "0.8")

fig4_color <- c("#9cd49a", "#73c371", "#4daf4a", "#3b8639", "#295c27", "#163316")

FIG4.A <- new_fig4a_data %>% ggplot()+
  geom_point(aes(x=TE_mean_bio,y=SI_mean_bio,group=cut_off), size = 0.3)+
  facet_wrap(~cut_off,scales = "free_y")+
  stat_smooth(aes(x=TE_mean_bio,y=SI_mean_bio,group=cut_off), formula = "y ~ x", method = "lm",se = F,color='#80daeb')+
  scale_x_continuous(limits = c(6.5e-4, 13e-4), breaks = seq(0.0008, 0.0012, 0.0002), labels = c("8", "10", "12")) + 
  scale_y_continuous(limits = c(1, 13), breaks = seq(3, 12, 3), labels = seq(3, 12, 3))+
  labs(x="Translational error rate (e) (e-4)",y="Lifespan")+
  theme_classic()+
  theme(axis.text=element_text(size=9,face = "bold"),
        axis.title = element_text(size=10,face = "bold"),
        axis.text.x = element_text(vjust = 1,hjust = 1),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        legend.position = "none", 
        strip.text = element_text(size=8,color = "red",face = "bold"))
FIG4.A

#@###4.2 Figure 4B#######
shift <- 0.5
new_fig4b_data <- FIG4.B_For_plot_ob %>% 
  filter(cut_off != "0.2") %>% 
  filter(cut_off != "0.4") %>% 
  filter(cut_off != "0.6") %>% 
  filter(cut_off != "0.8") %>% 
  mutate(order = c(seq(0, 1, 0.2)))

#bootstrap-based SD
bootstrapSD <- c()
bootstrapMean <- c()
for (i in 1:nrow(new_fig4b_data)) {
  myData <- new_fig4b_data$data[[i]] %>% select(TE_mean_bio, SI_mean_bio) %>% t() %>% as_tibble()
  myCorr <- c()
  for (t in 1:100) {
    set.seed(i*100+t)
    myBootstrap <- sample(myData, ncol(myData), replace = T) %>% t() %>% as_tibble()
    myCorr <- c(myCorr, cor.test(myBootstrap$V1, myBootstrap$V2, method = "spearman")$estimate)
  }
  bootstrapSD <- c(bootstrapSD, sd(myCorr))
  bootstrapMean <- c(bootstrapMean, mean(myCorr))
}
new_fig4b_data <- new_fig4b_data %>% mutate(bootstrap_mean = bootstrapMean, bootstrap_SD = bootstrapSD)

#add bonferroni p
new_fig4b_data <- new_fig4b_data %>% mutate(adjust_p = p.adjust(p.value, method = "bonferroni"))

FIG4.B_bar_plot <- new_fig4b_data %>% ggplot()+
  geom_bar(aes(x=order,y=bootstrap_mean*20+4+shift, fill = as.factor(cut_off)), 
           stat = "identity",position =position_dodge(width = 0.15),width = 0.12)+
  geom_bar(aes(x=order,y=log10(adjust_p)*4-4-shift, fill = as.factor(cut_off)), 
           stat = "identity",position =position_dodge(width = 0.15),width = 0.02)+
  geom_errorbar(aes(x = order, ymax = (bootstrap_mean + bootstrap_SD)*20+4+shift, 
                    ymin = ifelse((bootstrap_mean - bootstrap_SD)*20+4+shift >= 0.55, 
                                  (bootstrap_mean - bootstrap_SD)*20+4+shift, 0.55)), 
                stat = "identity",position = position_dodge(width = 0.15),width = 0.08, color = "black") + 
  geom_point(aes(x=order,y=log10(adjust_p)*4-4-shift, fill = as.factor(cut_off)), 
             position = position_dodge(width = 0.15), size=18, color = fig4_color, shape=21)+
  geom_abline(slope = 0,intercept=log10(0.05)*4-4-shift,color="red",linewidth=1.2,linetype=2)+
  annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymax = shift, ymin = -shift,fill = "white") +
  annotate(geom = "text", y = 0, x = new_fig4b_data$order %>% unique(), 
           label = new_fig4b_data$cut_off %>% unique(), size = 9/.pt, color = "black", fontface = "bold") +
  scale_y_continuous(limits = c(-15, 15), 
                     breaks=c(-12-shift, -8-shift, -4-shift, 0-shift, 0+shift, 4+shift, 8+shift ,12+shift), 
                     labels = c(expression(10^-2), expression(10^-1), 0, expression(10^1), 
                                "-0.2", "0", "0.2", "0.4"))+
  scale_x_continuous(breaks = c(seq(0, 1, 0.2)), labels = c("0%", "10%", "30%", "50%", "70%", "90%"))+
  scale_fill_manual(values = fig4_color)+
  scale_color_manual(values = fig4_color)+
  ylab("Spearman’s correlation between translational\nerror rate and lifespan\nlog10 P\trho")+
  theme_classic()+
  theme(axis.text=element_text(size=9,face = "bold", color = "black"),
        axis.title = element_text(size=10,face = "bold", color = "black"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank()
  )
FIG4.B_bar_plot

# save
ggsave(paste("./Fig.RR3_", Sys.Date(), ".pdf", sep = ""), FIG4.A, width = 9, height = 6, units = "cm")
ggsave(paste("./Fig.RR3_", Sys.Date(), ".pdf", sep = ""), FIG4.B_bar_plot, width = 18, height = 15, units = "cm")

