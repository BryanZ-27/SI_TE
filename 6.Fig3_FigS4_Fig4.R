###figure 3b & figure 3c
# library
library(dplyr)
library(ggplot2)
library(ggpmisc)
library(colorspace)

# load data
load("./5.TE_techrep_data.Rdata")

# MAIN#######
FIG.3.B_For_plot_ob <- TE_techrep_clean_list$TE_11_techrep_clean_list$TE_11_techrep_clean_data

FIG.3.B <- 
  FIG.3.B_For_plot_ob %>%
  ggplot() +
  geom_hex(aes(x=`F`,y=`R`),bins = 11)+
  geom_smooth(aes(x=`F`,y=`R`), method = "lm", formula = y ~ x + 0, se = F, color = "#027ee1")+
  scale_fill_gradientn(colors = sequential_hcl(n = 10, palette = "Grays", rev = T)[2:10]) + 
  annotate("text", label = "paste(italic(y) == 7.16 * italic(x), ', ', italic(R)^2 == 0.98)", 
           x = 2500000, y = 50000000, parse = T, size = 8/.pt) + 
  labs(x="Firefly luminescent signal\n(wild type) (e+6)",y="Renilla luminescent signal\n(wild type) (e+7)")+
  scale_x_continuous(breaks = c(1e+6,2e+6,3e+6,4e+6,5e+6,6e+6),
                     labels = c(1,2,3,4,5,6))+
  scale_y_continuous(breaks = c(1e+7,2e+7,3e+7,4e+7,5e+7),
                     labels = c(1,2,3,4,5))+
  theme_classic()+theme(#plot.title = element_text(hjust = 0.5,size = 16, face = "bold"),
    axis.text=element_text(size=9,face = "bold"),
    axis.title = element_text(size=10,face = "bold"),
    panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
    legend.background = element_blank(),legend.title =element_blank(),
    legend.text = element_text(size = 8))

FIG.3.C_For_plot_ob <- TE_techrep_clean_list$TE_13_techrep_clean_list$TE_13_techrep_clean_data 

FIG.3.C <- 
  FIG.3.C_For_plot_ob %>%
  ggplot() +
  geom_hex(aes(x=`F`,y=`R`),bins = 11)+
  geom_smooth(aes(x=`F`,y=`R`), method = "lm", formula = y ~ x + 0, se = F, color = "#027ee1")+
  scale_fill_gradientn(colors = sequential_hcl(n = 10, palette = "Grays", rev = T)[2:10]) + 
  annotate("text", label = "paste(italic(y) == 7.77 %*% 10^3 * italic(x), ', ', italic(R)^2 == 0.95)", 
           x = 4000, y = 60000000, parse = T, size = 8 / .pt) + 
  labs(x="Firefly luminescent signal\n(mutant type) (e+3)",y="Renilla luminescent signal\n(mutant type) (e+7)")+
  scale_x_continuous(breaks = c(0,2e+3,4e+3,6e+3,8e+3,1e+4),
                     labels = c(0,2,4,6,8,10))+
  scale_y_continuous(breaks = c(1e+7,2e+7,3e+7,4e+7,5e+7,6e+7),
                     labels = c(1,2,3,4,5,6))+
  theme_classic()+theme(#plot.title = element_text(hjust = 0.5,size = 16, face = "bold"),
    axis.text=element_text(size=9,face = "bold"),
    axis.title = element_text(size=10,face = "bold"),
    panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
    legend.background = element_blank(),legend.title =element_blank(), 
    legend.text = element_text(size = 8))

# save
ggsave(paste("./Fig.3.B_", Sys.Date(), ".pdf", sep = ""), FIG.3.B, width = 9, height = 6, units = "cm")
ggsave(paste("./Fig.3.C_", Sys.Date(), ".pdf", sep = ""), FIG.3.C, width = 9, height = 6, units = "cm")


###figure 3d
# library
library(stringr)

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

# load data
#TE data ######
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

#BY RM TE data ######
BY_11_biorep_data <- TE_env$TE_11_biorep_clean_list$TE_11_for_biorep_clean_bio_data %>%
  filter(strain == "BY4716")
BY_13_biorep_data <- TE_env$TE_13_biorep_clean_list$TE_13_for_biorep_clean_bio_data %>%
  filter(strain == "BY4716")
BY_TE_biorep_data <- BY_13_biorep_data %>% 
  full_join(BY_11_biorep_data,by = c("strain","panel", "seed_volume", "Sample_volume", "lysate_type", 
                                     "Time_to_shape_lysate", "Reagent_volume")) %>%
  mutate(strain = ifelse(strain == "BY4716","BY",strain)) %>%
  mutate(TE = F_R_ratio.x/F_R_ratio.y)
BY_TE_for_strain_data <- BY_TE_biorep_data %>%
  group_by(strain) %>%
  dplyr::summarise(sd_TE=sd(TE,na.rm = T),TE_mean_bio=mean(TE,na.rm=T) )

RM_11_biorep_data <-  TE_env$TE_11_biorep_clean_list$TE_11_for_biorep_clean_bio_data %>%
  filter(strain == "RM11")
RM_11_for_strain_data <- RM_11_biorep_data %>% 
  group_by(strain, mutant_type, panel, seed_volume,Sample_volume, lysate_type, Time_to_shape_lysate) %>%
  dplyr::summarise(sd_F_R_ratio=sd(F_R_ratio),F_R_ratio_mean_bio=mean(F_R_ratio ) )
RM_13_biorep_data <-  TE_env$TE_13_biorep_clean_list$TE_13_for_biorep_clean_bio_data %>%
  filter(strain == "RM11")
RM_TE_for_biorep_data <- RM_13_biorep_data %>% 
  left_join(RM_11_for_strain_data, c("strain",  "panel", "seed_volume", "Sample_volume", "lysate_type", 
                                     "Time_to_shape_lysate")) %>%
  mutate(strain = ifelse(strain == "RM11","RM",strain))  %>%
  mutate(TE = F_R_ratio/F_R_ratio_mean_bio)

RM_TE_for_strain_data <- RM_TE_for_biorep_data %>% group_by(strain) %>% 
  dplyr::summarise(sd_TE=sd(TE,na.rm = T),TE_mean_bio=mean(TE,na.rm=T) )

# MAIN#######
TE_for_strians_260_data <- TE_271_and_11_13_RM11point %>%
  filter(!is.na(TE_mean_bio))

FIG.3.D.E_For_plot_ob <- TE_for_strians_260_data

FIG.3.D <- 
  FIG.3.D.E_For_plot_ob %>%
  ggplot()+
  stat_bin(aes(TE_mean_bio), geom = "bar", bins = 30, fill = colorRampPalette(colors = c("#98c0e0", "#1c405d"))(30))+
  # scale_x_continuous(breaks = seq(0, 20, by = 2.5))+
  geom_vline(aes(xintercept = 0.000813,color="BY"),key_glyph="path",linewidth=1,show.legend = T,linetype=2)+
  geom_errorbar(aes(xmin=0.000813-0.0000865,xmax=0.000813+0.0000865,y = 32),size=1,width=2,color="#edb850")+
  geom_vline(aes(xintercept = 0.000993,color="RM"),key_glyph="path",linewidth=1,show.legend = T,linetype=2)+
  geom_errorbar(aes(xmin=0.000993-0.0000618,xmax=0.000993+0.0000618,y = 32),size=1,width=2,color="#624c7c",)+
  scale_color_manual(name="",values = c(BY="#edb850",RM="#624c7c"))+
  coord_cartesian(xlim = c(4E-4, 1.4E-3))+
  labs(x="",y="Count")+
  theme_classic()+theme(#plot.title = element_text(hjust = 0.5,size = 16, face = "bold"),
    axis.text=element_text(size=9,face = "bold"),
    axis.title = element_text(size=10,face = "bold"),
    axis.text.x = element_blank(),axis.ticks.x = element_blank(),
    panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
    legend.position = c(0.15,0.80),legend.text = element_text(size=8))

FIG.3.E <- 
  FIG.3.D.E_For_plot_ob %>%
  bind_rows(BY_TE_for_strain_data,RM_TE_for_strain_data) %>%
  ungroup() %>% ggplot()+
  geom_point(aes(x=TE_mean_bio,y=reorder(strain,TE_mean_bio), 
                 color = TE_mean_bio), size = 0.01)+
  geom_errorbar(aes(xmin=TE_mean_bio-sd_TE,xmax=TE_mean_bio+sd_TE,y=reorder(strain,TE_mean_bio), 
                    color = TE_mean_bio), width = 0, linewidth = 0.1)+
  geom_point(mapping = aes(x=TE_mean_bio,y=reorder(strain,TE_mean_bio)),color=c(BY="#edb850",RM="#624c7c"),
             data = BY_TE_for_strain_data %>% bind_rows(RM_TE_for_strain_data),show.legend = F, size = 2)+
  geom_errorbar(aes(xmin=TE_mean_bio-sd_TE,xmax=TE_mean_bio+sd_TE,y=reorder(strain,TE_mean_bio)),
                color=c(BY="#edb850",RM="#624c7c"),size = 1.1,width=15,
                data = BY_TE_for_strain_data %>% bind_rows(RM_TE_for_strain_data),show.legend = F)+
  geom_vline(aes(xintercept = TE_mean_bio),
             color=c(BY="#edb850",RM="#624c7c"),linewidth = 1,linetype=2,
             data = BY_TE_for_strain_data %>% bind_rows(RM_TE_for_strain_data),show.legend = F)+
  scale_color_gradient(low = "#98c0e0", high = "#1c405d")+
  scale_x_continuous(limits = c(0.0003, 0.0015), breaks = seq(0.0006, 0.0015, by = 0.0003), 
                     labels = c("6", "9", "12", "15")) + 
  coord_cartesian(xlim = c(4E-4, 1.4E-3))+
  labs(x="Translational error rate (e)",y="Strains")+
  theme_classic()+theme(#plot.title = element_text(hjust = 0.5,size = 16, face = "bold"),
    axis.text=element_text(size=9,face = "bold"),
    axis.title = element_text(size=10,face = "bold"),
    axis.text.y = element_blank(),axis.ticks.y =element_blank(), 
    panel.grid.major=element_blank(),panel.grid.minor=element_blank(), 
    legend.position = "none")

library(patchwork)
FIG.3.D/FIG.3.E -> Fig.3.D_E

# save
ggsave(paste("./Fig.3.D_E_", Sys.Date(), ".pdf", sep = ""), Fig.3.D_E, width = 9, height = 12, units = "cm")


###figure S4
# library
library(ggbreak)

# load data
load("./6.FigS4_data.Rdata")

# MAIN#######
#PLOT#####
FIG.S4.A_For_plot_ob <- Alignment_score_in_genotype_verMT_data %>%
  filter(!strain %in% c("BY4716","RM11","A0141","A1035","A1029",
                        "A1024","A0190","A0645","A0863","A0306",
                        "A1052","A1085","A1044",
                        "A0110","A0431","A1084","A0479","A0456")) %>% 
  filter(CHR != "chrIV")
fig.s4a_summ <- FIG.S4.A_For_plot_ob %>% 
  group_by(strain) %>% 
  dplyr::summarise(mean_score = mean(Alignment_score, na.rm = T), 
                   sd_score = sd(Alignment_score, na.rm = T), 
                   rate = sum(Alignment_score >= 0.9, na.rm = T)/n()) %>% 
  ungroup() %>% arrange(mean_score) %>% 
  mutate(no = seq(21.72, 999.9889, 3.7761))

FIG.S4.B_For_plot_ob <- True_plus_fake_alignment_score_data_frame_wihout_chr4 %>%  
  filter(!strain %in% c("BY4716","RM11","A0141","A1035","A1029",
                        "A1024","A0190","A0645","A0863","A0306",
                        "A1052","A1085","A1044",
                        "A0110","A0431","A1084","A0479","A0456")) %>% 
  dplyr::summarise_at(vars(Alignment_score:Alignment_score_1000), mean, na.rm = TRUE) %>%
  melt() %>%.[-1,]

FigureS4.A  <- 
  ggplot()+
  geom_point(aes(mean_score, no, color = mean_score), fig.s4a_summ, size = 0.1)+
  geom_errorbar(aes(xmin = mean_score - sd_score, 
                    xmax = ifelse(mean_score + sd_score > 1, 1, mean_score + sd_score), 
                    y = no, color = mean_score), fig.s4a_summ, width = 0.1, linewidth = 0.1)+
  #geom_point(aes(rate, no), fig.s4a_summ, 
  #           size = 0.4, shape = 17, color = sequential_hcl(n = 519, palette = "Grays", rev = T)[260:519])+
  #geom_density(aes(value), FIG.S4.B_For_plot_ob, color = "#8dc1ec")+
  geom_freqpoly(aes(value), FIG.S4.B_For_plot_ob, color = "#8dc1ec") + 
  scale_color_gradient(low = "#2075bc", high = "#124168")+
  scale_x_continuous(name = "The fraction of reads consistent with 
  the genotype data among all reads 
                     aligned to a segregating site", 
                     limits = c(0.4, 1), breaks = seq(0.4, 1, 0.1), 
                     #                   sec.axis = sec_axis(~., breaks = seq(0.4, 1, 0.1), 
                     #                                       labels = paste(seq(40, 100, 10), "%", sep = ""), 
                     #                                       name = "The proportion of alignment score >= 0.95")
  )+
  scale_y_continuous(name = "Count", 
                     sec.axis = sec_axis(~., name = "Strains"))+
  #labs(x="",y="Strain")+
  #scale_y_continuous(limits = c(0, 10000), breaks = seq(0, 10000, 2000), labels = seq(0, 1, 0.2))+
  #scale_y_break(c(300, 9000))+
  theme_classic()+theme(#plot.title = element_text(hjust = 0.5,size = 16, face = "bold"),
    axis.text=element_text(size = 9, face = "bold", color = "black"),
    axis.title = element_text(size = 10, face = "bold", color = "black"),
    axis.text.y.right = element_blank(),
    #axis.text.x =element_blank(),axis.title.x=element_blank(), 
    panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
    legend.position = "none")

#save
ggsave(paste("./Fig.S4A_", Sys.Date(), ".pdf", sep = ""), FigureS4.A, width = 8, height = 9, units = "cm")


###figure 4
# library
library(tidyr)
library(purrr)
library(doParallel)
library(generics)

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
SI_TE_cor_diff_SI_cutoff_output_list <-  mclapply(mc.cores = 16, X = c(seq(from=0,to=0.9,by=0.1)) ,FUN = function(cutoff){
  output <- TE_SI_data %>% 
    dplyr::mutate(N=n()) %>%
    arrange(SI_mean_bio) %>%
    filter(row_number() >=N*cutoff) %>%
    nest() %>%
    mutate(res = purrr::map(data, .f = ~ tidy(cor.test(~TE_mean_bio +SI_mean_bio, data = .,method="spearman")))) %>%
    unnest(res) %>%
    mutate(observation_number = purrr::map(data,.f = function(data){nrow(data)})) %>%
    mutate(cut_off = cutoff) %>%
    unnest(observation_number)
  return(output)
})

#@## 2.2 point list correspond to cor list###
SI_TE_point_plot_diff_SI_cutoff_output_list <- mclapply(X = c(seq(from=0,to=0.9,by=0.1)) ,FUN = function(cutoff){
  output <- TE_SI_data %>% 
    dplyr::mutate(N=n()) %>%
    arrange(SI_mean_bio) %>%
    filter(row_number() >=N*cutoff) %>%
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
  geom_bar(aes(x=order,y=-bootstrap_mean*20+shift, fill = as.factor(cut_off)), 
           stat = "identity",position =position_dodge(width = 0.15),width = 0.12)+
  geom_bar(aes(x=order,y=log10(adjust_p)*4-shift, fill = as.factor(cut_off)), 
           stat = "identity",position =position_dodge(width = 0.15),width = 0.02)+
  geom_errorbar(aes(x = order, ymin = -(bootstrap_mean - bootstrap_SD)*20+shift, 
                    ymax = ifelse(-(bootstrap_mean + bootstrap_SD)*20+shift > 0.55, 
                    -(bootstrap_mean + bootstrap_SD)*20+shift, 0.55)), 
                stat = "identity",position = position_dodge(width = 0.15),width = 0.08, color = "black") + 
  geom_point(aes(x=order,y=log10(adjust_p)*4-shift, fill = as.factor(cut_off)), 
             position = position_dodge(width = 0.15), size=18, color = fig4_color, shape=21)+
  geom_abline(slope = 0,intercept=log10(0.05)*4-shift,color="red",linewidth=1.2,linetype=2)+
  annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymax = shift, ymin = -shift,fill = "white") +
  annotate(geom = "text", y = 0, x = new_fig4b_data$order %>% unique(), 
           label = new_fig4b_data$cut_off %>% unique(), size = 9/.pt, color = "black", fontface = "bold") +
  scale_y_continuous(limits = c(-15, 15), 
                     breaks=c(-12-shift, -8-shift, -4-shift, 0-shift, 0+shift, 4+shift, 8+shift ,12+shift), 
                     labels = c(expression(10^-3), expression(10^-2), expression(10^-1), 0, 
                                0, "-0.2", "-0.4", "-0.6"))+
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

# save
ggsave(paste("./Fig.4A_", Sys.Date(), ".pdf", sep = ""), FIG4.A, width = 9, height = 6, units = "cm")
ggsave(paste("./Fig.4B_", Sys.Date(), ".pdf", sep = ""), FIG4.B_bar_plot, width = 18, height = 15, units = "cm")

