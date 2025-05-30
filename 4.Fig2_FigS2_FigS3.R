###figure 2b
# library
library(dplyr)
library(ggplot2)

# load data
# SI_rm_wt_gc_dpsg_wm_ol15_list 
# SI data
load("../data/4.SI_data.Rdata")

# MAIN#######
FIG.2.B_For_plot_ob <- SI_rm_wt_gc_dp_dpsg_wm_ol15_list[["SI_for_strain_clean_SI_biorep_ori_and_curve_rm_wt_gc_dp_data_output_list"]]$SI_for_strain_clean_SI_biorep_ori_and_curve_rm_na_sd_day2_min_0.2_0.5_shift_0.3_day2 %>%
  semi_join(SI_rm_wt_gc_dp_dpsg_wm_ol15_list$Strains_group_duplicate_for_strain_data_list_output$SI_for_strain_clean_SI_biorep_ori_and_curve_rm_na_sd_day2_min_0.2_0.5_shift_0.3_day2,by = c("strains"))

fig.2b_data <- FIG.2.B_For_plot_ob %>%
  filter(strains_group == "3H7H12-6A1A9") %>%
  left_join( FIG.2.B_For_plot_ob %>%
               filter(strains_group == "6A1-A12(6A1-A3rep)"),by = c("strains", "tube_type_YPD", "tube_type_SC", "glu_concentration"))

FIG.2.B <- 
  fig.2b_data %>%
  ggplot()+
  geom_point(aes(x=SI_mean_bio.x,y=SI_mean_bio.y))+
  geom_smooth(aes(x=SI_mean_bio.x,y=SI_mean_bio.y),method = "lm",se = F, color = "#4daf4a")+
  #stat_cor(aes(x=SI_mean_bio.x,y=SI_mean_bio.y),method = "pearson",size=6)+
  annotate("text", label = "paste(italic(R) == 0.92, ', ', italic(P) == 0.0029)", x = 5, y = 12.5, parse = T, size = 8/.pt) + 
  geom_errorbar(aes(x=SI_mean_bio.x,y=SI_mean_bio.y,
                    xmin=SI_mean_bio.x-sd_SI.x,xmax=SI_mean_bio.x+sd_SI.x))+
  geom_errorbar(aes(x=SI_mean_bio.x,y=SI_mean_bio.y,
                    ymin=SI_mean_bio.y-sd_SI.y,ymax=SI_mean_bio.y+sd_SI.y))+
  scale_x_continuous(limits = c(0, 15)) + 
  scale_y_continuous(limits = c(0, 15)) + 
  labs(x="Lifespan in batch 1",
       y="Lifespan in batch 2")+
  theme_classic()+theme(#plot.title = element_text(hjust = 0.5,size = 16, face = "bold"),
    axis.text=element_text(size=9,face = "bold"),
    axis.title = element_text(size=10,face = "bold"),
    axis.text.x = element_text(),
    panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
    # legend.position = c(.1,.85),legend.background = element_rect(colour = "black",size = 1.2,fill = "lightyellow"),legend.title = element_text(size=12,face = "bold"))
  )

# save
ggsave(paste("./Fig.2.B_", Sys.Date(), ".pdf", sep = ""), FIG.2.B, width = 9, height = 6, units = "cm")


###figure 2c
# load data
# SI_rm_wt_gc_dpsg_wm_ol15_list 
# SI data
load("../data/4.SI_data.Rdata")

# MAIN#######
#1.strain SI data with different glu_concentration####
Glu_0.5_data <- SI_rm_wt_gc_dp_dpsg_wm_ol15_list$Glu_concentration_for_biorep_data_output_list$SI_for_biorep_clean_SI_biorep_ori_and_curve_day2_min_0.2_0.5_shift_0.3_day2
Glu_2_data <- SI_rm_wt_gc_dp_dpsg_wm_ol15_list$SI_for_biorep_clean_SI_biorep_ori_and_curve_rm_wt_gc_dp_dpsg_wm_ol15_data_output_list$SI_for_biorep_clean_SI_biorep_ori_and_curve_day2_min_0.2_0.5_shift_0.3_day2 %>%
  filter(strains == "A0626")

#Plot#
FIG.2.C_For_plot_ob <- Glu_0.5_data %>% bind_rows(Glu_2_data) %>%
  mutate(glu_concentration=factor(glu_concentration,levels = c(0.5,2)))

fig.2c_color <- c("#41933e", "#73c371")

FIG.2.C <- 
  ggplot(FIG.2.C_For_plot_ob,aes(glu_concentration,SI,color=glu_concentration,fill=glu_concentration))+
  #plot
  geom_bar(stat="summary",fun=mean,position="dodge")+
  geom_jitter(width = 0.2, size = 1.5)+
  stat_summary(fun.data = 'mean_sd', geom = "errorbar", width = 0.5)+
  #t-test
  stat_compare_means(aes(label = "p.format"),method="t.test",
                     method.args = list(alternative = "greater"), comparisons = list(c("0.5", "2")), size = 8/.pt)+
  #scale
  scale_y_continuous(limits = c(0, 13), breaks = seq(0, 12, 3), labels = seq(0, 12, 3))+
  scale_color_manual(values = rep("black", 2))+
  scale_fill_manual(values = fig.2c_color)+
  #label
  labs(x="Glucose concentration (%)",y="Lifespan")+
  #theme
  theme_classic()+theme(axis.text=element_text(size=9),
                        axis.title = element_text(size=10,face = "bold"),
                        axis.text.x = element_text(color = "black"),
                        axis.text.y = element_text(color = "black"),
                        panel.grid.major=element_blank(),
                        panel.grid.minor=element_blank(),
                        legend.position = "none")

# save
ggsave(paste("./Fig.2.C_", Sys.Date(), ".pdf", sep = ""), FIG.2.C, width = 9, height = 6, units = "cm")


###figure 2d
# load data
# SI_rm_wt_gc_dpsg_wm_ol15_list 
# SI data
load("../data/4.SI_data.Rdata")

BY_SI_for_strain_in_two_strains_group_data <- SI_rm_wt_gc_dp_dpsg_wm_ol15_list$BY_for_strain_data_output_list$SI_for_strain_clean_SI_biorep_ori_and_curve_rm_na_sd_day2_min_0.2_0.5_shift_0.3_day2

BY_SI_for_PLOT <- BY_SI_for_strain_in_two_strains_group_data %>% ungroup() %>%
  group_by(strains) %>%
  dplyr::summarise(sd_SI=sd(SI_mean_bio),SI_mean_bio=mean(SI_mean_bio)) %>%
  mutate(strains=ifelse(strains %in% c("BY4716"),"BY",strains))

RM_SI_for_strain_in_two_strains_group_data <- SI_rm_wt_gc_dp_dpsg_wm_ol15_list$RM_for_strain_data_output_list$SI_for_strain_clean_SI_biorep_ori_and_curve_rm_na_sd_day2_min_0.2_0.5_shift_0.3_day2

RM_SI_for_PLOT <- RM_SI_for_strain_in_two_strains_group_data%>% ungroup() %>%
  group_by(strains) %>%
  dplyr::summarise(sd_SI=sd(SI_mean_bio),SI_mean_bio=mean(SI_mean_bio))%>%
  mutate(strains=ifelse(strains %in% c("RM11"),"RM",strains))

# MAIN#######
#1.get the SI phenotype data 
FIG2.D.E_For_plot_ob <- SI_rm_wt_gc_dp_dpsg_wm_ol15_list$SI_for_strain_clean_SI_biorep_ori_and_curve_rm_wt_gc_dp_dpsg_wm_ol15_data_output_list$SI_for_strain_clean_SI_biorep_ori_and_curve_rm_na_sd_day2_min_0.2_0.5_shift_0.3_day2

#PLOT#####
fig.2de_data <- FIG2.D.E_For_plot_ob %>%
  ungroup() %>% 
  arrange(SI_mean_bio) %>% 
  mutate(bin = c(rep(1:3, each = 26), rep(4:27, each = 27), rep(28:30, each = 26)))

FIG2.D <- 
  fig.2de_data %>% 
  ggplot()+
  stat_bin(aes(SI_mean_bio), geom = "bar", bins = 30, fill = colorRampPalette(colors = c("#aadaa8", "#295c27"))(30))+
  geom_vline(aes(xintercept = 8.33,color="BY"),key_glyph="path",linewidth=1,show.legend = T,linetype=2)+
  geom_errorbar(aes(xmin=8.33-0.538,xmax=8.33+0.538,y = 60),size=1,width=2,color="#edb850")+
  geom_vline(aes(xintercept = 11.6,color="RM"),key_glyph="path",linewidth=1,show.legend = T,linetype=2)+
  geom_errorbar(aes(xmin=11.6-0.880,xmax=11.6+0.880,y = 60),size=1,width=2,color="#624c7c",)+
  scale_color_manual(name="",values = c(BY="#edb850",RM="#624c7c"))+
  coord_cartesian(xlim = c(0, 15))+
  labs(x="",y="Count")+
  theme_classic()+theme(#plot.title = element_text(hjust = 0.5,size = 16, face = "bold"),
    axis.text=element_text(size=9,face = "bold", color = "black"),
    axis.title = element_text(size=10,face = "bold"),
    axis.text.x = element_blank(),axis.ticks.x = element_blank(),
    panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
    legend.position = c(0.15,0.80),legend.text = element_text(size=8))
###
##@##2.Figure.2.E ####
FIG2.E <- 
  FIG2.D.E_For_plot_ob %>%
  bind_rows(BY_SI_for_PLOT,RM_SI_for_PLOT) %>%
  ggplot()+
  geom_point(aes(x=SI_mean_bio,y=reorder(strains,SI_mean_bio), 
                 color = SI_mean_bio), size = 0.01)+
  geom_errorbar(aes(xmin=SI_mean_bio-sd_SI,xmax=SI_mean_bio+sd_SI,y=reorder(strains,SI_mean_bio), 
                    color = SI_mean_bio), width = 0, linewidth = 0.1)+
  geom_point(mapping = aes(x=SI_mean_bio,y=reorder(strains,SI_mean_bio)),color=c(BY="#edb850",RM="#624c7c"),
             data = BY_SI_for_PLOT %>% bind_rows(RM_SI_for_PLOT),show.legend = F, size = 2)+
  geom_errorbar(aes(xmin=SI_mean_bio-sd_SI,xmax=SI_mean_bio+sd_SI,y=reorder(strains,SI_mean_bio)),
                color=c(BY="#edb850",RM="#624c7c"),size = 1.1,width=30,
                data = BY_SI_for_PLOT %>% bind_rows(RM_SI_for_PLOT),show.legend = F)+
  geom_vline(aes(xintercept = SI_mean_bio),
             color=c(BY="#edb850",RM="#624c7c"),linewidth = 1,linetype=2,
             data = BY_SI_for_PLOT %>% bind_rows(RM_SI_for_PLOT),show.legend = F)+
  coord_cartesian(xlim = c(0, 15))+
  scale_color_gradient(low = "#aadaa8", high = "#295c27")+
  labs(x="Lifespan",y="Strains")+
  theme_classic()+theme(#plot.title = element_text(hjust = 0.5,size = 16, face = "bold"),
    axis.text=element_text(size=9,face = "bold", color = "black"),
    axis.title = element_text(size=10,face = "bold"),
    axis.text.y = element_blank(),axis.ticks.y =element_blank(), 
    panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
    legend.position = "none")

library(patchwork)
FIG2.D/FIG2.E -> Fig2.D_E

# save
ggsave(paste("./Fig.2.D-E_", Sys.Date(), ".pdf", sep = ""), Fig2.D_E, width = 9, height = 12, units = "cm")


###figure S2
# library
library(reshape2)
library(ggbreak)

# load data
load("../data/4.FigS2_data.Rdata")

# MAIN#######
#PLOT#####
fig.s2a_summ <- Alignment_score_in_genotype_verMT_data %>% 
  group_by(strain) %>% 
  dplyr::summarise(mean_score = mean(Alignment_score, na.rm = T), 
                   sd_score = sd(Alignment_score, na.rm = T), 
                   rate = sum(Alignment_score >= 0.90, na.rm = T)/n()) %>% 
  ungroup() %>% arrange(mean_score) %>% 
  mutate(no = seq(19, 498.9997, 53.3333))
FIG.S2.B_For_plot_ob <- Summarize_Alignment_score_data %>% melt %>%.[-1,]
FigureS2.A  <- ggplot()+
  geom_point(aes(mean_score, no, color = mean_score), fig.s2a_summ, size = 2)+
  geom_errorbar(aes(xmin = mean_score - sd_score, 
                    xmax = ifelse(mean_score + sd_score > 1, 1, mean_score + sd_score), 
                    y = no, color = mean_score), fig.s2a_summ, width = 10, linewidth = 0.5)+
  geom_freqpoly(aes(value), FIG.S2.B_For_plot_ob, color = "#aadaa8", binwidth = 0.02) + 
  scale_color_gradient(low = "#4daf4a", high = "#295c27")+
  scale_x_continuous(name = "The fraction of reads consistent with 
  the genotype data among all reads 
                     aligned to a segregating site", 
                     limits = c(0.4, 1), breaks = seq(0.4, 1, 0.1))+
  scale_y_continuous(name = "Count", 
                     sec.axis = sec_axis(~., name = "Strains"))+
  theme_classic()+theme(axis.text=element_text(size = 9, face = "bold", color = "black"),
    axis.title = element_text(size = 10, face = "bold", color = "black"),
    axis.text.y.right = element_blank(),
    panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
    legend.position = "none")

#save
ggsave(paste("./Fig.S2A_", Sys.Date(), ".pdf", sep = ""), FigureS2.A, width = 8, height = 9, units = "cm")


###figure s3
# library
library(ggpubr)

# load data
# SI data 
load("../data/4.SI_data.Rdata")

# MAIN#######
#PLOT######
SI_data_with_different_argument <- SI_rm_wt_gc_dp_dpsg_wm_ol15_list$SI_for_strain_clean_SI_biorep_ori_and_curve_rm_wt_gc_dp_dpsg_wm_ol15_data_output_list %>% bind_rows(.id = "Argument") %>%
  mutate(Argument_1 = gsub(pattern = ".*sd_(.*)","\\1",x = Argument)) %>%
  dcast(strains+strains_group+tube_type_YPD+tube_type_SC+glu_concentration ~ Argument_1,value.var = ("SI_mean_bio"))

SI_cor_matrix_with_different_argument <- SI_data_with_different_argument[,-c(1:5)] %>% 
  corrr::correlate() %>%corrr::stretch()

#B###
FigureS3.B_For_plot_ob <- SI_cor_matrix_with_different_argument %>%
  filter(x%in% c("day2_min_0.2_0.5_shift_0.3_day2")) %>%
  mutate(r = ifelse(is.na(r),1,r)) %>%
  mutate(Calculate_Method = gsub("(.*)_(\\d\\.\\d+_\\d\\.\\d+)_shift_(\\d\\.\\d)_day2","\\1",y),
         OD_range = gsub("(.*)_(\\d\\.\\d+_\\d\\.\\d+)_shift_(\\d\\.\\d)_day2","\\2",y),
         OD_time_shift = as.numeric(gsub("(.*)_(\\d\\.\\d+_\\d\\.\\d+)_shift_(\\d\\.\\d)_day2","\\3",y))) %>% 
  mutate(`Correlation with the star` = r) %>% 
  mutate(Calculate_Method = gsub("day2_min", "1", Calculate_Method), 
         Calculate_Method = gsub("median", "2", Calculate_Method), 
         Calculate_Method = gsub("day2_mean", "3", Calculate_Method))

FigureS3.B_For_plot_ob %>%
  ggplot(aes(x=OD_range, y=OD_time_shift, fill=`Correlation with the star`,group=interaction(Calculate_Method)))+
  geom_tile()+
  scale_fill_gradient(low = "#FFFFDF", high = "#D73027")+
  facet_wrap(~Calculate_Method,labeller =as_labeller(c(`1`="Day2 min",`2`="Median",`3`="Day2 mean")))+
  labs(x="OD range",y="Time shift")+
  scale_x_discrete(labels=c("[0.2, 0.5]","[0.25, 0.75]","[0.3, 0.8]"))+
  theme_classic()+
  theme(axis.text= element_text(size=9,face = "bold", color = "black"),
    axis.title = element_text(size=10,face = "bold", color = "black"),
    axis.text.x = element_blank(),
    panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
    legend.background = element_blank(),
    legend.title = element_blank(),
    strip.text = element_text(size=8,color = "black",face = "bold")) -> Fig.S3.B

# save
ggsave(paste("./Fig.S3.B_", Sys.Date(), ".pdf", sep = ""), Fig.S3.B, width = 12, height = 8-0.5, units = "cm")

