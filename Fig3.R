###library
library(dplyr)
library(ggplot2)
library(ggpubr)
library(patchwork)


###load data
load("./Data_for_Fig3.RData")


###MAIN
##Fig3.b
fig.3b_color <- c("#41933e", "#73c371")

FIG.3b <- ggplot(FIG.3b_For_plot_ob,aes(glu_concentration,SI,color=glu_concentration,fill=glu_concentration))+
  geom_bar(stat="summary",fun=mean,position="dodge")+
  geom_jitter(width = 0.2, size = 1.5)+
  stat_summary(fun.data = 'mean_sd', geom = "errorbar", width = 0.5)+
  stat_compare_means(aes(label = "p.format"),method="t.test",
                     method.args = list(alternative = "greater"), comparisons = list(c("0.5", "2")), size = 8/.pt)+
  scale_y_continuous(limits = c(0, 13), breaks = seq(0, 12, 3), labels = seq(0, 12, 3))+
  scale_color_manual(values = rep("black", 2))+
  scale_fill_manual(values = fig.3b_color)+
  labs(x="Glucose concentration (%)",y="Lifespan")+
  theme_classic()+theme(axis.text=element_text(size=9),
                        axis.title = element_text(size=10,face = "bold"),
                        axis.text.x = element_text(color = "black"),
                        axis.text.y = element_text(color = "black"),
                        panel.grid.major=element_blank(),
                        panel.grid.minor=element_blank(),
                        legend.position = "none")

##Fig3.c
FIG.3.c <- fig.3c_data %>%
  ggplot()+
  geom_point(aes(x=SI_mean_bio.x,y=SI_mean_bio.y))+
  geom_smooth(aes(x=SI_mean_bio.x,y=SI_mean_bio.y),method = "lm",se = F, color = "#4daf4a")+
  annotate("text", label = "paste(italic(R) == 0.92, ', ', italic(P) == 0.0029)", x = 5, y = 12.5, parse = T, size = 8/.pt) + 
  geom_errorbar(aes(x=SI_mean_bio.x,y=SI_mean_bio.y,
                    xmin=SI_mean_bio.x-sd_SI.x,xmax=SI_mean_bio.x+sd_SI.x))+
  geom_errorbar(aes(x=SI_mean_bio.x,y=SI_mean_bio.y,
                    ymin=SI_mean_bio.y-sd_SI.y,ymax=SI_mean_bio.y+sd_SI.y))+
  labs(x="Lifespan in batch 1",
       y="Lifespan in batch 2")+
  theme_classic()+theme(axis.text=element_text(size=9,face = "bold"),
                        axis.title = element_text(size=10,face = "bold"),
                        axis.text.x = element_text(),
                        panel.grid.major=element_blank(),panel.grid.minor=element_blank())

##Fig3.d
FIG3.D1 <- fig.3d1_data %>% 
  ggplot()+
  stat_bin(aes(SI_mean_bio), geom = "bar", bins = 30, fill = colorRampPalette(colors = c("#aadaa8", "#295c27"))(30))+
  geom_vline(aes(xintercept = 8.33,color="BY"),key_glyph="path",linewidth=1,show.legend = T,linetype=2)+
  geom_errorbar(aes(xmin=8.33-0.538,xmax=8.33+0.538,y = 60),size=1,width=2,color="#edb850")+
  geom_vline(aes(xintercept = 11.6,color="RM"),key_glyph="path",linewidth=1,show.legend = T,linetype=2)+
  geom_errorbar(aes(xmin=11.6-0.880,xmax=11.6+0.880,y = 60),size=1,width=2,color="#624c7c",)+
  scale_color_manual(name="",values = c(BY="#edb850",RM="#624c7c"))+
  coord_cartesian(xlim = c(0, 15))+
  labs(x="",y="Count")+
  theme_classic()+theme(axis.text=element_text(size=9,face = "bold", color = "black"),
    axis.title = element_text(size=10,face = "bold"),
    axis.text.x = element_blank(),axis.ticks.x = element_blank(),
    panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
    legend.position = c(0.15,0.80),legend.text = element_text(size=8))

BY_SI_for_PLOT <- fig.3d2_data %>% filter(strains == "BY")
RM_SI_for_PLOT <- fig.3d2_data %>% filter(strains == "RM")
FIG3.D2 <- fig.3d2_data %>%
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
  theme_classic()+theme(axis.text=element_text(size=9,face = "bold", color = "black"),
    axis.title = element_text(size=10,face = "bold"),
    axis.text.y = element_blank(),axis.ticks.y =element_blank(), 
    panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
    legend.position = "none")

FIG3.D1/FIG3.D2 -> Fig3.D

###save
ggsave(paste("./Fig3.b_", Sys.Date(), ".pdf", sep = ""), FIG.3b, width = 9, height = 6, units = "cm")
ggsave(paste("./Fig3.c_", Sys.Date(), ".pdf", sep = ""), FIG.3.c, width = 9, height = 6, units = "cm")
ggsave(paste("./Fig3.d_", Sys.Date(), ".pdf", sep = ""), Fig3.D, width = 9, height = 12, units = "cm")