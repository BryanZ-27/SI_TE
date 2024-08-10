###library
library(dplyr)
library(ggplot2)
library(colorspace)
library(patchwork)



###load data
load("./Data_for_Fig3.RData")



###main
#figure 3.b
FIG.3.B <- fig.3b_data %>%
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
  theme_classic()+theme(axis.text=element_text(size=9,face = "bold"),
                        axis.title = element_text(size=10,face = "bold"),
                        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
                        legend.background = element_blank(),legend.title =element_blank(),
                        legend.text = element_text(size = 8))

#figure 3.c
FIG.3.C <- fig.3c_data %>%
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
  theme_classic()+theme(axis.text=element_text(size=9,face = "bold"),
                        axis.title = element_text(size=10,face = "bold"),
                        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
                        legend.background = element_blank(),legend.title =element_blank(), 
                        legend.text = element_text(size = 8))

#figure 3.d1
FIG.3.D1 <- fig.3d_data %>%
  ggplot()+
  stat_bin(aes(TE_mean_bio), geom = "bar", bins = 30, fill = colorRampPalette(colors = c("#98c0e0", "#1c405d"))(30))+
  geom_vline(aes(xintercept = 0.000813,color="BY"),key_glyph="path",linewidth=1,show.legend = T,linetype=2)+
  geom_errorbar(aes(xmin=0.000813-0.0000865,xmax=0.000813+0.0000865,y = 32),size=1,width=2,color="#edb850")+
  geom_vline(aes(xintercept = 0.000993,color="RM"),key_glyph="path",linewidth=1,show.legend = T,linetype=2)+
  geom_errorbar(aes(xmin=0.000993-0.0000618,xmax=0.000993+0.0000618,y = 32),size=1,width=2,color="#624c7c",)+
  scale_color_manual(name="",values = c(BY="#edb850",RM="#624c7c"))+
  coord_cartesian(xlim = c(4E-4, 1.4E-3))+
  labs(x="",y="Count")+
  theme_classic()+theme(axis.text=element_text(size=9,face = "bold"),
                        axis.title = element_text(size=10,face = "bold"),
                        axis.text.x = element_blank(),axis.ticks.x = element_blank(),
                        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
                        legend.position = c(0.15,0.80),legend.text = element_text(size=8))

#figure 3.d2
FIG.3.D2 <- fig.3d_data %>%
  bind_rows(BY_TE_for_strain_data,RM_TE_for_strain_data) %>%
  ungroup() %>%
  ggplot()+
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
  scale_x_continuous(limits = c(0.0003, 0.0015), breaks = seq(0.0006, 0.0015, by = 0.0003), labels = c("6", "9", "12", "15")) + 
  coord_cartesian(xlim = c(4E-4, 1.4E-3))+
  labs(x="Translational error rate (e)",y="Strains")+
  theme_classic()+theme(axis.text=element_text(size=9,face = "bold"),
                        axis.title = element_text(size=10,face = "bold"),
                        axis.text.y = element_blank(),axis.ticks.y =element_blank(), 
                        panel.grid.major=element_blank(),panel.grid.minor=element_blank(), 
                        legend.position = "none")

FIG.3.D1/FIG.3.D2 -> Fig.3.D



###save
ggsave(paste("./Fig.3.B_", Sys.Date(), ".pdf", sep = ""), FIG.3.B, width = 9, height = 6, units = "cm")
ggsave(paste("./Fig.3.C_", Sys.Date(), ".pdf", sep = ""), FIG.3.C, width = 9, height = 6, units = "cm")
ggsave(paste("./Fig.3.D_", Sys.Date(), ".pdf", sep = ""), Fig.3.D, width = 9, height = 12, units = "cm")
