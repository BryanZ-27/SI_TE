###library
library(dplyr)
library(ggplot2)



###load data
load("./Data_for_FigS5.RData")



###main
#figure S5A
FIGS5.A <- fig.s5a_data %>% 
  ggplot()+
  geom_rect(data = chr_rect[seq(1,15,by=2),],mapping = aes(xmin = x1,xmax=x2, fill=c), ymin = -Inf,ymax = Inf,show.legend = F)+
  geom_point(aes(x=Marker_id,y=value,color=variable), size = 0.1)+
  geom_hline(yintercept = c(3.30),color = c("#4daf4a"),linetype = 2)+
  geom_hline(yintercept = c(3.51),color = c("#2075bc"),linetype = 2)+
  labs(x="Genome position",y="LOD score")+
  scale_x_continuous(breaks = chr_marker_id_mean/2)+
  scale_fill_manual(values=c("#BEBEBE"))+
  scale_color_manual(name="",values=c("#4daf4a","#2075bc"),label=c("Lifespan","Translation error rate"))+
  theme_classic()+
  theme(axis.text=element_text(size=9,face = "bold",color = "black"),
        axis.title = element_text(size=10,face = "bold"),
        axis.text.x = element_text(angle = 20),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        legend.text = element_text(size=12,face = "bold"),legend.position = "none")

#figure S5B
FigureS5.B <- fig.s5b_data %>% 
  ggplot()+
  geom_point(aes(x=TE_mean_bio,y=SI_mean_bio), size = 0.2)+
  stat_cor(aes(x=TE_mean_bio,y=SI_mean_bio),method = "spearman",cor.coef.name = "rho", size = 8/.pt)+
  geom_smooth(aes(x=TE_mean_bio,y=SI_mean_bio), formula = "y ~ x",method = "lm",se = F, color = "#4daf4a")+
  labs(x="Translational error rate (e) (x10-4)",y="Lifespan")+
  scale_x_continuous(limit = c(6.6e-4, 12.6e-4), breaks = seq(7e-4, 11e-4, 2e-4), labels = seq(7, 11, 2))+
  scale_y_continuous(limits = c(1.9, 14), breaks = seq(3, 12, 3), labels = seq(3, 12, 3)) + 
  theme_classic()+theme(axis.text=element_text(size=9,face = "bold", color = "black"),
                        axis.title = element_text(size=10,face = "bold", color = "black"),
                        axis.title.x = element_blank(),axis.title.y = element_blank(),
                        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
                        legend.position = "none")

#figure S5C
FIGS5.C <- fig.s5c_data %>%
  na.omit() %>%
  ggplot(group=!!sym("7905399_chrXII_660371_C_T"))+
  geom_point(aes(x=!!sym("7905399_chrXII_660371_C_T"),y=pheno,color = !!sym("7905399_chrXII_660371_C_T")),position = position_jitter(w = 0.2, h = 0), 
             size = 0.3)+
  geom_errorbar(aes(ymin=fig.s5c_summ$m[1]-fig.s5c_summ$se[1],ymax=fig.s5c_summ$m[1]+fig.s5c_summ$se[1],x=1),
                size=0.3,width=0.2)+
  geom_errorbar(aes(ymin=fig.s5c_summ$m[2]-fig.s5c_summ$se[2],ymax=fig.s5c_summ$m[2]+fig.s5c_summ$se[2],x=2),
                size=0.3,width=0.2)+
  geom_segment(aes(y = fig.s5c_summ$m[1], yend = fig.s5c_summ$m[1], x = 0.75, xend = 1.25), size = 0.3)+
  geom_segment(aes(y = fig.s5c_summ$m[2], yend = fig.s5c_summ$m[2], x = 1.75, xend = 2.25), size = 0.3)+
  scale_x_discrete(labels= c("BY","RM"))+
  scale_color_manual(values = c("#e5d851", "#8c07bb"))+
  scale_y_continuous(limits=c(1.9, 13.9), breaks = seq(3, 12, 3), labels = seq(3, 12, 3))+
  labs(title = "7905399_chrXII_660371_C_T",x="Genotype",y="Lifespan")+
  theme_classic()+theme(axis.text=element_text(size=9,face = "bold",color = "black"),
                        axis.title = element_text(size=10,face = "bold",color = "black"),
                        axis.text.x = element_text(vjust = 1),
                        axis.title.x = element_blank(),plot.title = element_blank(),
                        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
                        legend.position = "none")

#figure S5D
FIGS5.D <- fig.s5d_data %>%
  na.omit() %>%
  ggplot(group=!!sym("9709121_chrXIV_461485_A_G"))+
  geom_point(aes(x=!!sym("9709121_chrXIV_461485_A_G"),y=pheno,color = !!sym("9709121_chrXIV_461485_A_G")),position = position_jitter(w = 0.2, h = 0), 
             size = 0.3)+
  geom_errorbar(aes(ymin=fig.s5d_summ$m[1]-fig.s5d_summ$se[1],ymax=fig.s5d_summ$m[1]+fig.s5d_summ$se[1],x=1),
                size=0.3,width=0.2)+
  geom_errorbar(aes(ymin=fig.s5d_summ$m[2]-fig.s5d_summ$se[2],ymax=fig.s5d_summ$m[2]+fig.s5d_summ$se[2],x=2),
                size=0.3,width=0.2)+
  geom_segment(aes(y = fig.s5d_summ$m[1], yend = fig.s5d_summ$m[1], x = 0.75, xend = 1.25), size = 0.3)+
  geom_segment(aes(y = fig.s5d_summ$m[2], yend = fig.s5d_summ$m[2], x = 1.75, xend = 2.25), size = 0.3)+
  scale_x_discrete(labels= c("BY","RM"))+
  scale_color_manual(values = c("#e5d851", "#8c07bb"))+
  scale_y_continuous(limits=c(1.9, 13.9), breaks = seq(3, 12, 3), labels = seq(3, 12, 3))+
  labs(title = "9709121_chrXIV_461485_A_G",x="Genotype",y="Lifespan")+
  theme_classic()+theme(axis.text=element_text(size=9,face = "bold",color = "black"),
                        axis.title = element_text(size=10,face = "bold",color = "black"),
                        axis.text.x = element_text(vjust = 1),
                        axis.title.x = element_blank(),plot.title = element_blank(),
                        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
                        legend.position = "none")

#figure S5E
FigureS5.E <- fig.s5e_data %>%
  ggplot()+
  geom_point(aes(x=TE_mean_bio,y=SI_mean_bio), size = 0.2)+
  stat_cor(aes(x=TE_mean_bio,y=SI_mean_bio),method = "spearman",cor.coef.name = "rho", size = 8/.pt)+
  geom_smooth(aes(x=TE_mean_bio,y=SI_mean_bio),method = "lm",se = F, color = "#4daf4a")+
  labs(x="Translational error rate (e) (x10-4)",y="Lifespan")+
  scale_x_continuous(limit = c(6.6e-4, 12.6e-4), breaks = seq(7e-4, 11e-4, 2e-4), labels = seq(7, 11, 2))+
  scale_y_continuous(limits = c(1.9, 14), breaks = seq(3, 12, 3), labels = seq(3, 12, 3)) + 
  theme_classic()+theme(axis.text=element_text(size=9,face = "bold", color = "black"),
                        axis.title = element_text(size=10,face = "bold", color = "black"),
                        axis.title.x = element_blank(),axis.title.y = element_blank(),
                        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
                        legend.position = "none")

#figure S5F
fig.s5f.plot <- fig.s5f_data %>%
  ggplot()+
  geom_rect(data = chr_rect[seq(1,15,by=2),],mapping = aes(xmin = x1,xmax=x2, fill=c), ymin = -Inf,ymax = Inf,show.legend = F)+
  geom_point(aes(x=Marker_id,y=value,color=variable), size = 0.1)+
  labs(x="Genome position",y="LOD score")+
  scale_x_continuous(breaks = chr_marker_id_mean)+
  scale_y_continuous(limits = c(0,3.5))+
  scale_fill_manual(values=c("#BEBEBE"))+
  scale_color_manual(name="",values=c("#4daf4a","#2075bc"),label=c("Lifespan","Translational error rate (e)"))+
  guides(color=guide_legend(override.aes = list(size=3)))+
  theme_classic()+
  theme(axis.text=element_text(size=9,face = "bold",color = "black"),
        axis.title = element_text(size=10,face = "bold"),
        axis.text.x = element_text(angle = 20),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        legend.text = element_text(size=8,face = "bold"),legend.position = "none")



###save
ggsave(paste("./Fig.S5.A_", Sys.Date(), ".pdf", sep = ""), FIGS5.A, width = 10.8+0.75, height = 6, units = "cm")
ggsave(paste("./Fig.S5.B_", Sys.Date(), ".pdf", sep = ""), FigureS5.B, width = 3.6-0.5, height = 6-0.75, units = "cm")
ggsave(paste("./Fig.S5.C_", Sys.Date(), ".pdf", sep = ""), FIGS5.C, width = 3.6, height = 6-1.25, units = "cm")
ggsave(paste("./Fig.S5.D_", Sys.Date(), ".pdf", sep = ""), FIGS5.D, width = 3.6, height = 6-1.25, units = "cm")
ggsave(paste("./Fig.S5.E_", Sys.Date(), ".pdf", sep = ""), FigureS5.E, width = 3.6-0.5, height = 6-0.75, units = "cm")
ggsave(paste("./Fig.S5.F_", Sys.Date(), ".pdf", sep = ""), fig.s5f.plot, width = 10.8+0.75, height = 6, units = "cm")
