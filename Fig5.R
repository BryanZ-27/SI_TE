###library
library(dplyr)
library(ggplot2)



###load data
load("./Data_for_Fig5.RData")



###main
#figure A and B
FIG5.A <- fig.5ab_data %>%
  ggplot()+
  geom_rect(data = chr_rect[seq(1,15,by=2),],mapping = aes(xmin = x1,xmax=x2, fill=c), ymin = -Inf,ymax = Inf,show.legend = F)+
  geom_point(aes(x=Marker_id,y=value,color=variable), size = 0.1)+
  geom_hline(yintercept = c(3.59,3.49),color = c("#4daf4a","#2075bc"),linetype = 2)+
  labs(x="Genome position",y="LOD score")+
  scale_x_continuous(breaks = chr_marker_id_mean)+
  scale_y_continuous(limits = c(0,4.5))+
  scale_fill_manual(values=c("#BEBEBE"))+
  scale_color_manual(name="",values=c("#4daf4a","#2075bc"),label=c("Lifespan","Translational error rate (e)"))+
  guides(color=guide_legend(override.aes = list(size=3)))+
  theme_classic()+
  theme(axis.text=element_text(size=9,face = "bold",color = "black"),
        axis.title = element_text(size=10,face = "bold"),
        axis.text.x = element_text(angle = 20),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        legend.text = element_text(size=8,face = "bold"),legend.position = "none")

FIG5.B <- fig.5ab_data %>%
  filter(chr == 10) %>%
  mutate(Marker_pos = as.numeric(gsub("\\d+_chr[IXV+]_(\\d+)_.*","\\1",Marker))) %>%
  filter(Marker_pos <= 700e3, Marker_pos >= 600e3) %>% 
  ggplot()+
  geom_rect(data = Interval_rect,mapping = aes(xmin = x1,xmax=x2, fill=c), ymin = -Inf,ymax = Inf,show.legend = F)+
  geom_point(aes(x=Marker_pos,y=value,color=variable), size = 0.4)+
  geom_hline(yintercept = c(3.59,3.49),color = c("#4daf4a","#2075bc"),linetype = 2)+
  labs(x="Genome position in chrX",y="LOD score")+
  scale_x_continuous(breaks = c(seq(0,600000,by=100000),641753,669427,700000),labels = c("0","100k","200k","300k","400k","500k","600k","641k","669k","700k"))+
  scale_y_continuous(limits = c(0,4.5))+
  scale_fill_manual(values=c("#BEBEBE"))+
  scale_color_manual(name="",values=c("#4daf4a","#2075bc"),label=c("Lifespan","Translational error rate (e)"))+
  guides(color=guide_legend(override.aes = list(size=3)))+
  theme_classic()+
  theme(axis.text=element_text(size=9,face = "bold", color = "black"),
        axis.title = element_text(size=10,face = "bold"),
        axis.text.x = element_text(),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        legend.text = element_text(size=8,face = "bold"),legend.position = "none")

#figure C
FIG5.C <- fig.5c_data %>%
  na.omit() %>%
  ggplot(group=!!sym("6487936_chrX_655475_C_T"))+
  geom_point(aes(x=!!sym("6487936_chrX_655475_C_T"),y=pheno,color = !!sym("6487936_chrX_655475_C_T")),position = position_jitter(w = 0.2, h = 0), 
             size = 0.3)+
  geom_errorbar(aes(ymin=fig.5c_summ$m[1]-fig.5c_summ$se[1],ymax=fig.5c_summ$m[1]+fig.5c_summ$se[1],x=1),
                size=0.3,width=0.2)+
  geom_errorbar(aes(ymin=fig.5c_summ$m[2]-fig.5c_summ$se[2],ymax=fig.5c_summ$m[2]+fig.5c_summ$se[2],x=2),
                size=0.3,width=0.2)+
  geom_segment(aes(y = fig.5c_summ$m[1], yend = fig.5c_summ$m[1], x = 0.75, xend = 1.25), size = 0.3)+
  geom_segment(aes(y = fig.5c_summ$m[2], yend = fig.5c_summ$m[2], x = 1.75, xend = 2.25), size = 0.3)+
  scale_x_discrete(labels= c("BY allele","RM allele"))+
  scale_color_manual(values = c("#e5d851", "#8c07bb"))+
  scale_y_continuous(limits=c(6.5,13), breaks = seq(7, 13, 2), labels = seq(7, 13, 2))+
  labs(title = "6487936_chrX_655475_C_T",x="Genotype",y="Lifespan")+
  theme_classic()+theme(axis.text=element_text(size=9,face = "bold",color = "black"),
                        axis.title = element_text(size=10,face = "bold",color = "black"),
                        axis.text.x = element_text(vjust = 1),
                        axis.title.x = element_blank(),plot.title = element_blank(),
                        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
                        legend.position = "none")

#figure D
FIG5.D <- fig.5d_data %>%
  na.omit() %>%
  ggplot(group=!!sym("6487936_chrX_655475_C_T"))+
  geom_point(aes(x=!!sym("6487936_chrX_655475_C_T"),y=pheno,color = !!sym("6487936_chrX_655475_C_T")),position = position_jitter(w = 0.2, h = 0), 
             size = 0.3)+
  geom_errorbar(aes(ymin=fig.5d_summ$m[1]-fig.5d_summ$se[1],ymax=fig.5d_summ$m[1]+fig.5d_summ$se[1],x=1),
                size=0.3,width=0.2)+
  geom_errorbar(aes(ymin=fig.5d_summ$m[2]-fig.5d_summ$se[2],ymax=fig.5d_summ$m[2]+fig.5d_summ$se[2],x=2),
                size=0.3,width=0.2)+
  geom_segment(aes(y = fig.5d_summ$m[1], yend = fig.5d_summ$m[1], x = 0.75, xend = 1.25), size = 0.3)+
  geom_segment(aes(y = fig.5d_summ$m[2], yend = fig.5d_summ$m[2], x = 1.75, xend = 2.25), size = 0.3)+
  scale_x_discrete(labels= c("BY allele","RM allele"))+
  scale_color_manual(values = c("#e5d851", "#8c07bb"))+
  scale_y_continuous(limits=c(0.00065,0.0013), breaks = seq(7e-4, 13e-4, 2e-4), labels = seq(7e-4, 13e-4, 2e-4))+
  labs(title = "6487936_chrX_655475_C_T",x="Genotype",y="Translational error rate")+
  theme_classic()+theme(axis.text=element_text(size=9,face = "bold",color = "black"),
                        axis.title = element_text(size=10,face = "bold",color = "black"),
                        axis.text.x = element_text(vjust = 1),
                        axis.title.x = element_blank(),plot.title = element_blank(),
                        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
                        legend.position = "none")

#figure E and F
FIG5.E <- fig.5e_data %>% 
  mutate(strains = ifelse(strains %in% c("BY4716"),"BY",ifelse(strains %in% c("RM11"),"RM",strains))) %>%
  dplyr::mutate(strains_sort = SI_Geneverify_gene_name_sort[[strains]]) %>%
  mutate(strains_col = case_when(strains %in% c("BY") ~ 1,
                                 strains %in% c("RM") ~ 99,
                                 strains %in% c("VPS70") ~ 3,
                                 TRUE ~ 2)) %>%
  ggplot(aes(x=reorder(strains,strains_sort),y=SI_mean_tube))+
  geom_boxplot(aes(fill = as.factor(strains_col)),show.legend = F, size = 0.4)+
  stat_compare_means(aes(label = ..p.format..),comparisons = SI_Geneverify_gene_comparisons_list,
                     method = "wilcox.test",method.args = list(alternative = "less"), size = 8/.pt)+
  scale_fill_manual(name="Strain",values=c("#e5d851","white","#ffb6c1","#8c07bb"))+
  scale_x_discrete(labels=c('BY', 'BY::ILM1-RM','BY::JHD2-RM','BY::STE24-RM','BY::VPS70-RM',
                            'BY::RSF2-YJR128W-RM','RM'))+
  labs(x="",y="Lifespan")+
  theme_classic()+theme(axis.text=element_text(size=9,face = "bold", color = "black"),
                        axis.title = element_text(size=10,face = "bold", color = "black"),
                        axis.text.x = element_blank(),
                        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
                        legend.text = element_text(size=15,face = "bold"),legend.title = element_text(size=15,face = "bold"))

FIG5.F <- fig.5f_data %>%
  ggplot(aes(x=reorder(strain,strain_sort),y=TE))+
  geom_boxplot(aes(fill = as.factor(strain_col)),show.legend = F, size = 0.4)+
  stat_compare_means(aes(label = ..p.format..),comparisons = TE_Geneverify_gene_comparisons_list,
                     method = "wilcox.test",method.args = list(alternative = "greater"), size = 8/.pt)+
  scale_fill_manual(name="Strain",values=c("#e5d851","#ffb6c1","#8c07bb"))+
  scale_y_continuous(breaks = c(seq(0.0006, 0.0015, 0.0003)), labels = c(seq(6, 15, 3)))+
  scale_x_discrete(labels=c('BY', 'BY::VPS70-RM','RM'))+
  labs(x="",y="Translational error rate")+
  theme_classic()+theme(axis.text=element_text(size=9,face = "bold", color = "black"),
                        axis.title = element_text(size=10,face = "bold", color = "black"),
                        axis.text.x = element_blank(),
                        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
                        legend.text = element_text(size=15,face = "bold"),legend.title = element_text(size=15,face = "bold"))

#figure G
#translation error rate
fig.5g1 <- ggplot(fig.5g1_data, aes(Strain, normal_TE, fill = Strain)) + 
  geom_boxplot(size = 0.4) + 
  geom_jitter(width = 0.1, size = 0.8) + 
  facet_wrap(~treatment, labeller = as_labeller(c(`1` = "DMSO", ConA = "ConA"))) + 
  scale_fill_manual(values = c("#dcd154", "#ffb6c1")) + 
  scale_y_continuous(limits = c(0.66, 1.3), breaks = c(0.7, 0.9, 1.1, 1.3))+
  theme_classic()+
  theme(axis.text= element_text(size=9,face = "bold", color = "black"),
        axis.title = element_text(size=10,face = "bold", color = "black"),
        axis.text.x = element_blank(), axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        legend.background = element_blank(),
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_blank(),
        strip.text = element_text(size=8,color = "black",face = "bold"))

#chronological lifespan
fig.5g2 <- ggplot(fig.5g2_data, aes(Strain, normal_SI, fill = Strain)) + 
  geom_boxplot(size = 0.4) + 
  geom_jitter(width = 0.1, size = 0.8) + 
  facet_wrap(~treatment, labeller = as_labeller(c(`1` = "DMSO", ConA = "ConA"))) + 
  scale_fill_manual(values = c("#dcd154", "#ffb6c1")) + 
  scale_y_continuous(limits = c(0.83, 1.7), breaks = seq(0.9, 1.7, 0.2))+
  theme_classic()+
  theme(axis.text= element_text(size=9,face = "bold", color = "black"),
        axis.title = element_text(size=10,face = "bold", color = "black"),
        axis.text.x = element_blank(), axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        legend.background = element_blank(),
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_blank(),
        strip.text = element_text(size=8,color = "black",face = "bold"))



###save
ggsave(paste("./Fig.5.A_", Sys.Date(), ".pdf", sep = ""), FIG5.A, width = 13+1, height = 6, units = "cm")
ggsave(paste("./Fig.5.B_", Sys.Date(), ".pdf", sep = ""), FIG5.B, width = 9, height = 6-0.5, units = "cm")
ggsave(paste("./Fig.5.C_", Sys.Date(), ".pdf", sep = ""), FIG5.C, width = 4.5, height = 6-1.5, units = "cm")
ggsave(paste("./Fig.5.D_", Sys.Date(), ".pdf", sep = ""), FIG5.D, width = 4.5, height = 6-1.5, units = "cm")
ggsave(paste("./Fig.5.E_", Sys.Date(), ".pdf", sep = ""), FIG5.E, width = 8, height = 6-1, units = "cm")
ggsave(paste("./Fig.5.F_", Sys.Date(), ".pdf", sep = ""), FIG5.F, width = 4, height = 6-0.5, units = "cm")
ggsave(paste("./Fig.5g1_", Sys.Date(), ".pdf", sep = ""), fig.5g1, width = 6-0.5, height = 6-0.5, units = "cm")
ggsave(paste("./Fig.5g2_", Sys.Date(), ".pdf", sep = ""), fig.5g2, width = 6-0.5, height = 6-0.5, units = "cm")
