###library
library(ggplot2)



###load data
load("./Data_for_FigS2.RData")



###main
FigureS2  <- ggplot()+
  geom_point(aes(mean_score, no, color = mean_score), fig.s2a_summ, size = 2)+
  geom_errorbar(aes(xmin = mean_score - sd_score, 
                    xmax = ifelse(mean_score + sd_score > 1, 1, mean_score + sd_score), 
                    y = no, color = mean_score), fig.s2a_summ, width = 10, linewidth = 0.5)+
  geom_freqpoly(aes(value), fig.s2b_data, color = "#aadaa8", binwidth = 0.02) + 
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



###save
ggsave(paste("./Fig.S2_", Sys.Date(), ".pdf", sep = ""), FigureS2, width = 8, height = 9, units = "cm")
