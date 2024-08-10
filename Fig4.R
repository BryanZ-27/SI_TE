###library
library(dplyr)
library(ggplot2)



###load data
load("./Data_for_Fig4.RData")


###main
fig4_color <- c("#9cd49a", "#73c371", "#4daf4a", "#3b8639", "#295c27", "#163316")

FIG4.A <- fig.4a_data %>%
  ggplot()+
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

shift <- 0.5

FIG4.B <- fig.4b_data %>% 
  ggplot()+
  geom_bar(aes(x=order,y=-estimate*20+shift, fill = as.factor(cut_off)), 
           stat = "identity",position =position_dodge(width = 0.15),width = 0.12)+
  geom_bar(aes(x=order,y=log10(p.value)*4-shift, fill = as.factor(cut_off)), 
           stat = "identity",position =position_dodge(width = 0.15),width = 0.02)+
  geom_point(aes(x=order,y=log10(p.value)*4-shift, fill = as.factor(cut_off)), 
             position = position_dodge(width = 0.15), size=18, color = fig4_color, shape=21)+
  geom_abline(slope = 0,intercept=log10(0.05)*4-shift,color="red",linewidth=1.2,linetype=2)+
  annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymax = shift, ymin = -shift,fill = "white") +
  annotate(geom = "text", y = 0, x = fig.4b_data$order %>% unique(), 
           label = fig.4b_data$cut_off %>% unique(), size = 9/.pt, color = "black", fontface = "bold") +
  scale_y_continuous(limits = c(-14, 14), 
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
        axis.line.x = element_blank())



###save
ggsave(paste("./Fig.4.A_", Sys.Date(), ".pdf", sep = ""), FIG4.A, width = 9, height = 6, units = "cm")
ggsave(paste("./Fig.4.B_", Sys.Date(), ".pdf", sep = ""), FIG4.B, width = 18, height = 15, units = "cm")
