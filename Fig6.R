###library
library(ggplot2)
library(ggbreak)



###load data
load("./Data_for_Fig6.RData")



###main
fig.6a <- ggplot(fig.6a_data, aes(BP, BONF)) + 
  annotate("rect", xmin = 4916e04, xmax = 8993e04, ymin = 1, ymax = 1.05, fill = "#666666") + 
  annotate("rect", xmin = 4916e04, xmax = 8993e04, ymin = 0.5, ymax = 1, fill = "#858585") + 
  annotate("rect", xmin = 4916e04, xmax = 8993e04, ymin = 0.1, ymax = 0.5, fill = "#a3a3a3") + 
  annotate("rect", xmin = 4916e04, xmax = 8993e04, ymin = 0.05, ymax = 0.1, fill = "#c2c2c2") + 
  annotate("rect", xmin = 4916e04, xmax = 8993e04, ymin = 0, ymax = 0.05, fill = "#e0e0e0") + 
  geom_point(aes(), size = 0.5) + 
  geom_hline(yintercept = -0.05, color = "#31c831") + 
  scale_x_continuous(limits = c(4916e04, 8993e04)) + 
  scale_x_break(c(4924e04, 8986e04)) + 
  scale_y_reverse(breaks = c(0.05, 0.1, 0.5, 1)) + 
  labs(x="BP",y="Bonferroni adjusted P")+
  theme_classic()+
  theme(axis.text = element_blank(), 
        axis.title = element_blank(), 
        legend.position = "none")
for (i in 1:19){
  fig.6a <- fig.6a+
    annotate("rect", xmin = exons_NAALAD2[i*2-1], xmax = exons_NAALAD2[i*2], ymin = -0.01, ymax = -0.09, fill = "#134e13")+
    annotate("rect", xmin = exons_FOLH1[i*2-1], xmax = exons_FOLH1[i*2], ymin = -0.01, ymax = -0.09, fill = "#134e13")
}

fig.6b <- ggplot(fig.6b_data, aes(BP, BONF)) + 
  annotate("rect", xmin = 4916e04, xmax = 8993e04, ymin = 1, ymax = 1.05, fill = "#666666") + 
  annotate("rect", xmin = 4916e04, xmax = 8993e04, ymin = 0.5, ymax = 1, fill = "#858585") + 
  annotate("rect", xmin = 4916e04, xmax = 8993e04, ymin = 0.1, ymax = 0.5, fill = "#a3a3a3") + 
  annotate("rect", xmin = 4916e04, xmax = 8993e04, ymin = 0.05, ymax = 0.1, fill = "#c2c2c2") + 
  annotate("rect", xmin = 4916e04, xmax = 8993e04, ymin = 0, ymax = 0.05, fill = "#e0e0e0") + 
  geom_point(aes(), size = 0.5) + 
  geom_hline(yintercept = -0.05, color = "#31c831") + 
  scale_x_continuous(limits = c(4916e04, 8993e04)) + 
  scale_x_break(c(4924e04, 8986e04)) + 
  scale_y_reverse(breaks = c(0.05, 0.1, 0.5, 1)) + 
  labs(x="BP",y="Bonferroni adjusted P")+
  theme_classic()+
  theme(axis.text = element_blank(), 
        axis.title = element_blank(), 
        legend.position = "none")
for (i in 1:19){
  fig.6b <- fig.6b+
    annotate("rect", xmin = exons_NAALAD2[i*2-1], xmax = exons_NAALAD2[i*2], ymin = -0.01, ymax = -0.09, fill = "#134e13")+
    annotate("rect", xmin = exons_FOLH1[i*2-1], xmax = exons_FOLH1[i*2], ymin = -0.01, ymax = -0.09, fill = "#134e13")
}

fig.6c <- ggplot(fig.6c_data, aes(Allele, lifespan)) + 
  geom_boxplot(aes(fill = Allele), size = 0.5, color = "black")+
  geom_jitter(width = 0.2, size = 0.8)+
  scale_fill_manual(values = c("#8ecf8c", "#357833"))+
  labs(x="Allele",y="Lifespan") +
  theme_classic()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
                        axis.text = element_blank(), axis.title = element_blank(), 
                        legend.position = "none")



###save
ggsave(paste("./Fig.6.A_", Sys.Date(), ".pdf", sep = ""), fig.6a, width = 9+1, height = 5+1, units = "cm")
ggsave(paste("./Fig.6.B_", Sys.Date(), ".pdf", sep = ""), fig.6b, width = 9+1, height = 5+1, units = "cm")
ggsave(paste("./Fig.6.C_", Sys.Date(), ".pdf", sep = ""), fig.6c, width = 5.9+0.1, height = 10.8+0.2, units = "cm")
