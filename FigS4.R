###library
library(dplyr)
library(ggplot2)



###load data
load("./Data_for_FigS4.RData")



###main
Fig.S4 <- fig.s4_data %>%
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
        strip.text = element_text(size=8,color = "black",face = "bold"))



###save
ggsave(paste("./Fig.S4_", Sys.Date(), ".pdf", sep = ""), Fig.S4, width = 12, height = 8-0.5, units = "cm")
