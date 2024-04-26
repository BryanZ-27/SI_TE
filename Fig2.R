###library
library(dplyr)
library(ggplot2)


###load data
load("./Data_for_Fig2.RData")


###main
#adjust y-axis
data_for_fig2_adjusted <- data_for_fig2_final %>%
  mutate(adjust_rho = -rho*20) %>% 
  group_by(P) %>%
  summarise(mean_rho = mean(rho), sd_rho = sd(rho), 
            adjust_mean_rho = mean(adjust_rho) + 0.5, adjust_sd_rho = sd(adjust_rho), 
            mean_p.value = mean(p.value), sd_p.value = sd(p.value)) %>%
  ungroup() %>% 
  mutate(adjust_mean_p.value = log10(mean_p.value)*3-0.5)

#set colors
color_for_fig2 <- c("#90cb9c", "#70bd79", "#57aa56", "#4d8544", "#375b33")

#draw main figure
fig2.e <- data_for_fig2_adjusted %>% 
  ggplot() +
  geom_bar(aes(x = P, y = adjust_mean_rho, fill = as.factor(P)), 
           stat = "identity", position = position_dodge(width = 0.15), width = 0.10) +
  geom_bar(aes(x = P,y = adjust_mean_p.value, fill = as.factor(P)), 
           stat = "identity", position = position_dodge(width = 0.15), width = 0.01) +
  geom_errorbar(aes(x = P, ymin = adjust_mean_rho, ymax = adjust_mean_rho + adjust_sd_rho, color = as.factor(P)), 
                stat = "identity", position = position_dodge(width = 0.15), width = 0.10) + 
  geom_point(aes(x = P, y = adjust_mean_p.value, fill = as.factor(P)), 
             position = position_dodge(width = 0.15), size = 30, color = color_for_fig2, shape = 21) +
  geom_abline(slope = 0, intercept = log10(0.05)*3-0.5, color = "red", linewidth = 1.2, linetype = 2) +
  annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymax = 0.5, ymin = -0.5, fill = "white") +
  annotate(geom = "text", y = 0, x = data_for_fig2_final$P %>% unique(), label = data_for_fig2_final$P %>% unique(), size = 4.5) +
  annotate(geom = "text", label = "alpha==1.5", x = 0.1, y = 12, size = 4.5, parse = TRUE) +
  annotate(geom = "text", label = "~italic(D)==0.06", x = 0.1, y = 11, size = 4.5, parse = TRUE) +
  annotate(geom = "text", label = "~~~~~italic(U)==2%*%10^-3", x = 0.1, y = 10, size = 4.5, parse = TRUE) +
  annotate(geom = "text", label = "~~~~~italic(L)==5%*%10^-4", x = 0.1, y = 9, size = 4.5, parse = TRUE) +
  ylab("Spearman’s correlation between translational\nerror rate and lifespan\nP value(log10 scale)\trho") +
  scale_fill_manual(values = color_for_fig2) +
  scale_color_manual(values = color_for_fig2) +
  scale_y_continuous(limits = c(-14, 14), 
                     breaks = c(-12.5, -9.5, -6.5, -3.5, -0.5, 0.5, 3.5, 6.5, 9.5, 12.5), 
                     labels = c(expression(10^-4), expression(10^-3), expression(10^-2), expression(10^-1), 0, 
                                0, "-0.15", "-0.30", "-0.45", "-0.60")) +
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size=16, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank())

#draw zoom1 figure
fig2.e_zoom1 <- data_for_fig2_zoom1 %>% 
  ggplot() +
  geom_point(aes(x = Fake_E, y = Fake_T), size = 1) +
  geom_smooth(aes(x = Fake_E, y = Fake_T), method = "lm", se = F, color = color_for_fig2[1]) +
  scale_x_continuous(limits = c(0.0005, 0.0020), breaks = c(seq(0.0005, 0.0020, 0.0005)), labels = c(seq(5, 20, 5)))+
  scale_y_continuous(breaks = c(seq(0, 10, 5)), labels = c(seq(0, 10, 5)), limits = c(0, 10))+
  labs(x = "", y = "Lifespan") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5,size = 18, face = "bold"), 
        axis.text = element_text(size = 14, face = "bold"), 
        axis.title = element_text(size = 16, face = "bold"), 
        axis.text.x = element_text(vjust = 1, hjust = 1), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

#draw zoom2 figure
fig2.e_zoom2 <- data_for_fig2_zoom2 %>% 
  ggplot() +
  geom_point(aes(x = Fake_E, y = Fake_T), size = 1) +
  geom_smooth(aes(x = Fake_E, y = Fake_T), method = "lm", se = F, color = color_for_fig2[3]) +
  scale_x_continuous(limits = c(0.0005, 0.0020), breaks = c(seq(0.0005, 0.0020, 0.0005)), labels = c(seq(5, 20, 5)))+
  scale_y_continuous(breaks = c(seq(0, 10, 5)), labels = c(seq(0, 10, 5)), limits = c(0, 10))+
  labs(x = "", y = "Lifespan") +
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(vjust = 1, hjust = 1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

#draw zoom3 figure
fig2.e_zoom3 <- data_for_fig2_zoom3 %>% 
  ggplot()+
  geom_point(aes(x = Fake_E, y = Fake_T), size = 1) +
  geom_smooth(aes(x = Fake_E, y = Fake_T), method = "lm", se = F, color = color_for_fig2[5]) +
  scale_x_continuous(limits = c(0.0005, 0.0020), breaks = c(seq(0.0005, 0.0020, 0.0005)), labels = c(seq(5, 20, 5))) +
  scale_y_continuous(breaks = c(seq(0, 10, 5)), labels = c(seq(0, 10, 5)), limits = c(0, 10)) +
  labs(x="", y="Lifespan") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5,size = 18, face = "bold"),
        axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(vjust = 1, hjust = 1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


###save
ggsave(paste("./Fig2.e_", Sys.Date(), ".pdf", sep = ""), fig2.e, width = 32, height = 24, units = "cm")
ggsave(paste("./Fig2.e_zoom1_", Sys.Date(), ".pdf", sep = ""), fig2.e_zoom1, width = 8, height = 6, units = "cm")
ggsave(paste("./Fig2.e_zoom2_", Sys.Date(), ".pdf", sep = ""), fig2.e_zoom2, width = 8, height = 6, units = "cm")
ggsave(paste("./Fig2.e_zoom3_", Sys.Date(), ".pdf", sep = ""), fig2.e_zoom3, width = 8, height = 6, units = "cm")