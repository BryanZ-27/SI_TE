###library
library(dplyr)
library(ggplot2)
library(cowplot)



###load data
load("./Data_for_FigS1.RData")



###main
##set colors
figs1_color <- c("#97cd96", "#70bb6f", "#4da74a", "#3b8239", "#285a2f")


##s1 alter D
#adjust y-axis
figs1_alterD_data <- figs1_alterD_data_final %>%
  filter(Ed_cutoff != 0.08) %>% 
  mutate(adjust_rho = -rho*30) %>% 
  group_by(Ed_cutoff, P) %>%
  summarise(mean_rho = mean(rho), sd_rho = sd(rho), 
            adjust_mean_rho = mean(adjust_rho)+2, adjust_sd_rho = sd(adjust_rho), 
            mean_p.value = mean(p.value), sd_p.value = sd(p.value)) %>%
  ungroup() %>% 
  mutate(adjust_mean_p.value = log10(mean_p.value)*3-2)

#draw figures with different D
fig.s1.D <- figs1_alterD_data %>% 
  ggplot() +
  geom_bar(aes(x = P, y = adjust_mean_rho, fill = as.factor(P)), 
           stat = "identity", position = position_dodge(width = 0.10), width = 0.10) +
  geom_bar(aes(x = P, y = adjust_mean_p.value, fill = as.factor(P)), 
           stat = "identity", position = position_dodge(width = 0.10), width = 0.01) +
  geom_errorbar(aes(x = P, ymin = adjust_mean_rho - adjust_sd_rho, 
                    ymax = adjust_mean_rho + adjust_sd_rho, color = as.factor(P)), 
                stat = "identity", position = position_dodge(width = 0.10), width = 0.05) + 
  geom_point(aes(x = P, y = adjust_mean_p.value, fill = as.factor(P)), 
             position = position_dodge(width = 0.10), size=2, color = rep(figs1_color, 4), shape = 21) +
  geom_abline(slope = 0, intercept = log10(0.05)*3-2, color = "red", linewidth = 0.8, linetype = 2) +
  geom_abline(slope = 0, intercept = 0, color = "white", linewidth = 5) +
  facet_grid(rows = vars(Ed_cutoff)) +
  geom_text(aes(label = paste("italic(D)==", Ed_cutoff, sep = "")), 
            x = -0.05, y = 20, parse = T, hjust = 0, size = 8/.pt) +
  scale_fill_manual(values = figs1_color) +
  scale_color_manual(values = figs1_color) +
  scale_x_continuous(limits = c(-0.1, 0.9), breaks = seq(0, 0.8, 0.2), labels = c("0%", "20%", "40%", "60%", "80%")) +
  scale_y_continuous(limits = c(-23, 24), 
                     breaks=c(-20, -14, -8, -2, 2, 8, 14, 20), 
                     labels = c(expression(10^-6), expression(10^-4), expression(10^-2), 0, 
                                0, -6/30, -12/30, -18/30)) +
  theme_classic() +
  theme(axis.text = element_text(size = 9, face = "bold", color = "black"),
        axis.title = element_blank(),
        axis.text.x = element_text(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        strip.text = element_blank(), 
        strip.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none")


##s1 alter alpha
#adjust y-axis
figs1_alterA_data <- figs1_alterA_data_final %>% 
  filter(Alpha != 1.8) %>% 
  mutate(adjust_rho = -rho*30) %>% 
  group_by(Alpha, P) %>%
  summarise(mean_rho = mean(rho), sd_rho = sd(rho), 
            adjust_mean_rho = mean(adjust_rho)+2, adjust_sd_rho = sd(adjust_rho), 
            mean_p.value = mean(p.value), sd_p.value = sd(p.value)) %>%
  ungroup() %>% 
  mutate(adjust_mean_p.value = log10(mean_p.value)*3-2)

#draw figures with different alpha
fig.s1.A <- figs1_alterA_data %>% 
  ggplot()+
  geom_bar(aes(x = P, y = adjust_mean_rho, fill = as.factor(P)), 
           stat = "identity", position = position_dodge(width = 0.10), width = 0.10) +
  geom_bar(aes(x = P, y = adjust_mean_p.value, fill = as.factor(P)), 
           stat = "identity", position = position_dodge(width = 0.10), width = 0.01) +
  geom_errorbar(aes(x = P, ymin = adjust_mean_rho - adjust_sd_rho, 
                    ymax = adjust_mean_rho + adjust_sd_rho, color = as.factor(P)), 
                stat = "identity", position = position_dodge(width = 0.10), width = 0.05) + 
  geom_point(aes(x = P, y = adjust_mean_p.value, fill = as.factor(P)), 
             position = position_dodge(width = 0.10), size = 2, color = rep(figs1_color, 4), shape=21) +
  geom_abline(slope = 0, intercept = log10(0.05)*3-2, color = "red", linewidth = 0.8, linetype = 2) +
  geom_abline(slope = 0, intercept = 0, color = "white", linewidth = 5) +
  facet_grid(rows = vars(Alpha)) +
  geom_text(aes(label = paste("alpha==", Alpha, sep = "")), x = -0.05, y = 20, parse = T, hjust = 0, size = 8/.pt) +
  scale_fill_manual(values = figs1_color) +
  scale_color_manual(values = figs1_color) +
  scale_x_continuous(limits = c(-0.1,0.9), breaks = seq(0, 0.8, 0.2), labels = c("0%", "20%", "40%", "60%", "80%")) +
  scale_y_continuous(limits = c(-23, 24), 
                     breaks = c(-20, -14, -8, -2, 2, 8, 14, 20), 
                     labels = c(expression(10^-6), expression(10^-4), expression(10^-2), 0, 
                                0, -6/30, -12/30, -18/30)) +
  theme_classic() +
  theme(axis.text=element_text(size=9,face = "bold", color = "black"),
        axis.text.x = element_text(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        strip.text = element_blank(), 
        strip.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none")


##s1 alter U
#adjust y-axis
figs1_alterU_data <- figs1_alterU_data_final %>%
  mutate(adjust_rho = -rho*30) %>% 
  group_by(Em, P) %>%
  summarise(mean_rho = mean(rho), sd_rho = sd(rho), 
            adjust_mean_rho = mean(adjust_rho)+2, adjust_sd_rho = sd(adjust_rho), 
            mean_p.value = mean(p.value), sd_p.value = sd(p.value)) %>%
  ungroup() %>% 
  mutate(adjust_mean_p.value = log10(mean_p.value)*3-2)

#draw figures with different U
fig.s1.U <- figs1_alterU_data %>% 
  ggplot() +
  geom_bar(aes(x = P, y = adjust_mean_rho, fill = as.factor(P)), 
           stat = "identity", position = position_dodge(width = 0.10), width = 0.10) +
  geom_bar(aes(x = P, y = adjust_mean_p.value, fill = as.factor(P)), 
           stat = "identity", position = position_dodge(width = 0.10), width = 0.01) +
  geom_errorbar(aes(x = P, ymin = adjust_mean_rho - adjust_sd_rho, 
                    ymax = adjust_mean_rho + adjust_sd_rho, color = as.factor(P)), 
                stat = "identity", position = position_dodge(width = 0.10), width = 0.05) + 
  geom_point(aes(x = P, y = adjust_mean_p.value, fill = as.factor(P)), 
             position = position_dodge(width = 0.10), size = 2, color = rep(figs1_color, 4), shape = 21) +
  geom_abline(slope = 0, intercept = log10(0.05)*3-2, color = "red", linewidth = 0.8, linetype = 2) +
  geom_abline(slope = 0, intercept = 0, color = "white", linewidth = 5) +
  facet_grid(rows = vars(Em)) +
  geom_text(aes(label = paste("italic(U)==", Em, sep = "")), x = -0.05, y = 20, parse = T, hjust = 0, size = 8/.pt) +
  scale_fill_manual(values = figs1_color) +
  scale_color_manual(values = figs1_color) +
  scale_x_continuous(limits = c(-0.1,0.9), breaks = seq(0, 0.8, 0.2), labels = c("0%", "20%", "40%", "60%", "80%")) +
  scale_y_continuous(limits = c(-23, 24), 
                     breaks = c(-20, -14, -8, -2, 2, 8, 14, 20), 
                     labels = c(expression(10^-6), expression(10^-4), expression(10^-2), 0, 
                                0, -6/30, -12/30, -18/30)) +
  theme_classic() +
  theme(axis.text = element_text(size = 9, face = "bold", color = "black"),
        axis.text.x = element_text(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        strip.text = element_blank(),
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")


##s1 alter L
#adjust y-axis
figs1_alterL_data <- figs1_alterL_data_final %>%
  filter(El != 7e-4) %>% 
  mutate(adjust_rho = -rho*30) %>% 
  group_by(El, P) %>%
  summarise(mean_rho = mean(rho), sd_rho = sd(rho), 
            adjust_mean_rho = mean(adjust_rho)+2, adjust_sd_rho = sd(adjust_rho), 
            mean_p.value = mean(p.value), sd_p.value = sd(p.value)) %>%
  ungroup() %>% 
  mutate(adjust_mean_p.value = log10(mean_p.value)*3-2)

#draw figures with different L
fig.s1.L <- figs1_alterL_data %>% 
  ggplot() +
  geom_bar(aes(x = P, y = adjust_mean_rho, fill = as.factor(P)), 
           stat = "identity", position = position_dodge(width = 0.10), width = 0.10) +
  geom_bar(aes(x = P, y = adjust_mean_p.value, fill = as.factor(P)), 
           stat = "identity", position = position_dodge(width = 0.10), width = 0.01) +
  geom_errorbar(aes(x = P, ymin = adjust_mean_rho - adjust_sd_rho, 
                    ymax = adjust_mean_rho + adjust_sd_rho, color = as.factor(P)), 
                stat = "identity", position = position_dodge(width = 0.10), width = 0.05) + 
  geom_point(aes(x = P, y = adjust_mean_p.value, fill = as.factor(P)), 
             position = position_dodge(width = 0.10), size = 2, color = rep(figs1_color, 4), shape = 21) +
  geom_abline(slope = 0, intercept = log10(0.05)*3-2, color = "red", linewidth = 0.8, linetype = 2) +
  geom_abline(slope = 0, intercept = 0, color = "white", linewidth = 5) +
  facet_grid(rows = vars(El)) +
  geom_text(aes(label = paste("italic(L)==", El, sep = "")), x = -0.05, y = 20, parse = T, hjust = 0, size = 8/.pt) +
  scale_fill_manual(values = figs1_color) +
  scale_color_manual(values = figs1_color) +
  scale_x_continuous(limits = c(-0.1,0.9), breaks = seq(0, 0.8, 0.2), labels = c("0%", "20%", "40%", "60%", "80%")) +
  scale_y_continuous(limits = c(-23, 24), 
                     breaks=c(-20, -14, -8, -2, 2, 8, 14, 20), 
                     labels = c(expression(10^-6), expression(10^-4), expression(10^-2), 0, 
                                0, -6/30, -12/30, -18/30))+
  theme_classic() +
  theme(axis.text = element_text(size = 9, face = "bold", color = "black"),
        axis.text.x = element_text(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        strip.text = element_blank(),
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")


##combine figures
fig.s1 <- ggdraw()+
  draw_plot(fig.s1.D, x = 0, y = 0, width = 0.25, height = 1) +
  draw_plot(fig.s1.A, x = 0.25, y = 0, width = 0.25, height = 1) +
  draw_plot(fig.s1.U, x = 0.50, y = 0, width = 0.25, height = 1) +
  draw_plot(fig.s1.L, x = 0.75, y = 0, width = 0.25, height = 1)



###save
ggsave(paste("./FigS1_", Sys.Date(), ".pdf", sep = ""), fig.s1, width = 17, height = 22, units = "cm")
