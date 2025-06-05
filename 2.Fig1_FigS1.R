###library
library(ggplot2)
library(ggpubr)
library(cowplot)


###load
load("./2.cor_filter_P_data_200.Rdata")


###main
##fig1e
fig1e_data_original <- cor_filter_P_data_200 %>% 
  ungroup() %>% 
  filter(Alpha==1.5,Ed_cutoff == 0.06,El==5*10^-4, Em == 2e-3, P == 0) %>% 
  select(1:8)

fig1e_data_final <- tibble(Ed_cutoff = 0.06,Alpha=1.5,e0 = 0,El=5*10^-4, Em = 2e-3, P = 0, 
                           seed_number = NA, rho = NA, p.value = NA)
for(i in 1 : nrow(fig1e_data_original)){
  df <- fig1e_data_original$data[i]
  rho <- cor.test(df[[1]]$Fake_T, df[[1]]$Fake_E, method = "spearman")$estimate
  p.value <- cor.test(df[[1]]$Fake_T, df[[1]]$Fake_E, method = "spearman")$p.value
  fig1e_data_final <- rbind(fig1e_data_final, tibble(Ed_cutoff = 0.06,Alpha=1.5,e0 = 0,El=5*10^-4, Em = 2e-3, P = 0, 
                                                     seed_number = i, rho = rho, p.value = p.value)) %>% 
    filter(!is.na(rho))
  for (ratio in seq(0.2, 0.8, 0.2)) {
    new_data <- df[[1]] %>% filter(Fake_T >= quantile(Fake_T, ratio))
    new_ratio <- tibble(Ed_cutoff = 0.06,Alpha=1.5,e0 = 0,El=5*10^-4, Em = 2e-3, P = ratio, seed_number = i, 
                        rho = cor.test(new_data$Fake_T, new_data$Fake_E, method = "spearman")$estimate, 
                        p.value = cor.test(new_data$Fake_T, new_data$Fake_E, method = "spearman")$p.value)
    fig1e_data_final <- rbind(fig1e_data_final, new_ratio)
  }
}

fig1e_data <- fig1e_data_final %>%
  mutate(adjust_rho = -rho*20) %>% 
  group_by(P) %>%
  dplyr::summarise(mean_rho = mean(rho), sd_rho = sd(rho), 
                   adjust_mean_rho = mean(adjust_rho) + 0.5, adjust_sd_rho = sd(adjust_rho), 
                   mean_p.value = mean(p.value)*5, sd_p.value = sd(p.value)) %>%
  ungroup() %>% 
  mutate(adjust_mean_p.value = log10(mean_p.value)*3-0.5)

fig1e_color <- c("#99cc99", "#66cc66", "#669933", "#339933", "#336633")

FIG.1.E <- fig1e_data %>% ggplot()+
  geom_bar(aes(x=P,y=adjust_mean_rho, fill = as.factor(P)), 
           stat = "identity",position =position_dodge(width = 0.15),width = 0.10)+
  geom_bar(aes(x=P,y=adjust_mean_p.value, fill = as.factor(P)), 
           stat = "identity",position =position_dodge(width = 0.15),width = 0.01)+
  geom_errorbar(aes(x = P, ymin = adjust_mean_rho - adjust_sd_rho, ymax = adjust_mean_rho + adjust_sd_rho), 
                stat = "identity",position = position_dodge(width = 0.15),width = 0.10, color = "black") + 
  geom_point(aes(x=P,y=adjust_mean_p.value, fill = as.factor(P)), 
             position = position_dodge(width = 0.15), size=30, color = fig1e_color, shape=21)+
  geom_abline(slope = 0,intercept=log10(0.05)*3-0.5,color="red",linewidth=1.2,linetype=2)+
  annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymax = 0.5, ymin = -0.5,fill = "white") +
  annotate(geom = "text", y = 0, x = fig1e_data$P %>% unique(), label = fig1e_data$P %>% unique(), size = 4.5) +
  annotate(geom = "text",label="alpha==1.5",x = 0.1,y = 12,size=4.5,parse = TRUE)+
  annotate(geom = "text",label="~italic(D)==0.06",x = 0.1,y = 11,size=4.5,parse = TRUE)+
  annotate(geom = "text",label="~~~~~italic(U)==2%*%10^-3",x = 0.1,y = 10,size=4.5,parse = TRUE)+
  annotate(geom = "text",label="~~~~~italic(L)==5%*%10^-4",x = 0.1,y = 9,size=4.5,parse = TRUE)+
  ylab("Spearman’s correlation between translational\nerror rate and lifespan\nP value(log10 scale)\trho")+
  ###plot modify###
  scale_fill_manual(values = fig1e_color)+
  scale_color_manual(values = fig1e_color)+
  scale_y_continuous(limits = c(-13, 13), 
                     breaks=c(-12.5, -9.5, -6.5, -3.5, -0.5, 0.5, 3.5, 6.5, 9.5 ,12.5), 
                     labels = c(expression(10^-4), expression(10^-3), 
                                expression(10^-2), expression(10^-1), 0, 
                                0, "-0.15", "-0.30", "-0.45", "-0.60"))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5,size = 18, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title = element_text(size=16,face = "bold"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank()
  )

ggsave(paste("./Fig.1.E_", Sys.Date(), ".pdf", sep = ""), FIG.1.E, width = 32, height = 24, units = "cm")

#zoom
zoom1_data <- fig1e_data_original %>% 
  filter(seed_number == 618) %>% 
  pull(data) %>% as.data.frame()

zoom1_plot <- zoom1_data %>% ggplot()+
  geom_point(aes(x=Fake_E,y=Fake_T), size = 1)+
  geom_smooth(aes(x=Fake_E,y=Fake_T),method = "lm",se = F, color = fig1e_color[1])+
  scale_x_continuous(limits = c(0.0005, 0.0020), breaks = c(seq(0.0005, 0.0020, 0.0005)), labels = c(seq(5, 20, 5)))+
  scale_y_continuous(breaks = c(seq(0, 10, 5)), labels = c(seq(0, 10, 5)), limits = c(3, 13))+
  labs(x="",
       y="Lifespan")+
  theme_classic()+theme(plot.title = element_text(hjust = 0.5,size = 18, face = "bold"),
                        axis.text=element_text(size=14,face = "bold"),
                        axis.title = element_text(size=16,face = "bold"),
                        axis.text.x = element_text(vjust = 1,hjust = 1),
                        panel.grid.major=element_blank(),panel.grid.minor=element_blank())

zoom2_data <- fig1e_data_original %>% 
  filter(seed_number == 618) %>% 
  pull(data) %>% as.data.frame() %>% filter(Fake_T >= quantile(Fake_T, 0.4))

zoom2_plot <- zoom2_data %>% ggplot()+
  geom_point(aes(x=Fake_E,y=Fake_T), size = 1)+
  geom_smooth(aes(x=Fake_E,y=Fake_T),method = "lm",se = F, color = fig1e_color[3])+
  scale_x_continuous(limits = c(0.0005, 0.0020), breaks = c(seq(0.0005, 0.0020, 0.0005)), labels = c(seq(5, 20, 5)))+
  scale_y_continuous(breaks = c(seq(0, 10, 5)), labels = c(seq(0, 10, 5)), limits = c(3, 13))+
  labs(x="",y="Lifespan")+
  theme_classic()+theme(plot.title = element_text(hjust = 0.5,size = 18, face = "bold"),
                        axis.text=element_text(size=14,face = "bold"),
                        axis.title = element_text(size=16,face = "bold"),
                        axis.text.x = element_text(vjust = 1,hjust = 1),
                        panel.grid.major=element_blank(),panel.grid.minor=element_blank())

zoom3_data <- fig1e_data_original %>% 
  filter(seed_number == 618) %>% 
  pull(data) %>% as.data.frame() %>% filter(Fake_T >= quantile(Fake_T, 0.8))

zoom3_plot <- zoom3_data %>% ggplot()+
  geom_point(aes(x=Fake_E,y=Fake_T), size = 1)+
  geom_smooth(aes(x=Fake_E,y=Fake_T),method = "lm",se = F, color = fig1e_color[5])+
  scale_x_continuous(limits = c(0.0005, 0.0020), breaks = c(seq(0.0005, 0.0020, 0.0005)), labels = c(seq(5, 20, 5)))+
  scale_y_continuous(breaks = c(seq(0, 10, 5)), labels = c(seq(0, 10, 5)), limits = c(3, 13))+
  labs(x="",
       y="Lifespan")+
  theme_classic()+theme(plot.title = element_text(hjust = 0.5,size = 18, face = "bold"),
                        axis.text=element_text(size=14,face = "bold"),
                        axis.title = element_text(size=16,face = "bold"),
                        axis.text.x = element_text(vjust = 1,hjust = 1),
                        panel.grid.major=element_blank(),panel.grid.minor=element_blank())

ggsave(paste("./Fig.1.E.1_", Sys.Date(), ".pdf", sep = ""), zoom1_plot, width = 8, height = 6, units = "cm")
ggsave(paste("./Fig.1.E.2_", Sys.Date(), ".pdf", sep = ""), zoom2_plot, width = 8, height = 6, units = "cm")
ggsave(paste("./Fig.1.E.3_", Sys.Date(), ".pdf", sep = ""), zoom3_plot, width = 8, height = 6, units = "cm")

##FigS1
figs1_color <- c("#97cd96", "#70bb6f", "#4da74a", "#3b8239", "#285a2f")

#s1 alter D
figs1_alterD_data_original <- cor_filter_P_data_200 %>% 
  ungroup() %>% 
  filter(Ed_cutoff %in% c(0.03, 0.045, 0.06, 0.075, 0.09), Alpha==1.5,El==5*10^-4, Em == 2e-3, P == 0) %>% 
  dplyr::select(1:8)

figs1_alterD_data_final <- tibble(Ed_cutoff = NA,Alpha=1.5,e0 = 0,El=5*10^-4, Em = 2e-3, P = 0, 
                                  seed_number = NA, rho = NA, p.value = NA)

for (i in 1 : nrow(figs1_alterD_data_original)) {
  ed <- figs1_alterD_data_original$Ed_cutoff[i]
  df <- figs1_alterD_data_original$data[i]
  rho <- cor.test(df[[1]]$Fake_T, df[[1]]$Fake_E, method = "spearman")$estimate
  p.value <- cor.test(df[[1]]$Fake_T, df[[1]]$Fake_E, method = "spearman")$p.value
  figs1_alterD_data_final <- rbind(figs1_alterD_data_final, tibble(Ed_cutoff = ed,Alpha=1.5,e0 = 0,El=5*10^-4, Em = 2e-3, P = 0, 
                                                                   seed_number = i, rho = rho, p.value = p.value)) %>% 
    filter(!is.na(seed_number))
  for (ratio in seq(0.2, 0.8, 0.2)) {
    new_data <- df[[1]] %>% filter(Fake_T >= quantile(Fake_T, ratio))
    new_ratio <- tibble(Ed_cutoff = ed,Alpha=1.5,e0 = 0,El=5*10^-4, Em = 2e-3, P = ratio, seed_number = i, 
                        rho = cor.test(new_data$Fake_T, new_data$Fake_E, method = "spearman")$estimate, 
                        p.value = cor.test(new_data$Fake_T, new_data$Fake_E, method = "spearman")$p.value)
    figs1_alterD_data_final <- rbind(figs1_alterD_data_final, new_ratio)
  }
}
minP <- min(filter(figs1_alterD_data_final, p.value != 0)$p.value)
figs1_alterD_data <- figs1_alterD_data_final %>% 
  mutate(p.value = ifelse(p.value == 0, minP, p.value)) %>% 
  mutate(adjust_rho = -rho*30) %>% 
  group_by(Ed_cutoff, P) %>%
  dplyr::summarise(mean_rho = mean(rho), sd_rho = sd(rho), 
                   adjust_mean_rho = mean(adjust_rho) + 2, adjust_sd_rho = sd(adjust_rho), 
                   mean_p.value = ifelse(mean(p.value)*5 > 1, 1, mean(p.value)*5), sd_p.value = sd(p.value)) %>%
  ungroup() %>% 
  mutate(mean_p.value = ifelse(mean_p.value < 1e-3, 1e-3, mean_p.value), adjust_mean_p.value = log10(mean_p.value)*6-2)

fig.s1.D <- figs1_alterD_data %>% ggplot()+
  geom_bar(aes(x=P,y=adjust_mean_rho, fill = as.factor(P)), 
           stat = "identity",position =position_dodge(width = 0.10),width = 0.10)+
  geom_bar(aes(x=P,y=adjust_mean_p.value, fill = as.factor(P)), 
           stat = "identity",position =position_dodge(width = 0.10),width = 0.01)+
  geom_errorbar(aes(x = P, ymin = adjust_mean_rho - adjust_sd_rho, ymax = adjust_mean_rho + adjust_sd_rho), 
                color = "black", stat = "identity",position = position_dodge(width = 0.1),width = 0.05) + 
  geom_point(aes(x=P,y=adjust_mean_p.value, fill = as.factor(P)), 
             position = position_dodge(width = 0.15), size=3, color = rep(figs1_color, 5), shape=21)+
  geom_abline(slope = 0,intercept=log10(0.05)*6-2,color="red",linewidth=0.8,linetype=2)+
  facet_grid(rows = vars(Ed_cutoff))+
  geom_text(aes(label=paste("italic(D)==",Ed_cutoff,sep = "")),x = 0.1,y = 20,parse = T, size = 8/.pt)+
  ###plot modify###
  scale_fill_manual(values = figs1_color)+
  scale_color_manual(values = figs1_color)+
  scale_x_continuous(limits = c(-0.1,0.9), breaks = seq(0, 0.8, 0.2), labels = c("0%", "20%", "40%", "60%", "80%"))+
  scale_y_continuous(limits = c(-22, 22), 
                     breaks=c(-18-2, -12-2, -6-2, -2, 2, 6+2, 12+2, 18+2), 
                     labels = c(expression(10^-3), expression(10^-2), 
                                expression(10^-1), 0, 
                                0, "-0.2", "-0.4", "-0.6"))+
  theme_classic()+theme(axis.text=element_text(size=9,face = "bold", color = "black"),
                        axis.title = element_blank(),
                        axis.text.x = element_text(),
                        axis.title.x = element_blank(),
                        axis.text.y = element_blank(),
                        strip.text = element_blank(),strip.background = element_blank(),
                        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
                        legend.position = "none"
  )

#s1 alter alpha
figs1_alterA_data_original <- cor_filter_P_data_200 %>% 
  ungroup() %>% 
  filter(Ed_cutoff == 0.06,El==5*10^-4, Em == 2e-3, P == 0, Alpha %in% c(1.25, 1.375, 1.5, 1.625, 1.75)) %>% 
  dplyr::select(1:8)

figs1_alterA_data_final <- tibble(Ed_cutoff = 0.06,Alpha=NA,e0 = 0,El=5*10^-4, Em = 2e-3, P = 0, 
                                  seed_number = NA, rho = NA, p.value = NA)

for (i in 1 : nrow(figs1_alterA_data_original)) {
  a <- figs1_alterA_data_original$Alpha[i]
  df <- figs1_alterA_data_original$data[i]
  rho <- cor.test(df[[1]]$Fake_T, df[[1]]$Fake_E, method = "spearman")$estimate
  p.value <- cor.test(df[[1]]$Fake_T, df[[1]]$Fake_E, method = "spearman")$p.value
  figs1_alterA_data_final <- rbind(figs1_alterA_data_final, tibble(Ed_cutoff = 0.06,Alpha=a,e0 = 0,El=5*10^-4, Em = 2e-3, P = 0, 
                                                                   seed_number = i, rho = rho, p.value = p.value)) %>% 
    filter(!is.na(seed_number))
  for (ratio in seq(0.2, 0.8, 0.2)) {
    new_data <- df[[1]] %>% filter(Fake_T >= quantile(Fake_T, ratio, na.rm = T))
    new_ratio <- tibble(Ed_cutoff = 0.06,Alpha=a,e0 = 0,El=5*10^-4, Em = 2e-3, P = ratio, seed_number = i, 
                        rho = cor.test(new_data$Fake_T, new_data$Fake_E, method = "spearman")$estimate, 
                        p.value = cor.test(new_data$Fake_T, new_data$Fake_E, method = "spearman")$p.value)
    figs1_alterA_data_final <- rbind(figs1_alterA_data_final, new_ratio)
  }
}

minP <- min(filter(figs1_alterA_data_final, p.value != 0)$p.value)
figs1_alterA_data <- figs1_alterA_data_final %>% 
  mutate(p.value = ifelse(p.value == 0, minP, p.value)) %>% 
  mutate(adjust_rho = -rho*30) %>% 
  group_by(Alpha, P) %>%
  dplyr::summarise(mean_rho = mean(rho), sd_rho = sd(rho), 
                   adjust_mean_rho = mean(adjust_rho) + 2, adjust_sd_rho = sd(adjust_rho), 
                   mean_p.value = ifelse(mean(p.value)*5 > 1, 1, mean(p.value)*5), sd_p.value = sd(p.value)) %>%
  ungroup() %>% 
  mutate(mean_p.value = ifelse(mean_p.value < 1e-3, 1e-3, mean_p.value), adjust_mean_p.value = log10(mean_p.value)*6-2)

fig.s1.A <- figs1_alterA_data %>% ggplot()+
  geom_bar(aes(x=P,y=adjust_mean_rho, fill = as.factor(P)), 
           stat = "identity",position =position_dodge(width = 0.10),width = 0.10)+
  geom_bar(aes(x=P,y=adjust_mean_p.value, fill = as.factor(P)), 
           stat = "identity",position =position_dodge(width = 0.10),width = 0.01)+
  geom_errorbar(aes(x = P, ymin = adjust_mean_rho - adjust_sd_rho, ymax = adjust_mean_rho + adjust_sd_rho), 
                color = "black", stat = "identity",position = position_dodge(width = 0.1),width = 0.05) + 
  geom_point(aes(x=P,y=adjust_mean_p.value, fill = as.factor(P)), 
             position = position_dodge(width = 0.15), size=3, color = rep(figs1_color, 5), shape=21)+
  geom_abline(slope = 0,intercept=log10(0.05)*6-2,color="red",linewidth=0.8,linetype=2)+
  facet_grid(rows = vars(Alpha))+
  geom_text(aes(label=paste("alpha==",Alpha,sep = "")),x = -0.05,y = 20,parse = T, hjust = 0, size = 8/.pt)+
  ###plot modify###
  scale_fill_manual(values = figs1_color)+
  scale_color_manual(values = figs1_color)+
  scale_x_continuous(limits = c(-0.1,0.9), breaks = seq(0, 0.8, 0.2), labels = c("0%", "20%", "40%", "60%", "80%"))+
  scale_y_continuous(limits = c(-22, 22), 
                     breaks=c(-18-2, -12-2, -6-2, -2, 2, 6+2, 12+2, 18+2), 
                     labels = c(expression(10^-3), expression(10^-2), 
                                expression(10^-1), 0, 
                                0, "-0.2", "-0.4", "-0.6"))+
  theme_classic()+theme(axis.text=element_text(size=9,face = "bold", color = "black"),
                        axis.text.x = element_text(),
                        axis.title.x = element_blank(),
                        axis.title.y = element_blank(),
                        axis.text.y = element_blank(),
                        strip.text = element_blank(),strip.background = element_blank(),
                        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
                        legend.position = "none"
  )

#s1 alter U
figs1_alterU_data_original <- cor_filter_P_data_200 %>% 
  ungroup() %>% 
  filter(Ed_cutoff == 0.06,Alpha == 1.5,El==5*10^-4, P == 0, Em %in% c(1e-3, 1.5e-3, 2e-3, 2.5e-3, 3e-3)) %>% 
  dplyr::select(1:8)

figs1_alterU_data_final <- tibble(Ed_cutoff = 0.06,Alpha=1.5,e0 = 0,El=5*10^-4, Em = NA, P = 0, 
                                  seed_number = NA, rho = NA, p.value = NA)

for (i in 1 : nrow(figs1_alterU_data_original)) {
  u <- figs1_alterU_data_original$Em[i]
  df <- figs1_alterU_data_original$data[i]
  rho <- cor.test(df[[1]]$Fake_T, df[[1]]$Fake_E, method = "spearman")$estimate
  p.value <- cor.test(df[[1]]$Fake_T, df[[1]]$Fake_E, method = "spearman")$p.value
  figs1_alterU_data_final <- rbind(figs1_alterU_data_final, tibble(Ed_cutoff = 0.06,Alpha=1.5,e0 = 0,El=5*10^-4, Em = u, P = 0, 
                                                                   seed_number = i, rho = rho, p.value = p.value)) %>% 
    filter(!is.na(seed_number))
  for (ratio in seq(0.2, 0.8, 0.2)) {
    new_data <- df[[1]] %>% filter(Fake_T >= quantile(Fake_T, ratio))
    new_ratio <- tibble(Ed_cutoff = 0.06,Alpha=1.5,e0 = 0,El=5*10^-4, Em = u, P = ratio, seed_number = i, 
                        rho = cor.test(new_data$Fake_T, new_data$Fake_E, method = "spearman")$estimate, 
                        p.value = cor.test(new_data$Fake_T, new_data$Fake_E, method = "spearman")$p.value)
    figs1_alterU_data_final <- rbind(figs1_alterU_data_final, new_ratio)
  }
}

minP <- min(filter(figs1_alterU_data_final, p.value != 0)$p.value)
figs1_alterU_data <- figs1_alterU_data_final %>% 
  mutate(p.value = ifelse(p.value == 0, minP, p.value)) %>% 
  mutate(adjust_rho = -rho*30) %>% 
  group_by(Em, P) %>%
  dplyr::summarise(mean_rho = mean(rho), sd_rho = sd(rho), 
                   adjust_mean_rho = mean(adjust_rho) + 2, adjust_sd_rho = sd(adjust_rho), 
                   mean_p.value = ifelse(mean(p.value)*5 > 1, 1, mean(p.value)*5), sd_p.value = sd(p.value)) %>%
  ungroup() %>% 
  mutate(mean_p.value = ifelse(mean_p.value < 1e-3, 1e-3, mean_p.value), adjust_mean_p.value = log10(mean_p.value)*6-2)

fig.s1.U <- figs1_alterU_data %>% ggplot()+
  geom_bar(aes(x=P,y=adjust_mean_rho, fill = as.factor(P)), 
           stat = "identity",position =position_dodge(width = 0.10),width = 0.10)+
  geom_bar(aes(x=P,y=adjust_mean_p.value, fill = as.factor(P)), 
           stat = "identity",position =position_dodge(width = 0.10),width = 0.01)+
  geom_errorbar(aes(x = P, ymin = adjust_mean_rho - adjust_sd_rho, ymax = adjust_mean_rho + adjust_sd_rho), 
                color = "black", stat = "identity",position = position_dodge(width = 0.1),width = 0.05) + 
  geom_point(aes(x=P,y=adjust_mean_p.value, fill = as.factor(P)), 
             position = position_dodge(width = 0.15), size=3, color = rep(figs1_color, 5), shape=21)+
  geom_abline(slope = 0,intercept=log10(0.05)*6-2,color="red",linewidth=0.8,linetype=2)+
  facet_grid(rows = vars(Em))+
  geom_text(aes(label=paste("italic(U)==",Em,sep = "")),x = -0.05,y = 20,parse = T, hjust = 0, size = 8/.pt)+
  ###plot modify###
  scale_fill_manual(values = figs1_color)+
  scale_color_manual(values = figs1_color)+
  scale_x_continuous(limits = c(-0.1,0.9), breaks = seq(0, 0.8, 0.2), labels = c("0%", "20%", "40%", "60%", "80%"))+
  scale_y_continuous(limits = c(-22, 22), 
                     breaks=c(-18-2, -12-2, -6-2, -2, 2, 6+2, 12+2, 18+2), 
                     labels = c(expression(10^-3), expression(10^-2), 
                                expression(10^-1), 0, 
                                0, "-0.2", "-0.4", "-0.6"))+
  theme_classic()+theme(axis.text=element_text(size=9,face = "bold", color = "black"),
                        axis.text.x = element_text(),
                        axis.title.x = element_blank(),
                        axis.title.y = element_blank(),
                        axis.text.y = element_blank(),
                        strip.text = element_blank(),strip.background = element_blank(),
                        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
                        legend.position = "none"
  )

#s1 alter L
figs1_alterL_data_original <- cor_filter_P_data_200 %>% 
  ungroup() %>% 
  filter(Ed_cutoff == 0.06,Alpha == 1.5,Em==0.002, P == 0, El %in% c(2.5e-4, 3.75e-4, 5e-4, 6.25e-4, 7.5e-4)) %>% 
  dplyr::select(1:8)

figs1_alterL_data_final <- tibble(Ed_cutoff = 0.06,Alpha=1.5,e0 = 0,El=NA, Em = 0.002, P = 0, 
                                  seed_number = NA, rho = NA, p.value = NA)

for (i in 1 : nrow(figs1_alterL_data_original)) {
  l <- figs1_alterL_data_original$El[i]
  df <- figs1_alterL_data_original$data[i]
  rho <- cor.test(df[[1]]$Fake_T, df[[1]]$Fake_E, method = "spearman")$estimate
  p.value <- cor.test(df[[1]]$Fake_T, df[[1]]$Fake_E, method = "spearman")$p.value
  figs1_alterL_data_final <- rbind(figs1_alterL_data_final, tibble(Ed_cutoff = 0.06,Alpha=1.5,e0 = 0,El=l, Em = 0.002, P = 0, 
                                                                   seed_number = i, rho = rho, p.value = p.value)) %>% 
    filter(!is.na(seed_number))
  for (ratio in seq(0.2, 0.8, 0.2)) {
    new_data <- df[[1]] %>% filter(Fake_T >= quantile(Fake_T, ratio, na.rm = T))
    new_ratio <- tibble(Ed_cutoff = 0.06,Alpha=1.5,e0 = 0,El=l, Em = 0.002, P = ratio, seed_number = i, 
                        rho = cor.test(new_data$Fake_T, new_data$Fake_E, method = "spearman")$estimate, 
                        p.value = cor.test(new_data$Fake_T, new_data$Fake_E, method = "spearman")$p.value)
    figs1_alterL_data_final <- rbind(figs1_alterL_data_final, new_ratio)
  }
}

minP <- min(filter(figs1_alterL_data_final, p.value != 0)$p.value)
figs1_alterL_data <- figs1_alterL_data_final %>% 
  mutate(p.value = ifelse(p.value == 0, minP, p.value)) %>% 
  mutate(adjust_rho = -rho*30) %>% 
  group_by(El, P) %>%
  dplyr::summarise(mean_rho = mean(rho), sd_rho = sd(rho), 
                   adjust_mean_rho = mean(adjust_rho) + 2, adjust_sd_rho = sd(adjust_rho), 
                   mean_p.value = ifelse(mean(p.value)*5 > 1, 1, mean(p.value)*5), sd_p.value = sd(p.value)) %>%
  ungroup() %>% 
  mutate(mean_p.value = ifelse(mean_p.value < 1e-3, 1e-3, mean_p.value), adjust_mean_p.value = log10(mean_p.value)*6-2)

fig.s1.L <- 
  figs1_alterL_data %>% ggplot()+
  geom_bar(aes(x=P,y=adjust_mean_rho, fill = as.factor(P)), 
           stat = "identity",position =position_dodge(width = 0.10),width = 0.10)+
  geom_bar(aes(x=P,y=adjust_mean_p.value, fill = as.factor(P)), 
           stat = "identity",position =position_dodge(width = 0.10),width = 0.01)+
  geom_errorbar(aes(x = P, ymin = adjust_mean_rho - adjust_sd_rho, ymax = adjust_mean_rho + adjust_sd_rho), 
                color = "black", stat = "identity",position = position_dodge(width = 0.1),width = 0.05) + 
  geom_point(aes(x=P,y=adjust_mean_p.value, fill = as.factor(P)), 
             position = position_dodge(width = 0.15), size=3, color = rep(figs1_color, 5), shape=21)+
  geom_abline(slope = 0,intercept=log10(0.05)*6-2,color="red",linewidth=0.8,linetype=2)+
  facet_grid(rows = vars(El))+
  geom_text(aes(label=paste("italic(L)==",El,sep = "")),x = -0.05,y = 20,parse = T, hjust = 0, size = 8/.pt)+
  ###plot modify###
  scale_fill_manual(values = figs1_color)+
  scale_color_manual(values = figs1_color)+
  scale_x_continuous(limits = c(-0.1,0.9), breaks = seq(0, 0.8, 0.2), labels = c("0%", "20%", "40%", "60%", "80%"))+
  scale_y_continuous(limits = c(-22, 22), 
                     breaks=c(-18-2, -12-2, -6-2, -2, 2, 6+2, 12+2, 18+2), 
                     labels = c(expression(10^-3), expression(10^-2), 
                                expression(10^-1), 0, 
                                0, "-0.2", "-0.4", "-0.6"))+
  theme_classic()+theme(axis.text=element_text(size=9,face = "bold", color = "black"),
                        axis.text.x = element_text(),
                        axis.title.x = element_blank(),
                        axis.title.y = element_blank(),
                        axis.text.y = element_blank(),
                        strip.text = element_blank(),strip.background = element_blank(),
                        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
                        legend.position = "none"
  )

#combine
fig.s1 <- cowplot::ggdraw()+
  cowplot::draw_plot(fig.s1.D,x=0, y=0, width = 0.25, height = 1)+
  cowplot::draw_plot(fig.s1.A,x=0.25, y=0, width = 0.25, height = 1)+
  cowplot::draw_plot(fig.s1.U,x=0.50, y=0, width = 0.25, height = 1)+
  cowplot::draw_plot(fig.s1.L,x=0.75, y=0, width = 0.25, height = 1)

ggsave(paste("./Fig.S1_", Sys.Date(), ".pdf", sep = ""), fig.s1, width = 17, height = 21, units = "cm")

