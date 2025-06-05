###figure 6
#library
library(dplyr)
#library(genio)
library(lfa)
library(gcatest)
library(ggplot2)
library(ggbreak)

#load
lifespan_all <- read.table("./8.lifespan_data.txt", header = T)
sex_all <- read.table("./8.sex_data.txt", header = T)

#system("~/src/plink_1.9/plink --bfile ./chr11 --extract ./snps.txt --geno 0.1 --hwe 0.0001 --maf 0.01 --make-bed --mind 0.1 --keep ./lifespan_data.txt --out ./chr11_2_genes")
#data <- read_plink("./chr11_2_genes")
#save(data_0.1, file = "./8.chr11_2_genes.Rdata")
load("./8.chr11_2_genes.Rdata")

#system("~/src/plink_1.9/plink --bfile ./chr11 --extract ./snps.txt --geno 0.1 --hwe 0.0001 --maf 0.01 --make-bed --mind 0.1 --keep ./lifespan_data_0.1.txt --out ./chr11_2_genes_0.1")
#data_0.1 <- read_plink("./chr11_2_genes_0.1")
#save(data_0.1, file = "./8.chr11_2_genes_0.1.Rdata")
load("./8.chr11_2_genes_0.1.Rdata")

#genotype-conditional association test
set.seed(111)
geno <- data$X
trait <- lifespan_all %>% filter(id %in% data$fam$id) %>% select(lifespan)
trait <- trait$lifespan
cov <- sex_all %>% filter(id %in% data$fam$id) %>% select(sex) %>% as.matrix()
LF <- lfa(geno, 3)
gof <- sHWE(geno, LF, B=2)
filtered <- gof < (1 / nrow(geno))
geno <- geno[!filtered,]
gcat_p <- gcat(geno, LF, trait, cov)
gcat_p

set.seed(111)
geno_0.1 <- data_0.1$X
trait_0.1 <- lifespan_all %>% filter(id %in% data_0.1$fam$id) %>% select(lifespan)
trait_0.1 <- trait_0.1$lifespan
cov_0.1 <- sex_all %>% filter(id %in% data_0.1$fam$id) %>% select(sex) %>% as.matrix()
LF_0.1 <- lfa(geno_0.1, 3)
gof_0.1 <- sHWE(geno_0.1, LF_0.1, B=2)
filtered_0.1 <- gof_0.1 < (1 / nrow(geno_0.1))
geno_0.1 <- geno_0.1[!filtered_0.1,]
gcat_p_0.1 <- gcat(geno_0.1, LF_0.1, trait_0.1, cov_0.1)
gcat_p_0.1

#figure 6a & figure 6b
mySNP <- tibble(id = names(gcat_p), p = gcat_p) %>% 
  mutate(bon_p = ifelse(p*length(gcat_p)*2 > 1, 1, p*length(gcat_p)*2)) %>% 
  left_join(data$bim, by = "id") %>% 
  mutate(gene = ifelse(pos > 5e7, "NAALAD2", "FOLH1"))
mySNP_0.1 <- tibble(id = names(gcat_p_0.1), p = gcat_p_0.1) %>% 
  mutate(bon_p = ifelse(p*length(gcat_p_0.1)*2 > 1, 1, p*length(gcat_p_0.1)*2)) %>% 
  left_join(data_0.1$bim, by = "id") %>% 
  mutate(gene = ifelse(pos > 5e7, "NAALAD2", "FOLH1"))

exons_NAALAD2 <- c(89867927, 89868008, 89868727, 89868838, 89880498, 89880684, 89882174, 89882275, 89883650, 89883775, 
                   89885466, 89885652, 89891313, 89891406, 89892407, 89892505, 89896117, 89896202, 89896478, 89896597, 
                   89896703, 89896785, 89902097, 89902160, 89903237, 89903304, 89906992, 89907083, 89909140, 89909230, 
                   89911021, 89911285, 89914788, 89914869, 89916084, 89916176, 89924726, 89924915)
exons_FOLH1 <- c(49168308, 49168497, 49170191, 49170283, 49175398, 49175479, 49175780, 49176044, 49178269, 49178359, 
                 49179504, 49179595, 49186257, 49186324, 49190747, 49190810, 49192747, 49192829, 49194909, 49195028, 
                 49196444, 49196529, 49197411, 49197509, 49204701, 49204794, 49207221, 49207407, 49208196, 49208321, 
                 49214345, 49214446, 49221807, 49221993, 49227619, 49227724, 49229844, 49229961)

fig6a_plot <- ggplot(mySNP, aes(pos, bon_p)) + 
  annotate("rect", xmin = 4916e04, xmax = 8993e04, ymin = 1, ymax = 1.05, fill = "#666666") + 
  annotate("rect", xmin = 4916e04, xmax = 8993e04, ymin = 0.5, ymax = 1, fill = "#858585") + 
  annotate("rect", xmin = 4916e04, xmax = 8993e04, ymin = 0.1, ymax = 0.5, fill = "#a3a3a3") + 
  annotate("rect", xmin = 4916e04, xmax = 8993e04, ymin = 0.05, ymax = 0.1, fill = "#c2c2c2") + 
  annotate("rect", xmin = 4916e04, xmax = 8993e04, ymin = 0, ymax = 0.05, fill = "#e0e0e0") + 
  geom_point(aes(color = gene), size = 2) + 
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
  fig6a_plot <- fig6a_plot +
    annotate("rect", xmin = exons_NAALAD2[i*2-1], xmax = exons_NAALAD2[i*2], ymin = -0.01, ymax = -0.09, fill = "#134e13")+
    annotate("rect", xmin = exons_FOLH1[i*2-1], xmax = exons_FOLH1[i*2], ymin = -0.01, ymax = -0.09, fill = "#134e13")
}
fig6a_plot

fig6b_plot <- ggplot(mySNP_0.1, aes(pos, bon_p)) + 
  annotate("rect", xmin = 4916e04, xmax = 8993e04, ymin = 1, ymax = 1.05, fill = "#666666") + 
  annotate("rect", xmin = 4916e04, xmax = 8993e04, ymin = 0.5, ymax = 1, fill = "#858585") + 
  annotate("rect", xmin = 4916e04, xmax = 8993e04, ymin = 0.1, ymax = 0.5, fill = "#a3a3a3") + 
  annotate("rect", xmin = 4916e04, xmax = 8993e04, ymin = 0.05, ymax = 0.1, fill = "#c2c2c2") + 
  annotate("rect", xmin = 4916e04, xmax = 8993e04, ymin = 0, ymax = 0.05, fill = "#e0e0e0") + 
  geom_point(aes(color = gene), size = 2) + 
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
  fig6b_plot <- fig6b_plot +
    annotate("rect", xmin = exons_NAALAD2[i*2-1], xmax = exons_NAALAD2[i*2], ymin = -0.01, ymax = -0.09, fill = "#134e13")+
    annotate("rect", xmin = exons_FOLH1[i*2-1], xmax = exons_FOLH1[i*2], ymin = -0.01, ymax = -0.09, fill = "#134e13")
}
fig6b_plot

#figure 6c
fig6c_data <- tibble(allele = geno_0.1["rs80078229",], lifespan = trait_0.1) %>% 
  mutate(allele = ifelse(allele == 0, "AA", "AG"))
lifespan1 <- fig6c_data %>% filter(allele == "AG")
lifespan2 <- fig6c_data %>% filter(allele == "AA")
sum(lifespan1$lifespan)/nrow(lifespan1) #79.9
sum(lifespan2$lifespan)/nrow(lifespan2) #78.95918
wilcox.test(lifespan1$lifespan, lifespan2$lifespan) #W = 393.5, p-value = 0.01992

fig6c_plot <- ggplot(fig6c_data, aes(allele, lifespan)) + 
  geom_boxplot(aes(fill = allele), size = 0.5, color = "black") +
  scale_fill_manual(values = c("#8ecf8c", "#357833"))+
  labs(x="Allele",y="Lifespan") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_blank(), 
        axis.title = element_blank(), 
        legend.position = "none")
fig6c_plot

#save
ggsave(paste("./Fig.6A_", Sys.Date(), ".pdf", sep = ""), fig6a_plot, width = 9+1, height = 5+1, units = "cm")
ggsave(paste("./Fig.6B_", Sys.Date(), ".pdf", sep = ""), fig6b_plot, width = 9+1, height = 5+1, units = "cm")
ggsave(paste("./Fig.6C_", Sys.Date(), ".pdf", sep = ""), fig6c_plot, width = 5.9+0.1, height = 10.8+0.2, units = "cm")

