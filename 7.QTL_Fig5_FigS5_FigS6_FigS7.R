###figure 5a & figure 5b
# library
library(dplyr)
library(qtl)
library(ggbio)
library(reshape2)

# source function
Transform_strain_name <- function(strainname,name_type){
  strainname_matrix <- matrix(data = c(1:96),nrow = 8,ncol = 12,byrow = TRUE,
                              dimnames = list(c("A","B","C","D","E","F","G","H"),c(1:12)))
  s <- str_split(strainname,pattern = "-",simplify = TRUE)
  plate <- as.numeric(s[,1])
  RN <- str_extract(s[,2],pattern = "\\D")
  CN <- as.numeric(str_extract(s[,2],pattern = "\\d+"))
  strain_number <- c()
  for (i in 1 : length(RN)) {
    strain_number <- c(strain_number, strainname_matrix[RN[i], CN[i]])
  }
  output_strainname <- paste0("A",str_pad(plate,width = 2,pad = "0"),str_pad(strain_number,width = 2,pad = "0"))
  return(output_strainname)
}
# load data
#1# SI data#####
load("./3.SI_data.Rdata")

SI_804_data <- SI_rm_wt_gc_dp_dpsg_wm_ol15_list$SI_for_strain_clean_SI_biorep_ori_and_curve_rm_wt_gc_dp_dpsg_wm_ol15_data_output_list$SI_for_strain_clean_SI_biorep_ori_and_curve_rm_na_sd_day2_min_0.2_0.5_shift_0.3_day2

#2# TE data#####
TE_env <- new.env()
load("./5.TE_data.Rdata",envir = TE_env)
load("./5.TE_11_biorep_data.Rdata",envir = TE_env)
load("./5.TE_13_biorep_data.Rdata",envir = TE_env)

TE_11_for_strain_data <- TE_env$TE_11_biorep_clean_list$TE_11_for_strain_clean_bio_have_biorep_data %>%
  dplyr::rename(sd_F_R_ratio_11=sd_F_R_ratio,F_R_ratio_mean_bio_11=F_R_ratio_mean_bio)

TE_13_for_strain_data <- TE_env$TE_13_biorep_clean_list$TE_13_for_strain_clean_bio_data %>%
  filter(!is.na(sd_F_R_ratio))%>%
  dplyr::rename(sd_F_R_ratio_13=sd_F_R_ratio,F_R_ratio_mean_bio_13=F_R_ratio_mean_bio)

#356 observation
TE_271_and_11_13_data <- TE_env$TE_data_list$TE_for_strain_have_biorep_data %>%
  full_join(TE_11_for_strain_data) %>%
  full_join(TE_13_for_strain_data, by = c("strain")) %>%
  filter(!strain %in% c("BY4716","10-G12","RM11")) %>% #remove the wrong strain in the first batch and RM11 strain .
  filter(!strain %in% c("1-A10","4-C7","4-G7","4-E8","10-G12","1-G6")) %>% #remove the wrong strain with wrong genoytpe
  rowwise() %>%
  mutate(strain_rename=Transform_strain_name(strain))

#353 observation
TE_271_and_11_13_RM3point <- TE_271_and_11_13_data %>% 
  filter(!strain_rename %in% c("A1052","A1085","A1044"))

#345 observation
TE_271_and_11_13_RM11point <- TE_271_and_11_13_RM3point %>%
  filter(!strain_rename %in% c("A0141","A1035","A1029","A1024")) %>%
  filter(!strain_rename %in% c("A0190","A0645","A0863","A0306")) 

#4# cross data############
load("./7.cross_data.Rdata")

# MAIN#######
#1.Construct TE_SI_data (use TE_data left join SI_data,remove NA) and the subset 0.5 fraction data######
TE_SI_data <- TE_271_and_11_13_RM11point %>%
  left_join(SI_804_data,by = c("strain_rename"="strains")) %>%
  na.omit() %>%
  ungroup()

TE_SI_data_SI_0.5 <- TE_SI_data %>%
  arrange(SI_mean_bio) %>%
  filter(row_number() >= dim(TE_SI_data)[1] * 0.5)

#2.Create cross object for following qtl mapping######
cross_TE_SI_data_SI_0.5 <- cross
cross_TE_SI_data_SI_0.5$pheno <- cross$pheno[,1:4] %>% left_join(TE_SI_data_SI_0.5,by = c("strains"="strain_rename"))

#3.qtl mapping in the subset data########
out.mr.cross_TE_SI_data_SI_0.5 <- scanone(cross_TE_SI_data_SI_0.5,pheno.col = c("SI_mean_bio","TE_mean_bio"),method = "mr")
#permutation result
#opermuate.mr.cross_TE_SI_data_SI_0.5 <- scanone(cross_TE_SI_data_SI_0.5,pheno.col = c("SI_mean_bio","TE_mean_bio"),
#                                                method = "mr",n.perm = 1000)
#> summary(opermuate.mr.cross_TE_SI_data_SI_0.5,alpha = c(0.05,0.1,0.2))
#LOD thresholds (1000 permutations)
#SI_mean_bio TE_mean_bio
#5%         3.59        3.49
#10%        3.28        3.23
#20%        2.99        2.89

#Plot####
#@###1.Figure5 A#####
FIG5.A_For_plot_ob <- out.mr.cross_TE_SI_data_SI_0.5 %>%
  rownames_to_column(var = "Marker") %>%
  rownames_to_column(var = "Marker_id") %>%
  mutate(Marker_id = as.numeric(Marker_id)) %>%
  melt(id.vars=c("Marker","Marker_id","chr","pos"))

#@#####1.1 chr rect####
chr_marker_id_length <- out.mr.cross_TE_SI_data_SI_0.5  %>% group_by(chr) %>% dplyr::mutate(chr_N=n()) %>% 
  filter(row_number()==1) %>% ungroup %>% mutate(v0=cumsum(chr_N)) %>% pull(v0) %>%c(0,.)
chr_marker_id_mean <- mclapply(1:length(chr_marker_id_length[-1]),
                               FUN = function(x){
                                 o <- mean(c(chr_marker_id_length[x],chr_marker_id_length[x+1]))
                                 return(o)}) %>%unlist()

names(chr_marker_id_mean) <-  rownames(out.mr.cross_TE_SI_data_SI_0.5) %>% gsub(".*_(chr\\D+)_.*","\\1",.) %>% unique()
chr_rect <- data.frame(x1=chr_marker_id_length[-length(chr_marker_id_length)],x2=chr_marker_id_length[-1],
                       y1=rep(-2,length(chr_marker_id_length[-1])),y2=rep(5,length(chr_marker_id_length[-1])),
                       c=(rep(c("#BEBEBE50",NA),8)),stringsAsFactors = F)
#@#####1.2 Lod value point plot in whole genome####
FIG5.A <- 
  FIG5.A_For_plot_ob %>%
  ggplot()+
  geom_rect(data = chr_rect[seq(1,15,by=2),],mapping = aes(xmin = x1,xmax=x2, fill=c), 
            ymin = -Inf,ymax = Inf,show.legend = F)+
  geom_point(aes(x=Marker_id,y=value,color=variable), size = 0.1)+
  geom_hline(yintercept = c(3.59,3.49),color = c("#4daf4a","#2075bc"),linetype = 2)+
  labs(x="Genome position",y="LOD score")+
  scale_x_continuous(breaks = chr_marker_id_mean)+
  scale_y_continuous(limits = c(0,4.5))+
  scale_fill_manual(values=c("#BEBEBE"))+
  scale_color_manual(name="",values=c("#4daf4a","#2075bc"),label=c("Lifespan","Translational error rate (e)"))+
  guides(color=guide_legend(override.aes = list(size=3)))+
  theme_classic()+
  theme(#plot.title = element_text(hjust = 0.5,size = 16, face = "bold"),
    axis.text=element_text(size=9,face = "bold",color = "black"),
    axis.title = element_text(size=10,face = "bold"),
    axis.text.x = element_text(angle = 20),
    panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
    legend.text = element_text(size=8,face = "bold"),legend.position = "none")

#@###2.Figure5 B#####
Interval_rect <- data.frame( x1=641753,x2=669427,y1=-2,y2=5,c="#BEBEBE50")

FIG5.B1_change <-
  FIG5.A_For_plot_ob %>%
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
  theme(#plot.title = element_text(hjust = 0.5,size = 16, face = "bold"),
    axis.text=element_text(size=9,face = "bold", color = "black"),
    axis.title = element_text(size=10,face = "bold"),
    axis.text.x = element_text(),
    panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
    legend.text = element_text(size=8,face = "bold"),legend.position = "none")

# save
ggsave(paste("./Fig.5.A_", Sys.Date(), ".pdf", sep = ""), FIG5.A, width = 13+1, height = 6, units = "cm")
ggsave(paste("./Fig.5.B_", Sys.Date(), ".pdf", sep = ""), FIG5.B1_change, width = 9, height = 6-0.5, units = "cm")


###figure 5c & figure 5d
# library
library(qtl)

# source function
Transform_strain_name <- function(strainname,name_type){
  strainname_matrix <- matrix(data = c(1:96),nrow = 8,ncol = 12,byrow = TRUE,
                              dimnames = list(c("A","B","C","D","E","F","G","H"),c(1:12)))
  s <- str_split(strainname,pattern = "-",simplify = TRUE)
  plate <- as.numeric(s[,1])
  RN <- str_extract(s[,2],pattern = "\\D")
  CN <- as.numeric(str_extract(s[,2],pattern = "\\d+"))
  strain_number <- c()
  for (i in 1 : length(RN)) {
    strain_number <- c(strain_number, strainname_matrix[RN[i], CN[i]])
  }
  output_strainname <- paste0("A",str_pad(plate,width = 2,pad = "0"),str_pad(strain_number,width = 2,pad = "0"))
  return(output_strainname)
}

PlotPXG_ggplot <- function(cross,marker,pheno.col,ylim=c(),ybre=c(),ylab = c(),ylabs="Survival integral"){
  PlotPXG_object <- qtl::plotPXG(x = cross,marker = marker,pheno.col = pheno.col,infer = F)
  SUMMARY_PlotPXG <- PlotPXG_object %>% 
    na.omit() %>%
    group_by(!!sym(marker)) %>%
    dplyr::summarise(m=mean(pheno),se=plotrix::std.error(pheno))
  
  output <- PlotPXG_object %>%
    na.omit() %>%
    ggplot(group=!!sym(marker))+
    geom_point(aes(x=!!sym(marker),y=pheno,color = !!sym(marker)),position = position_jitter(w = 0.2, h = 0), 
               size = 0.3)+
    geom_errorbar(aes(ymin=SUMMARY_PlotPXG$m[1]-SUMMARY_PlotPXG$se[1],ymax=SUMMARY_PlotPXG$m[1]+SUMMARY_PlotPXG$se[1],x=1),
                  size=0.3,width=0.2)+
    geom_errorbar(aes(ymin=SUMMARY_PlotPXG$m[2]-SUMMARY_PlotPXG$se[2],ymax=SUMMARY_PlotPXG$m[2]+SUMMARY_PlotPXG$se[2],x=2),
                  size=0.3,width=0.2)+
    geom_segment(aes(y = SUMMARY_PlotPXG$m[1], yend = SUMMARY_PlotPXG$m[1], x = 0.75, xend = 1.25), size = 0.3)+
    geom_segment(aes(y = SUMMARY_PlotPXG$m[2], yend = SUMMARY_PlotPXG$m[2], x = 1.75, xend = 2.25), size = 0.3)+
    scale_x_discrete(labels= c("BY allele","RM allele"))+
    scale_color_manual(values = c("#e5d851", "#8c07bb"))+
    scale_y_continuous(limits=ylim, breaks = ybre, labels = ylab)+
    labs(title = marker,x="Genotype",y=ylabs)+
    theme_classic()+theme(#plot.title = element_text(hjust = 0.5,size = 11, face = "bold"),
      axis.text=element_text(size=9,face = "bold",color = "black"),
      axis.title = element_text(size=10,face = "bold",color = "black"),
      axis.text.x = element_text(vjust = 1),
      axis.title.x = element_blank(),plot.title = element_blank(),
      panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
      legend.position = "none")
  return(output)
}

# load data
#1# SI data#####
load("./3.SI_data.Rdata")

SI_804_data <- SI_rm_wt_gc_dp_dpsg_wm_ol15_list$SI_for_strain_clean_SI_biorep_ori_and_curve_rm_wt_gc_dp_dpsg_wm_ol15_data_output_list$SI_for_strain_clean_SI_biorep_ori_and_curve_rm_na_sd_day2_min_0.2_0.5_shift_0.3_day2

#2# TE data#####
TE_env <- new.env()
load("./5.TE_data.Rdata",envir = TE_env)
load("./5.TE_11_biorep_data.Rdata",envir = TE_env)
load("./5.TE_13_biorep_data.Rdata",envir = TE_env)

TE_11_for_strain_data <- TE_env$TE_11_biorep_clean_list$TE_11_for_strain_clean_bio_have_biorep_data %>%
  dplyr::rename(sd_F_R_ratio_11=sd_F_R_ratio,F_R_ratio_mean_bio_11=F_R_ratio_mean_bio)

TE_13_for_strain_data <- TE_env$TE_13_biorep_clean_list$TE_13_for_strain_clean_bio_data %>%
  filter(!is.na(sd_F_R_ratio))%>%
  dplyr::rename(sd_F_R_ratio_13=sd_F_R_ratio,F_R_ratio_mean_bio_13=F_R_ratio_mean_bio)

#356 observation
TE_271_and_11_13_data <- TE_env$TE_data_list$TE_for_strain_have_biorep_data %>%
  full_join(TE_11_for_strain_data) %>%
  full_join(TE_13_for_strain_data, by = c("strain")) %>%
  filter(!strain %in% c("BY4716","10-G12","RM11")) %>% #remove the wrong strain in the first batch and RM11 strain .
  filter(!strain %in% c("1-A10","4-C7","4-G7","4-E8","10-G12","1-G6")) %>% #remove the wrong strain with wrong genoytpe
  rowwise() %>%
  mutate(strain_rename=Transform_strain_name(strain))

#353 observation
TE_271_and_11_13_RM3point <- TE_271_and_11_13_data %>% 
  filter(!strain_rename %in% c("A1052","A1085","A1044"))

#345 observation
TE_271_and_11_13_RM11point <- TE_271_and_11_13_RM3point %>%
  filter(!strain_rename %in% c("A0141","A1035","A1029","A1024")) %>%
  filter(!strain_rename %in% c("A0190","A0645","A0863","A0306")) 

#4# cross data############
load("./7.cross_data.Rdata")

# MAIN#######
#1.Construct TE_SI_data (use TE_data left join SI_data,remove NA) and the subset 0.5 fraction data######
TE_SI_data <- TE_271_and_11_13_RM11point %>%
  left_join(SI_804_data,by = c("strain_rename"="strains")) %>%
  na.omit() %>%
  ungroup()

TE_SI_data_SI_0.5 <- TE_SI_data %>%
  arrange(SI_mean_bio) %>%
  filter(row_number() >= dim(TE_SI_data)[1] * 0.5)

#2.Create cross object for following qtl mapping######
cross_TE_SI_data_SI_0.5 <- cross
cross_TE_SI_data_SI_0.5$pheno <- cross$pheno[,1:4] %>% left_join(TE_SI_data_SI_0.5,by = c("strains"="strain_rename"))

#3.qtl mapping in the subset data########
out.mr.cross_TE_SI_data_SI_0.5 <- scanone(cross_TE_SI_data_SI_0.5,pheno.col = c("SI_mean_bio","TE_mean_bio"),method = "mr")

summary(out.mr.cross_TE_SI_data_SI_0.5,format = "allpheno")

#4.choose one site for plot #####
#PLOT########
FIG5.C <- 
  PlotPXG_ggplot(cross = cross_TE_SI_data_SI_0.5,marker = "6487936_chrX_655475_C_T",pheno.col = "SI_mean_bio",
                 ylim = c(6.5,13),ybre = seq(7, 13, 2),ylab = seq(7, 13, 2),ylabs = "Lifespan")
FIG5.D <- 
  PlotPXG_ggplot(cross = cross_TE_SI_data_SI_0.5,marker = "6487936_chrX_655475_C_T",pheno.col = "TE_mean_bio",
                 ylim = c(0.00065,0.0013),ybre = seq(7e-4, 13e-4, 2e-4),ylab = seq(7, 13, 2),
                 ylabs = "Translational error rate")

# save
ggsave(paste("./Fig.5.C_", Sys.Date(), ".pdf", sep = ""), FIG5.C, width = 4.5, height = 6-1.5, units = "cm")
ggsave(paste("./Fig.5.D_", Sys.Date(), ".pdf", sep = ""), FIG5.D, width = 4.5, height = 6-1.5, units = "cm")


###figure 5e & figure 5f
# load data
GeneVerify_SI_TE_env <- new.env()
load("./7.SI_gene_verification_data.Rdata",envir = GeneVerify_SI_TE_env)
load("./7.TE_gene_verification_data.Rdata",envir = GeneVerify_SI_TE_env)

# MAIN#######
#1.SI gene verify data in chr 10#####
GeneVerify_SI_10_chr_data <- GeneVerify_SI_TE_env$SI_tuberep_clean_list$SI_tuberep_clean_day2_min_list$SI_for_biorep_clean_SI_tuberep_ori_and_curve_day2_min_data %>%
  filter(strains %in% c("STE24","ILM1","JHD2","YJR124C","VPS70","chrX658789")) %>%
  filter(strains != "YJR124C")

BY_RM_data_correspond_to_GeneVerify_SI_10_chr_data <- GeneVerify_SI_TE_env$SI_tuberep_clean_list$SI_tuberep_clean_day2_min_list$SI_for_biorep_clean_SI_tuberep_ori_and_curve_day2_min_data %>%
  filter(strains %in% c("BY4716","RM11")) %>%
  filter(strains_group %in% GeneVerify_SI_10_chr_data$strains_group[-c(1:3)])

GeneVerify_SI_10_chr_data_PLUS_BY_RM <- bind_rows(GeneVerify_SI_10_chr_data,
                                                  BY_RM_data_correspond_to_GeneVerify_SI_10_chr_data)

#2.TE gene verify data in chr 10#####
GeneVerify_TE_VPS_BY_data <- GeneVerify_SI_TE_env$TE_GeneVerify_data_list$TE_for_biorep_data %>% filter(strain %in% c("BY4716","VPS70","RM11"))

#PLOT#####
#@###1.SI GeneVerify biorep######
#@######1.1 SI gene sort#############
SI_Geneverify_gene_name_sort <- c("BY"=1,"ILM1"=2,"JHD2"=3,"STE24"=4,"VPS70"=5,"chrX658789"=6,"RM"=13)

SI_Geneverify_gene_comparisons_list <- mclapply(1:length(SI_Geneverify_gene_name_sort[-1]),FUN = function(x){
  out <- c("BY",names(SI_Geneverify_gene_name_sort[-1][x]))
  return(out)
})  

FIG5.E <- 
  GeneVerify_SI_10_chr_data_PLUS_BY_RM %>% 
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
  theme_classic()+theme(#plot.title = element_text(hjust = 0.5,size = 16, face = "bold"),
    axis.text=element_text(size=9,face = "bold", color = "black"),
    axis.title = element_text(size=10,face = "bold", color = "black"),
    #axis.text.x = element_text(angle = 20,vjust = 1,hjust = 1),
    axis.text.x = element_blank(),
    panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
    legend.text = element_text(size=15,face = "bold"),legend.title = element_text(size=15,face = "bold"))

#@###2.TE GeneVerify biorep######
#@######2.1 TE gene sort#############
TE_Geneverify_gene_name_sort <- c("BY"=1,"ADP1"=2,"POL4"=3,"CTO1"=4,"YCR016W"=5,"CWH43"=6,"YCR018C"=7,
                                  "MAK32"=8,"HSP30"=9,"YCR022C"=10,"YCR023C"=11,"VPS70"=12,"RM"=13)

TE_Geneverify_gene_name_sort_subset <- TE_Geneverify_gene_name_sort[c(1,12,13)]

TE_Geneverify_gene_comparisons_list <- mclapply(1:length(TE_Geneverify_gene_name_sort_subset[-1]),FUN = function(x){
  out <- c("BY",names(TE_Geneverify_gene_name_sort_subset[-1][x]))
  return(out)
})  

data <- GeneVerify_TE_VPS_BY_data%>% 
  mutate(strain = ifelse(strain %in% c("BY4716"),"BY",ifelse(strain %in% c("RM11"),"RM",strain))) %>%
  dplyr::mutate(strain_sort = TE_Geneverify_gene_name_sort[[strain]]) %>%
  mutate(strain_col = case_when(strain %in% c("BY") ~ 1,
                                strain %in% c("RM") ~ 99,
                                strain %in% c("VPS70") ~ 3,
                                TRUE ~ 2))

FIG5.F <- 
  data %>% ggplot(aes(x=reorder(strain,strain_sort),y=TE))+
  geom_boxplot(aes(fill = as.factor(strain_col)),show.legend = F, size = 0.4)+
  stat_compare_means(aes(label = ..p.format..),comparisons = TE_Geneverify_gene_comparisons_list,
                     method = "wilcox.test",method.args = list(alternative = "greater"), size = 8/.pt)+
  scale_fill_manual(name="Strain",values=c("#e5d851","#ffb6c1","#8c07bb"))+
  scale_y_continuous(breaks = c(seq(0.0006, 0.0015, 0.0003)), labels = c(seq(6, 15, 3)))+
  scale_x_discrete(labels=c('BY', 'BY::VPS70-RM','RM'))+
  labs(x="",y="Translational error rate")+
  theme_classic()+theme(#plot.title = element_text(hjust = 0.5,size = 16, face = "bold"),
    axis.text=element_text(size=9,face = "bold", color = "black"),
    axis.title = element_text(size=10,face = "bold", color = "black"),
    #axis.text.x = element_text(angle = 20,vjust = 1,hjust = 1),
    axis.text.x = element_blank(),
    panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
    legend.text = element_text(size=15,face = "bold"),legend.title = element_text(size=15,face = "bold"))

# save
ggsave(paste("./Fig.5.E_", Sys.Date(), ".pdf", sep = ""), FIG5.E, width = 8, height = 6-1, units = "cm")
ggsave(paste("./Fig.5.F_", Sys.Date(), ".pdf", sep = ""), FIG5.F, width = 4, height = 6-0.5, units = "cm")


###figure 5g
###library
library(dplyr)
library(ggplot2)

###loading
load("./7.VPS70_data_treated_with_ConA.Rdata")

###main
#translational error rate
mean_by_te <- group_by(TransError, treatment, strain) %>% 
  dplyr::summarise(mean = mean(TransError, na.rm = T)) %>% 
  ungroup() %>% filter(strain == "BY")

fig.5g_data <- filter(TransError, !is.na(TransError)) %>% 
  merge(mean_by_te, by = "treatment") %>% 
  mutate(Strain = strain.x, 
         normal_TE = TransError / mean, 
         treatment = gsub("DMSO", "1", treatment))

cona_by <- fig.5g_data %>% filter(treatment == "ConA", Strain == "BY")
cona_vps70 <- fig.5g_data %>% filter(treatment == "ConA", Strain == "BY::VPS70-RM")
t.test(cona_by$normal_TE, cona_vps70$normal_TE, alternative = "greater") #t = 0.032906, df = 8.1365, p-value = 0.4873

dmso_by <- fig.5g_data %>% filter(treatment == "1", Strain == "BY")
dmso_vps70 <- fig.5g_data %>% filter(treatment == "1", Strain == "BY::VPS70-RM")
t.test(dmso_by$TransError, dmso_vps70$TransError, alternative = "greater") #t = 3.2378, df = 6.4106, p-value = 0.008079

fig.5g <- 
  ggplot(fig.5g_data, aes(Strain, normal_TE, fill = Strain)) + 
  geom_boxplot(size = 0.4) + 
  geom_jitter(width = 0.1, size = 0.8) + 
  facet_wrap(~treatment, labeller = as_labeller(c(`1` = "DMSO", ConA = "ConA"))) + 
  scale_fill_manual(values = c("#dcd154", "#ffb6c1")) + 
  scale_y_continuous(limits = c(0.66, 1.3), breaks = c(0.7, 0.9, 1.1, 1.3))+
  #ylab("Translational error rate (e) (x10-4)")+
  theme_classic()+
  theme(#plot.title = element_text(hjust = 0.5,size = 18, face = "bold"),
    axis.text= element_text(size=9,face = "bold", color = "black"),
    axis.title = element_text(size=10,face = "bold", color = "black"),
    axis.text.x = element_blank(), axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
    legend.background = element_blank(),
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_blank(),
    strip.text = element_text(size=8,color = "black",face = "bold"))

#lifespan
mean_by_si <- group_by(SurvivalIntegral_finale, treatment, strain) %>% 
  dplyr::summarise(mean = mean(SI, na.rm = T)) %>% 
  ungroup() %>% filter(strain == "BY")

fig.5g2_data <- merge(SurvivalIntegral_finale, mean_by_si, by = "treatment") %>% 
  mutate(Strain = strain.x, 
         normal_SI = SI / mean, 
         treatment = gsub("DMSO", "1", treatment))

cona_by <- fig.5g2_data %>% filter(treatment == "ConA", Strain == "BY")
cona_vps70 <- fig.5g2_data %>% filter(treatment == "ConA", Strain == "VPS70")
t.test(cona_by$SI, cona_vps70$SI, alternative = "less") #t = -2.2697, df = 6.3252, p-value = 0.03073

dmso_by <- fig.5g2_data %>% filter(treatment == "1", Strain == "BY")
dmso_vps70 <- fig.5g2_data %>% filter(treatment == "1", Strain == "VPS70")
t.test(dmso_by$SI, dmso_vps70$SI, alternative = "less") #t = -6.5178, df = 7.8185, p-value = 0.0001021

fig.5g2 <- 
  ggplot(fig.5g2_data, aes(Strain, normal_SI, fill = Strain)) + 
  geom_boxplot(size = 0.4) + 
  geom_jitter(width = 0.1, size = 0.8) + 
  facet_wrap(~treatment, labeller = as_labeller(c(`1` = "DMSO", ConA = "ConA"))) + 
  scale_fill_manual(values = c("#dcd154", "#ffb6c1")) + 
  scale_y_continuous(limits = c(0.83, 1.7), breaks = seq(0.9, 1.7, 0.2))+
  #ylab("Translational error rate (e) (x10-4)")+
  theme_classic()+
  theme(#plot.title = element_text(hjust = 0.5,size = 18, face = "bold"),
    axis.text= element_text(size=9,face = "bold", color = "black"),
    axis.title = element_text(size=10,face = "bold", color = "black"),
    axis.text.x = element_blank(), axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
    legend.background = element_blank(),
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_blank(),
    strip.text = element_text(size=8,color = "black",face = "bold"))

###save
ggsave(paste("./Fig.5g_", Sys.Date(), ".pdf", sep = ""), fig.5g, width = 6-0.5, height = 6-0.5, units = "cm")
ggsave(paste("./Fig.5g2_", Sys.Date(), ".pdf", sep = ""), fig.5g2, width = 6-0.5, height = 6-0.5, units = "cm")


###figure S5
# load data
GeneVerify_SI_TE_env <- new.env()
load("./7.SI_gene_verification_data.Rdata",envir = GeneVerify_SI_TE_env)
load("./7.TE_gene_verification_data.Rdata",envir = GeneVerify_SI_TE_env)

# MAIN#######
GeneVerify_TE_chr3_data <-  GeneVerify_SI_TE_env$TE_GeneVerify_data_list$TE_for_biorep_data %>%
  filter(!strain %in% c("MAK32","POL4","YCR018C","VPS70"))

#ADP1,POL4,YCR016W,CWH43,YCR018C,MAK32,HSP30,YCR022C,YCR023C#
TE_Geneverify_gene_name_sort <- c("BY"=1,"ADP1"=2,"POL4"=3,"CTO1"=4,"YCR016W"=5,"CWH43"=6,"YCR018C"=7,
                                  "MAK32"=8,"HSP30"=9,"YCR022C"=10,"YCR023C"=11,"VPS70"=12,"RM"=13)

TE_Geneverify_gene_name_sort_subset <- TE_Geneverify_gene_name_sort[c(1:2,4:6,9:11,13)]

TE_Geneverify_gene_comparisons_list <- mclapply(1:length(TE_Geneverify_gene_name_sort_subset[-1]),FUN = function(x){
  out <- c("BY",names(TE_Geneverify_gene_name_sort_subset[-1][x]))
  return(out)
})  

#plot
fig.S5 <- 
  GeneVerify_TE_chr3_data %>% 
  mutate(strain = ifelse(strain %in% c("BY4716"),"BY",ifelse(strain %in% c("RM11"),"RM",strain))) %>%
  mutate(strain_sort = TE_Geneverify_gene_name_sort[[strain]]) %>%
  mutate(strain_col = case_when(strain %in% c("BY") ~ 1,
                                strain %in% c("RM") ~ 99,
                                strain %in% c("CWH43") ~ 3,
                                TRUE ~ 2)) %>%
  ggplot(aes(x=reorder(strain,strain_sort),y=TE))+
  geom_boxplot(aes(fill = as.factor(strain_col)),show.legend = F)+
  stat_compare_means(aes(label = ..p.format..),comparisons = TE_Geneverify_gene_comparisons_list,
                     method = "wilcox.test",method.args = list(alternative = "less"), size = 8/.pt)+
  scale_fill_manual(name="Strain",values=c("#ffdb00", "white","#e72019","#ff6d00"))+
  scale_x_discrete(labels=c('BY', 'BY::ADP1-RM','BY::CTO1-RM','BY::YCR016W-RM','BY::CWH43-RM',
                            'BY::HSP30-RM','BY::YCR022C-RM','BY::YCR023C-RM','RM'))+
  scale_y_continuous(breaks = c(seq(8e-4, 20e-4, 4e-4)), labels = c(seq(8, 20, 4)), limits = c(7e-4, 21e-4))+
  labs(x="",y="Translation error rate ")+
  theme_classic()+theme(#plot.title = element_text(hjust = 0.5,size = 16, face = "bold"),
    axis.text=element_text(size=9,face = "bold", color = "black"),
    axis.title = element_text(size=10,face = "bold", color = "black"),
    axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
    #axis.text.x = element_blank(),
    panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
    legend.text = element_text(size=15,face = "bold"),legend.title = element_text(size=15,face = "bold"))

#save
saveDir <- "./"
ggsave(paste(saveDir, "fig.S5_", Sys.Date(), ".pdf", sep = ""), fig.S5, width = 8, height = 10, units = "cm")


###figure S6
###library
library(dplyr)
library(ggplot2)

###loading
load("./7.VPS70_SI_data_treated_with_rapamycin.Rdata")
load("./7.VPS70_TE_data_treated_with_rapamycin.Rdata")

###main
####Normalize_lifespand####
#####DMSO_Rapamycin#####
SI_plotdata <- SI_data$SI_DMSO_Rapamycin

mean_by_SI <- group_by(SI_plotdata, condition, strain) %>% 
  dplyr::summarise(mean = mean(SI, na.rm = T)) %>% 
  ungroup() %>% filter(strain == "BY4716")

fig.s6_data <- SI_plotdata %>% 
  merge(mean_by_SI, by = "condition") %>% 
  mutate(Strain = strain.x, 
         normal_SI = SI / mean, 
         condition = gsub("DMSO", "1", condition))

rapa_by <- fig.s6_data %>% filter(condition == "Rapamycin", Strain == "BY4716")
rapa_vps70 <- fig.s6_data %>% filter(condition == "Rapamycin", Strain == "VPS70")
t.test(rapa_by $normal_SI, rapa_vps70$normal_SI, alternative = "less") #t = -0.053666, df = 12.049, p-value = 0.479

dmso_by <- fig.s6_data %>% filter(condition == "1", Strain == "BY4716")
dmso_vps70 <- fig.s6_data %>% filter(condition == "1", Strain == "VPS70")
t.test(dmso_by$SI, dmso_vps70$SI, alternative = "less") #t = -3.5602, df = 12.732, p-value = 0.001797

fig.s6_SI <- ggplot(fig.s6_data, aes(Strain, normal_SI, fill = Strain)) + 
  geom_boxplot(size = 0.4) + 
  geom_jitter(width = 0.1, size = 0.8) + 
  facet_wrap(~condition, labeller = as_labeller(c(`1` = "DMSO", Rapamycin = "Rapamycin"))) + 
  scale_fill_manual(values = c("#dcd154", "#ffb6c1")) + 
  scale_y_continuous(limits = c(0.6, 2.4), breaks = c(0.6, 1, 1.4,1.8,2.2))+
  #ylab("Translational error rate (e) (x10-4)")+
  theme_classic()+
  theme(#plot.title = element_text(hjust = 0.5,size = 18, face = "bold"),
    axis.text= element_text(size=9,face = "bold", color = "black"),
    axis.title = element_text(size=10,face = "bold", color = "black"),
    axis.text.x = element_blank(), axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
    legend.background = element_blank(),
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_blank(),
    strip.text = element_text(size=8,color = "black",face = "bold"), 
    legend.position = "none")
fig.s6_SI

#save
ggsave(paste("./fig_s6_SI_", Sys.Date(), ".pdf", sep = ""), fig.s6_SI, width = 8-0.5, height = 6-0.5, units = "cm")

####Normalize_TransError####
#####DMSO_Rapamycin#####
TE_plotdata <- TE_data$TE_DMSO_rapamycin

mean_by_te <- group_by(TE_plotdata, condition, strain) %>% 
  dplyr::summarise(mean = mean(TE, na.rm = T)) %>% 
  ungroup() %>% filter(strain == "BY")

fig.s62_data <- merge(TE_plotdata, mean_by_te, by = "condition") %>% 
  mutate(Strain = strain.x, 
         normal_TE = TE / mean, 
         condition = gsub("DMSO", "1", condition))

rapa_by <- fig.s62_data %>% filter(condition == "Rapamycin", Strain == "BY")
rapa_vps70 <- fig.s62_data %>% filter(condition == "Rapamycin", Strain == "VPS70")
t.test(rapa_by$TE, rapa_vps70$TE, alternative = "greater") #t = 0.10278, df = 18.125, p-value = 0.4596

dmso_by <- fig.s62_data %>% filter(condition == "1", Strain == "BY")
dmso_vps70 <- fig.s62_data %>% filter(condition == "1", Strain == "VPS70")
t.test(dmso_by$TE, dmso_vps70$TE, alternative = "greater") #t = 1.9578, df = 16.468, p-value = 0.03371

fig.s6_TE <- ggplot(fig.s62_data, aes(Strain, normal_TE, fill = Strain)) + 
  geom_boxplot(size = 0.4) + 
  geom_jitter(width = 0.1, size = 0.8) + 
  facet_wrap(~condition, labeller = as_labeller(c(`1` = "DMSO", Rapamycin = "Rapamycin"))) + 
  scale_fill_manual(values = c("#dcd154", "#ffb6c1")) + 
  scale_y_continuous(limits = c(0.4, 1.7), breaks = seq(0.4, 1.6, 0.3))+
  #ylab("Translational error rate (e) (x10-4)")+
  theme_classic()+
  theme(#plot.title = element_text(hjust = 0.5,size = 18, face = "bold"),
    axis.text= element_text(size=9,face = "bold", color = "black"),
    axis.title = element_text(size=10,face = "bold", color = "black"),
    axis.text.x = element_blank(), axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
    legend.background = element_blank(),
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_blank(),
    strip.text = element_text(size=8,color = "black",face = "bold"), 
    legend.position = "none")
fig.s6_TE

#save
ggsave(paste("./fig.s6_TE_", Sys.Date(), ".pdf", sep = ""), fig.s6_TE, width = 8-0.5, height = 6-0.5, units = "cm")


###figure S7a
# source function
Transform_strain_name <- function(strainname,name_type){
  strainname_matrix <- matrix(data = c(1:96),nrow = 8,ncol = 12,byrow = TRUE,dimnames = list(c("A","B","C","D","E","F","G","H"),
                                                                                             c(1:12)))
  s <- str_split(strainname,pattern = "-",simplify = TRUE)
  plate <- as.numeric(s[,1])
  RN <- str_extract(s[,2],pattern = "\\D")
  CN <- as.numeric(str_extract(s[,2],pattern = "\\d+"))
  strain_number <- c()
  for (i in 1 : length(RN)) {
    strain_number <- c(strain_number, strainname_matrix[RN[i], CN[i]])
  }
  output_strainname <- paste0("A",str_pad(plate,width = 2,pad = "0"),str_pad(strain_number,width = 2,pad = "0"))
  return(output_strainname)
}

# load data
#1.SI data ######
load("./3.SI_data.Rdata")

#2.cross data
load("./7.cross_data.Rdata")

#3.TE data 
TE_env <- new.env()
load("./5.TE_data.Rdata",envir = TE_env)
load("./5.TE_11_biorep_data.Rdata",envir = TE_env)
load("./5.TE_13_biorep_data.Rdata",envir = TE_env)

TE_11_for_strain_data <- TE_env$TE_11_biorep_clean_list$TE_11_for_strain_clean_bio_have_biorep_data %>%
  dplyr::rename(sd_F_R_ratio_11=sd_F_R_ratio,F_R_ratio_mean_bio_11=F_R_ratio_mean_bio)

TE_13_for_strain_data <- TE_env$TE_13_biorep_clean_list$TE_13_for_strain_clean_bio_data %>%
  filter(!is.na(sd_F_R_ratio))%>%
  dplyr::rename(sd_F_R_ratio_13=sd_F_R_ratio,F_R_ratio_mean_bio_13=F_R_ratio_mean_bio)

TE_271_and_11_13_data <- TE_env$TE_data_list$TE_for_strain_have_biorep_data %>%
  full_join(TE_11_for_strain_data) %>%
  full_join(TE_13_for_strain_data, by = c("strain")) %>%
  filter(!strain %in% c("BY4716","10-G12","RM11")) %>% #remove the wrong strain in the first batch and RM11 strain .
  filter(!strain %in% c("1-A10","4-C7","4-G7","4-E8","10-G12","1-G6")) %>% #remove the wrong strain with wrong genoytpe
  rowwise() %>%
  mutate(strain_rename=Transform_strain_name(strain))

te_data <- TE_271_and_11_13_data %>% dplyr::select(strain_rename, TE_mean_bio, sd_TE) %>% 
  filter(!is.na(TE_mean_bio)) %>%
  filter(!strain_rename %in% c("BY4716","RM11","A0141","A1035","A1029",
                               "A1024","A0190","A0645","A0863","A0306",
                               "A1052","A1085","A1044",
                               "A0110","A0431","A1084","A0479","A0456"))

# MAIN#######
#@###1.QTL mapping in SI #######
#@######1.1 create SI_804_cross####
SI_804_cross <- cross
SI_804_cross$pheno <- cross$pheno[,c(1:4)] %>%  
  left_join(SI_rm_wt_gc_dp_dpsg_wm_ol15_list$SI_for_strain_clean_SI_biorep_ori_and_curve_rm_wt_gc_dp_dpsg_wm_ol15_data_output_list$SI_for_strain_clean_SI_biorep_ori_and_curve_rm_na_sd_day2_min_0.2_0.5_shift_0.3_day2)

all_data <- left_join(SI_rm_wt_gc_dp_dpsg_wm_ol15_list$SI_for_strain_clean_SI_biorep_ori_and_curve_rm_wt_gc_dp_dpsg_wm_ol15_data_output_list$SI_for_strain_clean_SI_biorep_ori_and_curve_rm_na_sd_day2_min_0.2_0.5_shift_0.3_day2, 
                      te_data, by = c("strains"="strain_rename"))
all_cross <- cross
all_cross$pheno <- cross$pheno[,c(1:4)] %>%  left_join(all_data, by = "strains")
out.mr.all_data <- scanone(all_cross, method = "mr", pheno.col = c("SI_mean_bio", "TE_mean_bio"))
#out.mr.all_threshold <- scanone(all_cross, method = "mr", pheno.col = c("SI_mean_bio", "TE_mean_bio"), n.perm = 100)
#summary(out.mr.all_threshold, alpha = c(0.05, 0.1, 0.2))
#LOD thresholds (100 permutations)
#SI_mean_bio TE_mean_bio
#5%         3.30        3.51
#10%        3.15        3.37
#20%        2.95        2.86

#@######1.2  SI_804 qtl mapping####
out.mr.SI_804 <- scanone(SI_804_cross,method = "mr",pheno.col = c("SI_mean_bio"))

# outpermutation.mr.SI_804 <- scanone(SI_804_cross,method = "mr",pheno.col = c("SI_mean_bio"),n.perm = 1000)
# > summary(outpermutation.mr.SI_804,alpha = c(0.05,0.1,0.2))
# LOD thresholds (1000 permutations)
# lod
# 5%  3.59
# 10% 3.25
# 20% 2.95

#Plot####
#@###1.FigureS7 A#####
FIGS7.A_For_plot_ob <- out.mr.all_data %>%
  rownames_to_column(var = "Marker") %>%
  rownames_to_column(var = "Marker_id") %>%
  mutate(Marker_id = as.numeric(Marker_id)) %>%
  melt(id.vars=c("Marker","Marker_id","chr","pos"))

#@#####1.1 chr rect####
chr_marker_id_length <- FIGS7.A_For_plot_ob  %>% 
  group_by(chr) %>% dplyr::mutate(chr_N=n()) %>% filter(row_number()==1) %>% 
  ungroup %>% mutate(v0=cumsum(chr_N)) %>% pull(v0) %>%c(0,.)
chr_marker_id_mean <- mclapply(1:length(chr_marker_id_length[-1]),
                               FUN = function(x){
                                 o <- mean(c(chr_marker_id_length[x],chr_marker_id_length[x+1]))
                                 return(o)}) %>%unlist()

names(chr_marker_id_mean) <-  rownames(out.mr.SI_804) %>% 
  gsub(".*_(chr\\D+)_.*","\\1",.) %>% unique()
chr_rect <- data.frame(x1=chr_marker_id_length[-length(chr_marker_id_length)]/2,x2=chr_marker_id_length[-1]/2,
                       y1=rep(-2,length(chr_marker_id_length[-1])),y2=rep(5,length(chr_marker_id_length[-1])),
                       c=(rep(c("#BEBEBE50",NA),8)),stringsAsFactors = F)
#@#####1.2 Lod value point plot in whole genome####

FIGS7.A <- 
  FIGS7.A_For_plot_ob %>% 
  ggplot()+
  geom_rect(data = chr_rect[seq(1,15,by=2),],mapping = aes(xmin = x1,xmax=x2, fill=c), 
            ymin = -Inf,ymax = Inf,show.legend = F)+
  geom_point(aes(x=Marker_id,y=value,color=variable), size = 0.1)+
  geom_hline(yintercept = c(3.30),color = c("#4daf4a"),linetype = 2)+
  geom_hline(yintercept = c(3.51),color = c("#2075bc"),linetype = 2)+
  labs(x="Genome position",y="LOD score")+
  scale_x_continuous(breaks = chr_marker_id_mean/2)+
  #scale_y_continuous(limits = c(0,5))+
  scale_fill_manual(values=c("#BEBEBE"))+
  scale_color_manual(name="",values=c("#4daf4a","#2075bc"),label=c("Lifespan","Translation error rate"))+
  theme_classic()+
  theme(#plot.title = element_text(hjust = 0.5,size = 18, face = "bold"),
    axis.text=element_text(size=9,face = "bold",color = "black"),
    axis.title = element_text(size=10,face = "bold"),
    axis.text.x = element_text(angle = 20),
    panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
    legend.text = element_text(size=12,face = "bold"),legend.position = "none")

# save
ggsave(paste("./Fig.S7.A_", Sys.Date(), ".pdf", sep = ""), FIGS7.A, width = 10.8+0.75, height = 6, units = "cm")


###figure S7c & figure S7d
# source function
PlotPXG_ggplot <- function(cross,marker,pheno.col,ylim=c(),ybre=c(),ylab = c(),ylabs="Survival integral"){
  PlotPXG_object <- qtl::plotPXG(x = cross,marker = marker,pheno.col = pheno.col,infer = F)
  SUMMARY_PlotPXG <- PlotPXG_object %>% 
    na.omit() %>%
    group_by(!!sym(marker)) %>%
    dplyr::summarise(m=mean(pheno),se=plotrix::std.error(pheno))
  
  output <- PlotPXG_object %>%
    na.omit() %>%
    ggplot(group=!!sym(marker))+
    geom_point(aes(x=!!sym(marker),y=pheno,color = !!sym(marker)),position = position_jitter(w = 0.2, h = 0), 
               size = 0.3)+
    geom_errorbar(aes(ymin=SUMMARY_PlotPXG$m[1]-SUMMARY_PlotPXG$se[1],ymax=SUMMARY_PlotPXG$m[1]+SUMMARY_PlotPXG$se[1],x=1),
                  size=0.3,width=0.2)+
    geom_errorbar(aes(ymin=SUMMARY_PlotPXG$m[2]-SUMMARY_PlotPXG$se[2],ymax=SUMMARY_PlotPXG$m[2]+SUMMARY_PlotPXG$se[2],x=2),
                  size=0.3,width=0.2)+
    geom_segment(aes(y = SUMMARY_PlotPXG$m[1], yend = SUMMARY_PlotPXG$m[1], x = 0.75, xend = 1.25), size = 0.3)+
    geom_segment(aes(y = SUMMARY_PlotPXG$m[2], yend = SUMMARY_PlotPXG$m[2], x = 1.75, xend = 2.25), size = 0.3)+
    scale_x_discrete(labels= c("BY","RM"))+
    scale_color_manual(values = c("#e5d851", "#8c07bb"))+
    scale_y_continuous(limits=ylim, breaks = ybre, labels = ylab)+
    labs(title = marker,x="Genotype",y=ylabs)+
    theme_classic()+theme(#plot.title = element_text(hjust = 0.5,size = 11, face = "bold"),
      axis.text=element_text(size=9,face = "bold",color = "black"),
      axis.title = element_text(size=10,face = "bold",color = "black"),
      axis.text.x = element_text(vjust = 1),
      axis.title.x = element_blank(),plot.title = element_blank(),
      panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
      legend.position = "none")
  return(output)
}

# load data
# MAIN#######
#1.SI data ######
load("./3.SI_data.Rdata")

#2.cross data
load("./7.cross_data.Rdata")

# MAIN#######
#@###1.QTL mapping in SI #######
#@######1.1 create SI_804_cross####
SI_804_cross <- cross
SI_804_cross$pheno <- cross$pheno[,c(1:4)] %>%  
  left_join(SI_rm_wt_gc_dp_dpsg_wm_ol15_list$SI_for_strain_clean_SI_biorep_ori_and_curve_rm_wt_gc_dp_dpsg_wm_ol15_data_output_list$SI_for_strain_clean_SI_biorep_ori_and_curve_rm_na_sd_day2_min_0.2_0.5_shift_0.3_day2)

#@######1.2  SI_804 qtl mapping####
out.mr.SI_804 <- scanone(SI_804_cross,method = "mr",pheno.col = c("SI_mean_bio"))

#PLOT########
FIGS7.C <- 
  PlotPXG_ggplot(cross = SI_804_cross,marker = "7905399_chrXII_660371_C_T",pheno.col = "SI_mean_bio",
                 ylim = c(1.9, 13.9), ybre = seq(3, 12, 3), ylab = seq(3, 12, 3), ylabs = "Lifespan")
FIGS7.D <- 
  PlotPXG_ggplot(cross = SI_804_cross,marker = "9709121_chrXIV_461485_A_G",pheno.col = "SI_mean_bio",
                 ylim = c(1.9, 13.9), ybre = seq(3, 12, 3), ylab = seq(3, 12, 3), ylabs = "Lifespan")

# save
ggsave(paste("./Fig.S7.C_", Sys.Date(), ".pdf", sep = ""), FIGS7.C, width = 3.6, height = 6-1.25, units = "cm")
ggsave(paste("./Fig.S7.D_", Sys.Date(), ".pdf", sep = ""), FIGS7.D, width = 3.6, height = 6-1.25, units = "cm")


###figure S7f & figure S7b &figure S7e
#library
library(stringr)

# source function
Transform_strain_name <- function(strainname,name_type){
  strainname_matrix <- matrix(data = c(1:96),nrow = 8,ncol = 12,byrow = TRUE,
                              dimnames = list(c("A","B","C","D","E","F","G","H"),c(1:12)))
  s <- str_split(strainname,pattern = "-",simplify = TRUE)
  plate <- as.numeric(s[,1])
  RN <- str_extract(s[,2],pattern = "\\D")
  CN <- as.numeric(str_extract(s[,2],pattern = "\\d+"))
  strain_number <- c()
  for (i in 1 : length(RN)) {
    strain_number <- c(strain_number, strainname_matrix[RN[i], CN[i]])
  }
  output_strainname <- paste0("A",str_pad(plate,width = 2,pad = "0"),str_pad(strain_number,width = 2,pad = "0"))
  return(output_strainname)
}

# load data
#1.SI data output list 
load("./3.SI_data.Rdata")

#2.TE data 
TE_env <- new.env()
load("./5.TE_data.Rdata",envir = TE_env)
load("./5.TE_11_biorep_data.Rdata",envir = TE_env)
load("./5.TE_13_biorep_data.Rdata",envir = TE_env)

TE_11_for_strain_data <- TE_env$TE_11_biorep_clean_list$TE_11_for_strain_clean_bio_have_biorep_data %>%
  dplyr::rename(sd_F_R_ratio_11=sd_F_R_ratio,F_R_ratio_mean_bio_11=F_R_ratio_mean_bio)

TE_13_for_strain_data <- TE_env$TE_13_biorep_clean_list$TE_13_for_strain_clean_bio_data %>%
  filter(!is.na(sd_F_R_ratio))%>%
  dplyr::rename(sd_F_R_ratio_13=sd_F_R_ratio,F_R_ratio_mean_bio_13=F_R_ratio_mean_bio)

#356 observation
TE_271_and_11_13_data <- TE_env$TE_data_list$TE_for_strain_have_biorep_data %>%
  full_join(TE_11_for_strain_data) %>%
  full_join(TE_13_for_strain_data, by = c("strain")) %>%
  filter(!strain %in% c("BY4716","10-G12","RM11")) %>% #remove the wrong strain in the first batch and RM11 strain .
  filter(!strain %in% c("1-A10","4-C7","4-G7","4-E8","10-G12","1-G6")) %>% #remove the wrong strain with wrong genoytpe
  rowwise() %>%
  mutate(strain_rename=Transform_strain_name(strain))

#353 observation
TE_271_and_11_13_RM3point <- TE_271_and_11_13_data %>% 
  filter(!strain_rename %in% c("A1052","A1085","A1044"))

#345 observation
TE_271_and_11_13_RM11point <- TE_271_and_11_13_RM3point %>%
  filter(!strain_rename %in% c("A0141","A1035","A1029","A1024")) %>%
  filter(!strain_rename %in% c("A0190","A0645","A0863","A0306")) 

#4.Genotype VerMT data
Genotype_VerMT_env <- new.env()
load("./7.genotype_data.Rdata",envir = Genotype_VerMT_env)

#5.cross data
load("./7.cross_data.Rdata")

# MAIN#######
TE_SI_SNP_data <- TE_271_and_11_13_RM11point %>%
  ungroup%>%
  arrange(-TE_mean_bio) %>%
  # filter(row_number()>=5)%>%
  left_join(Genotype_VerMT_env$genotype_verMT_FULL_GENOTYPE_for_qtl_mapping_CSV_data%>%
              .[-1,] %>%
              dplyr::select(strains,`6490171_chrX_657710_C_T`,`7905399_chrXII_660371_C_T`,`9709121_chrXIV_461485_A_G`,
                            `1195948_chrIII_152546_A_G`),by = c("strain_rename"="strains")) %>%
  # filter(!is.na(`7905399_chrXII_660371_C_T`)) %>%
  left_join(SI_rm_wt_gc_dp_dpsg_wm_ol15_list$SI_for_strain_clean_SI_biorep_ori_and_curve_rm_wt_gc_dp_dpsg_wm_ol15_data_output_list$SI_for_strain_clean_SI_biorep_ori_and_curve_rm_na_sd_day2_min_0.2_0.5_shift_0.3_day2,
            by=c("strain_rename"="strains")) 

TE_SI_SNP_data_No2.highest <- TE_SI_SNP_data %>%
  na.omit() %>%
  filter(`7905399_chrXII_660371_C_T`==1)

TE_SI_SNP_data_No3.highest <- TE_SI_SNP_data %>%
  na.omit() %>%
  filter(`9709121_chrXIV_461485_A_G`==1)

TE_SI_SNP_data_No2and3.highest <- TE_SI_SNP_data %>%
  na.omit() %>%
  filter(`7905399_chrXII_660371_C_T`==1,`9709121_chrXIV_461485_A_G`==1)

TE_SI_SNP_data_No2and3.highest_cross <- cross
TE_SI_SNP_data_No2and3.highest_cross$pheno <- cross$pheno[,c(1:4)] %>% 
  left_join(TE_SI_SNP_data_No2and3.highest,by = c("strains"="strain_rename"))

out.mr.TE_SI_SNP_data_No2and3.highest <- scanone(TE_SI_SNP_data_No2and3.highest_cross,
                                                 pheno.col = c("SI_mean_bio","TE_mean_bio"),method = "mr")
#outpermutation.mr.TE_SI_SNP_data_No2and3.highest <- scanone(TE_SI_SNP_data_No2and3.highest_cross,
#                                                            pheno.col = c("SI_mean_bio","TE_mean_bio"),
#                                                            method = "mr", n.perm = 200)
#> summary(outpermutation.mr.TE_SI_SNP_data_No2and3.highest)
#LOD thresholds (200 permutations)
#SI_mean_bio TE_mean_bio
#5%         4.00        4.04
#10%        3.46        3.31

#te and si
outpermutation.mr.TE_SI_SNP_data_No2and3.highest <- out.mr.TE_SI_SNP_data_No2and3.highest

FIGS7.F <- outpermutation.mr.TE_SI_SNP_data_No2and3.highest %>%
  rownames_to_column(var = "Marker") %>%
  rownames_to_column(var = "Marker_id") %>%
  mutate(Marker_id = as.numeric(Marker_id)) %>%
  melt(id.vars=c("Marker","Marker_id","chr","pos"))

chr_marker_id_length <- outpermutation.mr.TE_SI_SNP_data_No2and3.highest  %>% 
  group_by(chr) %>% dplyr::mutate(chr_N=n()) %>% filter(row_number()==1) %>% 
  ungroup %>% mutate(v0=cumsum(chr_N)) %>% pull(v0) %>%c(0,.)
chr_marker_id_mean <- mclapply(1:length(chr_marker_id_length[-1]),
                               FUN = function(x){
                                 o <- mean(c(chr_marker_id_length[x],chr_marker_id_length[x+1]))
                                 return(o)}) %>%unlist()

names(chr_marker_id_mean) <-  rownames(outpermutation.mr.TE_SI_SNP_data_No2and3.highest) %>% 
  gsub(".*_(chr\\D+)_.*","\\1",.) %>% unique()
chr_rect <- data.frame(x1=chr_marker_id_length[-length(chr_marker_id_length)],x2=chr_marker_id_length[-1],
                       y1=rep(-2,length(chr_marker_id_length[-1])),y2=rep(5,length(chr_marker_id_length[-1])),
                       c=(rep(c("#BEBEBE50",NA),8)),stringsAsFactors = F)

fig.s7f.plot <- 
  FIGS7.F %>% ggplot()+
  geom_rect(data = chr_rect[seq(1,15,by=2),],mapping = aes(xmin = x1,xmax=x2, fill=c), ymin = -Inf,ymax = Inf,show.legend = F)+
  geom_point(aes(x=Marker_id,y=value,color=variable), size = 0.1)+
  #geom_hline(yintercept = c(3.59,3.49),color = c("#4daf4a","#2075bc"),linetype = 2)+
  labs(x="Genome position",y="LOD score")+
  scale_x_continuous(breaks = chr_marker_id_mean)+
  scale_y_continuous(limits = c(0,3.5))+
  scale_fill_manual(values=c("#BEBEBE"))+
  scale_color_manual(name="",values=c("#4daf4a","#2075bc"),label=c("Lifespan","Translational error rate (e)"))+
  guides(color=guide_legend(override.aes = list(size=3)))+
  theme_classic()+
  theme(#plot.title = element_text(hjust = 0.5,size = 16, face = "bold"),
    axis.text=element_text(size=9,face = "bold",color = "black"),
    axis.title = element_text(size=10,face = "bold"),
    axis.text.x = element_text(angle = 20),
    panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
    legend.text = element_text(size=8,face = "bold"),legend.position = "none")

#figure S7b & figure S7e
#PLOT######
FigureS7.B <- 
  TE_SI_SNP_data %>%
  filter(!is.na(SI_mean_bio)) %>% 
  filter(!is.na(TE_mean_bio)) %>% 
  ggplot()+
  geom_point(aes(x=TE_mean_bio,y=SI_mean_bio), size = 0.2)+
  #stat_cor(aes(x=TE_mean_bio,y=SI_mean_bio),method = "spearman",cor.coef.name = "rho", size = 8/.pt)+
  #rho = -0.032, P = 0.63
  geom_smooth(aes(x=TE_mean_bio,y=SI_mean_bio), formula = "y ~ x",method = "lm",se = F, color = "#4daf4a")+
  labs(x="Translational error rate (e) (x10-4)",y="Lifespan")+
  scale_x_continuous(limit = c(6.6e-4, 12.6e-4), breaks = seq(7e-4, 11e-4, 2e-4), labels = seq(7, 11, 2))+
  scale_y_continuous(limits = c(1.9, 14), breaks = seq(3, 12, 3), labels = seq(3, 12, 3)) + 
  theme_classic()+theme(#plot.title = element_text(hjust = 0.5,size = 16, face = "bold"),
    axis.text=element_text(size=9,face = "bold", color = "black"),
    axis.title = element_text(size=10,face = "bold", color = "black"),
    axis.title.x = element_blank(),axis.title.y = element_blank(),
    panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
    legend.position = "none")

FigureS7.E <- 
  TE_SI_SNP_data_No2and3.highest %>%
  ggplot()+
  geom_point(aes(x=TE_mean_bio,y=SI_mean_bio), size = 0.2)+
  #stat_cor(aes(x=TE_mean_bio,y=SI_mean_bio),method = "spearman",cor.coef.name = "rho", size = 8/.pt)+
  #rho = -0.4, P = 0.002
  geom_smooth(aes(x=TE_mean_bio,y=SI_mean_bio),method = "lm",se = F, color = "#4daf4a")+
  labs(x="Translational error rate (e) (x10-4)",y="Lifespan")+
  scale_x_continuous(limit = c(6.6e-4, 12.6e-4), breaks = seq(7e-4, 11e-4, 2e-4), labels = seq(7, 11, 2))+
  scale_y_continuous(limits = c(1.9, 14), breaks = seq(3, 12, 3), labels = seq(3, 12, 3)) + 
  theme_classic()+theme(#plot.title = element_text(hjust = 0.5,size = 16, face = "bold"),
    axis.text=element_text(size=9,face = "bold", color = "black"),
    axis.title = element_text(size=10,face = "bold", color = "black"),
    axis.title.x = element_blank(),axis.title.y = element_blank(),
    panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
    legend.position = "none")

# save
ggsave(paste("./Fig.S7.F_", Sys.Date(), ".pdf", sep = ""), fig.s7f.plot, width = 10.8+0.75, height = 6, units = "cm")
ggsave(paste("./Fig.S7.B_", Sys.Date(), ".pdf", sep = ""), FigureS7.B, width = 3.6-0.5, height = 6-0.75, units = "cm")
ggsave(paste("./Fig.S7.E_", Sys.Date(), ".pdf", sep = ""), FigureS7.E, width = 3.6-0.5, height = 6-0.75, units = "cm")

