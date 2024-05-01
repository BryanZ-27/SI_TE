#author gongwang yu
######################################
library(ggplot2)
library(dplyr)
library(readr)
library(parallel)
library(tidyr)
library(plyr)
library(ggpubr)
#for fig.1e
dfqsubs<-read.csv("detct_quant.csv")

dfqsubs<-mclapply(1:nrow(dfqsubs), function(i){
  dfBP<-dfqsubs[i,]%>%dplyr::select(c(starts_with("Base.intensity")))%>%t() %>%
    as.data.frame()%>%
    tibble::rownames_to_column()%>%
    tidyr::separate(col = "rowname",into = c("a","sample"),sep=".Exp..")%>%
    dplyr::select(2,3)
  colnames(dfBP)[2]<-"bp"

  dfDP<-dfqsubs[i,]%>%dplyr::select(c(starts_with("Modification.intensity")))%>%t()%>%
    as.data.frame()%>%
    tibble::rownames_to_column()%>%
    tidyr::separate(col = "rowname",into = c("a","sample"),sep=".Exp..")%>%
    dplyr::select(2,3)
  colnames(dfDP)[2]<-"dp"
  
  dfmerg<-merge(dfBP,dfDP)%>%dplyr::rowwise()%>%dplyr::mutate(dfratio=dp/sum(dp,bp,na.rm = T)) %>%
    dplyr::filter(dfratio<=0.01)
  
  young<-dfmerg%>%dplyr::filter(grepl("young",sample,ignore.case=T))%>%
    getElement("dfratio")%>%mean(na.rm=T)
  old<-dfmerg%>%dplyr::filter(grepl("old",sample,ignore.case=T))%>%
    getElement("dfratio")%>%mean(na.rm=T)
  
  return(data.frame(young_rate=young,old_rate=old))
  
},mc.cores = 50)%>%rbind.fill()%>%
  filter(young_rate >0 | old_rate >0)%>%na.omit()

misincorp <- dfqsubs %>% 
  dplyr::rename(young=young_rate,old=old_rate) %>% 
  mutate(old_increased = (old - young)/(old + young) ) %>%
  mutate(avg = (old + young)/2)

range(misincorp$avg)
dfdat=lapply(10 ** seq(-7,-3,by=1),function(thres){
  misincorp %>% 
    filter(avg > thres) %>%
    summarise(
      n = length(old_increased),
      n_increased = sum(old_increased > 0), 
      n_decreased = sum(old_increased < 0),
      frac_increased = n_increased / n,
      frac_decreased = n_decreased / n,
      mean_increase = mean(old_increased),
      se=sd(old_increased)/sqrt(n),
      mw.p = wilcox.test(young, old, paired=T,alternative = "l")$p.value) %>%
    mutate(thres = thres)
}) %>% rbind.fill()

dfplot<-dfdat%>%select(thres,frac_increased,frac_decreased,mean_increase,se,mw.p)%>%
  melt(id.vars =c("thres","mean_increase","se","mw.p") )

plot<-ggplot(dfplot,aes(x=factor(thres),y=ifelse(variable=="frac_increased",value,-value)) )+
  geom_bar( aes(fill=variable),stat = "identity")+
  geom_point(aes(x=factor(thres), y=ifelse(variable=="frac_increased",mean_increase,NA) ) )+
  geom_errorbar(aes(ymin=ifelse(variable=="frac_increased",mean_increase-se,NA) ,
                    ymax=ifelse(variable=="frac_increased",mean_increase+se,NA)),width=0.3)+
  geom_line(data=dfplot%>%filter(variable=="frac_increased"),
            aes(x=factor(thres), y=mean_increase),group=1 )+
  scale_fill_manual(name="",values = c("#F39B7FFF", "#8491B4FF"),
                    breaks = c("frac_increased","frac_decreased"),
                    labels = c("Increased","Decreased"))+
  scale_y_continuous( breaks = c(-0.5,-0.25,0,0.25,0.5), labels = abs, 
                      sec.axis = sec_axis(name = "Aged-to-young relative increase\nin translation error rates",
                                          trans=~.*1,breaks =c(-0.5,-0.25,0,0.25,0.5),  
                                          labels = c(-0.5,-0.25,0,0.25,0.5)),
                      expand = expansion(mult = c(0.1,0.1)))+
  geom_hline(aes(yintercept = 0),linetype=5,color="#4DBBD5FF")+
  labs(x="Threshold for translation error rates (>x)",y="Fraction of translation error sites")+
  scale_x_discrete(labels = c(expression(10**-7),expression(10**-6),expression(10**-5),
                              expression(10**-4),expression(10**-3)),
                   breaks = c(10**-7,10**-6,10**-5,10**-4,10**-3))+
  theme_test()+
  labs(title = "Human")+theme(plot.title=element_text(hjust=0.5))+
  theme(axis.text.y.left  = element_text(color =c(rep("#8491B4FF",2),"#4DBBD5FF",rep("#F39B7FFF",2))),
        axis.text.y.right = element_text(color =c(rep("black",2),"#4DBBD5FF",rep("black",2))),
        legend.position = "none",legend.key = element_blank())



###

############################
#fig.s1
amino_acids = c("A","C","D","E","F","G","H","K","I/L",
                "M","N","P","Q","R","S","T","V","W","Y")
AllCodons<-c()
allnucle<-c("U","C","A","G")
for (a in allnucle){
  a=a
  for (b in allnucle) {
    b=b
    for (c in allnucle) {
      c=c
      df=paste(a,b,c,sep = "")
      AllCodons<-c(AllCodons,df)
    }
  }
}
amino_acids_64 = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG' %>%
  str_split("")%>%unlist()
dfcodon_AA_coord<-data.frame(codon=AllCodons,AA=amino_acids_64)%>%
  mutate(AA_poll=ifelse((AA=="I" | AA=="L") ,"I/L",AA))%>%
  mutate(AAcor=paste(codon,AA_poll,sep = " to "))

allSubsType<-c()
for (a in AllCodons){
  a=a
  for (b in amino_acids) {
    b=b
    df=paste(a,b,sep = " to ")
    allSubsType<-c(allSubsType,df)
    
  }
}

dfsubs_compile<-mclapply(1:nrow(dfqsubs),function(i){
  eachrow<-dfqsubs[i,]
  protein_index<-sample(1:length(str_split(eachrow$proteins," ")%>%unlist()),1 )
  codon<-((str_split(eachrow$codons," ")%>%unlist())[protein_index])%>%gsub("T","U",.)
  destination<-eachrow$destination
  
  dfBP<-eachrow%>%dplyr::select(c(starts_with("Base.intensity")))%>%t() %>%
    as.data.frame()%>%
    tibble::rownames_to_column()%>%
    tidyr::separate(col = "rowname",into = c("a","sample"),sep=".Exp..")%>%
    dplyr::select(2,3)
  colnames(dfBP)[2]<-"bp"
  
  dfDP<-eachrow%>%dplyr::select(c(starts_with("Modification.intensity")))%>%t()%>%
    as.data.frame()%>%
    tibble::rownames_to_column()%>%
    tidyr::separate(col = "rowname",into = c("a","sample"),sep=".Exp..")%>%
    dplyr::select(2,3)
  colnames(dfDP)[2]<-"dp"
  
  dfmerg<-merge(dfBP,dfDP)%>%dplyr::rowwise()%>%dplyr::mutate(dfratio=dp/sum(dp,bp,na.rm = T)) %>%
    dplyr::filter(dfratio<=0.01)
  
  young<-dfmerg%>%dplyr::filter(grepl("young",sample,ignore.case=T))%>%
    getElement("dfratio")%>%mean(na.rm=T)
  old<-dfmerg%>%dplyr::filter(grepl("old",sample,ignore.case=T))%>%
    getElement("dfratio")%>%mean(na.rm=T)
  
  return(data.frame(codon,destination,young,old))
},mc.cores = 50)%>%rbind.fill()%>%filter(young>0 | old >0)%>%na.omit()%>%
  mutate(destination=ifelse((destination=="I" | destination=="L") ,"I/L",destination))%>%
  mutate(substype=paste(codon,destination,sep = " to "))

dfmtr_onecodonToOneAA<-mclapply(allSubsType, function(thistype){
  eachtype<-dfsubs_compile%>%filter(substype==thistype)
  codon=strsplit(thistype," to ")[[1]][1]
  destination=strsplit(thistype," to ")[[1]][2]
  if(nrow(eachtype)==0){return(data.frame(codon,destination,effeSize=NA,pvalue=NA))}
  pvalue<- (wilcox.test(eachtype$young,eachtype$old,alternative="l"))$p.value 
  effeSize<-(mean(eachtype$old,na.rm=T)-mean(eachtype$young,na.rm=T))/
    (mean(eachtype$old,na.rm=T)+mean(eachtype$young,na.rm=T))
  codon = eachtype$codon%>%unique()
  destination = eachtype$destination%>%unique()
  return(data.frame(codon,destination,effeSize,pvalue))
},mc.cores = 60)%>%rbind.fill() 


dfmtr_allCodnToOneAA<-mclapply(amino_acids, function(thistype){
  eachtype<-dfsubs_compile%>%filter(destination==thistype)
  if(nrow(eachtype)==0){return(data.frame(codon="allCodon",destination=thistype,effeSize=NA,pvalue=NA))}
  pvalue<- (wilcox.test(eachtype$young,eachtype$old,alternative="l"))$p.value 
  effeSize<-(mean(eachtype$old,na.rm=T)-mean(eachtype$young,na.rm=T))/
    (mean(eachtype$old,na.rm=T)+mean(eachtype$young,na.rm=T))
  return(data.frame(codon="allCodon", destination=thistype,effeSize,pvalue))
},mc.cores = 60)%>%rbind.fill()


dfmtr_oneCodonToallAA<-mclapply(AllCodons, function(thistype){
  eachtype<-dfsubs_compile%>%filter(codon==thistype)
  if(nrow(eachtype)==0){return(data.frame(codon=thistype,destination="allAA",effeSize=NA,pvalue=NA))}
  pvalue<- (wilcox.test(eachtype$young,eachtype$old,alternative="l"))$p.value 
  effeSize<-(mean(eachtype$old,na.rm=T)-mean(eachtype$young,na.rm=T))/
    (mean(eachtype$old,na.rm=T)+mean(eachtype$young,na.rm=T))
  return(data.frame(codon=thistype,destination="allAA",effeSize,pvalue))
},mc.cores = 60)%>%rbind.fill()


dfallMtr<-rbind(dfmtr_onecodonToOneAA,dfmtr_allCodnToOneAA)%>%rbind(dfmtr_oneCodonToallAA)%>%
  mutate(effeSize=ifelse(is.na(effeSize),0,effeSize))%>%
  mutate(substype=paste(codon,destination,sep = " to "))%>%
  mutate(nobio=ifelse((substype %in% dfcodon_AA_coord$AAcor)|(codon %in%c("UAA","UAG","UGA")) ,T,F))%>%
  mutate(effeSize=ifelse((substype %in% dfcodon_AA_coord$AAcor)|(codon %in%c("UAA","UAG","UGA")) ,NA,effeSize))%>%
  mutate(stars=ifelse(pvalue<0.05,"*",""))%>%
  filter(!(codon %in%c("UAA","UAG","UGA") ))



ggplot(dfallMtr,aes(x=destination,y=codon,fill=effeSize))+
  geom_tile(aes(fill=effeSize))+
  geom_text(aes(label=stars),color="black",size=4,
            hjust="middle",vjust="bootom")+
  scale_fill_gradient2("Aged-to-young relative increase\nin translation error rates",
                       na.value="grey",low = "#4DBBD5FF",high = "#E64B35FF",mid = "white",
                       guide=guide_colorbar(title.position ="top" ,direction="horizontal",title.hjust=0.5,
                                            barwidth =8,barheight = 1, ticks.colour="black",ticks.linewidth = 1 ))+
  labs(x="Misincorporation of amino acids",y="Original codon")+
  scale_y_discrete(breaks=dfallMtr$codon%>%unique()%>%sort(),labels=c("",(dfallMtr$codon%>%unique()%>%sort())[-1]))+               
  theme(
    axis.ticks.y = element_blank(),
    panel.background=element_blank())+
  theme_test()+
  theme(legend.position = "top")


