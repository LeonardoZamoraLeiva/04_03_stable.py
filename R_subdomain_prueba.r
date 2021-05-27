#!/usr/bin/env Rscript

library(ggplot2)
library(gggenes)
library(data.table)


setwd("C:/Users/Leo/PycharmProjects/Gene_Draw")

#test_file_2 <- subdomains_09_03_modified

antismash_colors<-c("#810e15","gold1","skyblue1","gray32","palegreen1","seagreen4","sienna","coral1","deeppink1","dodgerblue3","gray72","gray48")
                    
                    #"gold1","skyblue1","royalblue4","palegreen1","seagreen4","sienna","coral1","deeppink4","dodgerblue3","white","white")
                    
                    
                    #"darkorange3","lightpink1","royalblue4","palegreen1","seagreen4","cadetblue3","sienna4","gold1","hotpink3","white","white")
                    

dummies <- make_alignment_dummies(subdomains_09_03,aes(xmin = start, xmax = end, y = strain, id = gene_number),on = "gene50")
dummies['direction']<-1
dummies['antismash_smcog']<-'biosynthetic-additional'

tiff(file='./figures_out/09_04_test_04.tiff', units="in",width=15, height=5, res=250)

ggplot(subdomains_09_03, aes(xmin = start, xmax = end, y = strain,forward = direction)) +
  geom_gene_arrow(fill  = "white") +
  geom_subgene_arrow(data= subdomains_09_03,
                     aes(xmin=start,xmax=end,y=strain,fill=Aminoacid,
                         xsubmin=from,xsubmax=to),color='black',alpha=1)+ 
                     
  #geom_subgene_label(data= subdomains_09_03,
  #                   align = "left",
  #                   aes(xmin=start,xmax=end,y=strain,
  #                       xsubmin=from,xsubmax=to,label=Aminoacid),min.size = 0) +
  
  facet_wrap(~ strain, scales = "free", ncol = 1) +
  theme(legend.title = element_blank())+
  geom_blank(data = dummies) +
  scale_fill_manual(values=antismash_colors)+
  ylab('Strain')+
  xlab('Position (bp)')+
  theme_genes()
dev.off()