#!/usr/bin/env Rscript

library(ggplot2)
library(gggenes)
library(data.table)


setwd("C:/Users/Leo/PycharmProjects/Gene_Draw")

test_file <- read.csv('./out/mibig_close_to_g11c.csv')

antismash_colors<-c("#810e15","#f16d75","#808080","#2e8b57","#6495ed")

dummies <- make_alignment_dummies(test_file,aes(xmin = start, xmax = end, y = strain, id = gene_number),on = "gene1")
dummies['direction']<-1
dummies['antismash_smcog']<-'biosynthetic-additional'

tiff(file='./figures_out/mibig_close_to_g11c.tiff', units="in",width=15, height=5, res=250)

ggplot(test_file, aes(xmin = start, xmax = end, y = strain, fill = antismash_smcog, 
                          forward = direction)) +
  geom_gene_arrow() +
  facet_wrap(~ strain, scales = "free", ncol = 1) +
  theme(legend.title = element_blank())+
  geom_blank(data = dummies) +
  scale_fill_manual(values=antismash_colors)+
  ylab('Strain')+
  xlab('Position (bp)')+
  theme_genes()
dev.off()
