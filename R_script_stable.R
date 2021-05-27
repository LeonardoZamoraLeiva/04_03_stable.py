
library(ggplot2)
library(gggenes)
library(data.table)

#my_file_to_draw <- commandArgs(trailingOnly = TRUE)
#test_file<-fread(my_file_to_draw)
out_test_smcogs_version1 <- fread('./out_test_smcogs_version1.csv')

antismash_colors<-c("#f16d75","#f16d75","#810e15","#808080","#2e8b57","#6495ed")

dummies <- make_alignment_dummies(out_test_smcogs_version1,aes(xmin = start, xmax = end, y = strain, id = gene_number),on = "gene1")
dummies['direction']<-1
dummies['antismash_smcog']<-'Additional biosynthetic'

tiff(file='test_file.tiff',units="in", width=15, height=5, res=250)
ggplot(out_test_smcogs_version1, aes(xmin = start, xmax = end, y = strain, fill = antismash_smcog, 
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