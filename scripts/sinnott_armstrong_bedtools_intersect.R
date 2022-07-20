library(reshape2)
library(stringr)
library(data.table)
library(dplyr)
library(tidyr)
library(qvalue)


source("scripts/intact_gse.R")


tissues <- read.table("data/hukku_et_al/tissues.tab",header=F,sep = '\t')$V1


for (i in 198:201){

        system(paste0('bedtools intersect -a ',
                snakemake@input[[i]],
                ' -b ',snakemake@input[[197]],' -wa -wb > ',
                snakemake@output[[i-197]])
                )
}

print("bedtools intersection done")


molecular_gene_list <- lapply(c("Urate_all","IGF_1_all","Testosterone_male","Testosterone_female"),FUN = sig_genes_list,type = "molecular")

molecular_gene_hits <- data.frame()

for (i in 1:length(molecular_gene_list)){
  
  molecular_gene_hits <- rbind.data.frame(molecular_gene_hits, molecular_gene_list[[i]])
  
}


annotLookup <- read.table(snakemake@input[[202]],sep = '\t',header = T)


molecular_gene_hits %>% dplyr::rename(ensembl_gene_id = gene) %>% 
  left_join(annotLookup,by = "ensembl_gene_id")  %>%
  dplyr::select(trait,hgnc_symbol) %>%
  filter(complete.cases(.)) %>% 
  filter(hgnc_symbol != "") %>% 
  write.table(snakemake@output[[5]],sep='\t',col.names = T,row.names = F,quote = F)
