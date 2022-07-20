
library(reshape2)
library(stringr)
library(data.table)
library(dplyr)
library(tidyr)
library(qvalue)

source("scripts/intact_gse.R")

tissues <- read.table("data/hukku_et_al/tissues.tab",header=F,sep = '\t')$V1

#Generate Table 2

trait_tissue <- read.table("data/hukku_et_al/trait_tissue_nonmolecular.txt",header=F)$V1



#Top CAD-relevant pathways based on CTD


path_ids <- c("REACT:R-HSA-1430728",
	      "REACT:R-HSA-162582",
	      "REACT:R-HSA-168256",
 	      "REACT:R-HSA-556833",
	      "REACT:R-HSA-1280215",
	      "REACT:R-HSA-392499",
	      "REACT:R-HSA-73923",
	      "REACT:R-HSA-449147",
	      "REACT:R-HSA-109582",
	      "REACT:R-HSA-168249",
	      "REACT:R-HSA-1280218",
	      "REACT:R-HSA-174824")



cad_tissues <- c("Liver", "Adipose_Subcutaneous", "Artery_Coronary", "Adipose_Visceral_Omentum","Muscle_Skeletal")


ctd_genes <- read.table(snakemake@input[[197]],sep=',')

colnames(ctd_genes) <- c("hgnc_symbol","gene_id","pathway_name","pathway_id")


annotLookup <- read.table(snakemake@input[[199]],sep = '\t',header = T)



cad_res <- data.frame()

for (i in 1:length(cad_tissues)){
	  for (j in 1:length(path_ids)){
		      
		      cad_res <- rbind.data.frame(cad_res, cad_gse(path_ids[j], cad_tissues[i]))
  }
}

cad_res2 <- cad_res %>% 
	dplyr::select(Estimate,CI_Leftlim,CI_Rightlim,path,path_def,gene_set_size,CONVERGED,prop_rep,tissue) %>% arrange(tissue,path) %>% mutate(path_def = paste0(path_def," (",gene_set_size,")"))


cad_table <- cad_res2 %>% 
	dplyr::select(Estimate,path_def,CONVERGED,path,tissue,CI_Leftlim,CI_Rightlim) %>%
	  mutate(Estimate = round(Estimate,digits = 3),
		          CI_Leftlim = round(CI_Leftlim,digits = 3),
			           CI_Rightlim = round(CI_Rightlim,digits = 3)) %>% 
  mutate(path = str_sub(path,start = 7)) %>%
    mutate(Estimate_final = case_when(CI_Leftlim < 0 & CONVERGED == 1 ~ paste0(Estimate," (",CI_Leftlim,",",CI_Rightlim,")"),
	 CI_Leftlim > 0 & CONVERGED == 1 ~ paste0(Estimate,"* ","(",CI_Leftlim,",",CI_Rightlim,")"),  CONVERGED == 0 ~ "NI")) %>% 
    dplyr::select(-c(CI_Leftlim,CI_Rightlim,Estimate,CONVERGED)) %>% 
      spread(key = tissue, value = Estimate_final) 

colnames(cad_table) <- c("Pathway","Reactome ID","Subcutaneous Adipose","Visceral Adipose","Coronary Artery","Liver","Skeletal Muscle")


for (i in 3:ncol(cad_table)){
	  
	  cad_table[[i]] <- ifelse(cad_table[[i]] == "0 (-1.96,1.96)", "NI", cad_table[[i]])
  
}

cad_table %>% filter(!(Pathway %in% c("Hemostasis (359)","Signaling by Interleukins (287)"))) %>% 
	  write.table(snakemake@output[[1]],quote = F, row.names = F, col.names = T, sep = '\t')


 









#Table S4




non_molecular_table_list <- lapply(c("Heights","HDL","LDL","CAD"),FUN = sig_genes_fun,type = "nonmolecular")


non_molecular_table <- data.frame()

for (i in 1:length(non_molecular_table_list)){
	  
	  non_molecular_table <- rbind.data.frame(non_molecular_table, non_molecular_table_list[[i]])
  
}

non_molecular_table %>% mutate(post_hoc_col = paste0(post_hoc_n_hits, " (",post_hoc_unique_hits,")"))%>%
	  mutate(step_col = paste0(step_n_hits, " (",step_n_unique_hits,")")) %>%
	    mutate(linear_col = paste0(linear_n_hits, " (", linear_unique_hits,")")) %>%
	      rbind.data.frame(data.frame(trait = "Total",
					  post_hoc_n_hits = NA,
					  post_hoc_unique_hits = NA,
					  step_n_hits = NA,
					  step_n_unique_hits = NA,
					  linear_n_hits = NA,
				          linear_unique_hits = NA,
					  post_hoc_col = sum(non_molecular_table$post_hoc_n_hits),
					  step_col = sum(non_molecular_table$step_n_hits),
				          linear_col = sum(non_molecular_table$linear_n_hits))) %>%
  dplyr::select(trait,post_hoc_col,step_col,linear_col) %>%
    dplyr::rename("Trait" = trait,
		  "Hukku et al. 2022" = post_hoc_col,
	          "Step Prior" = step_col,
	          "Linear Prior" = linear_col) %>%
    write.table(snakemake@output[[2]],quote = F, row.names = F, col.names = T, sep = '\t')
  



