library(reshape2)
library(stringr)
library(data.table)
library(dplyr)
library(tidyr)
library(qvalue)


source("scripts/intact_gse.R")

tissues <- read.table("data/hukku_et_al/tissues.tab",header=F,sep = '\t')$V1


trait_tissue <- read.table("data/sinnott_armstrong_et_al/twas_coloc_sumstats/trait_tissue_molecular.txt",header=F)$V1


annotLookup <- fread(snakemake@input[[210]],sep = '\t',header = T)


#Table 3



traits <- c("Urate_all","IGF_1_all","Testosterone_male","Testosterone_female")
outdf <- data.frame()



for (i in 1:length(traits)){

	hit_genes <- read.table(snakemake@input[[i+200]],sep='\t')


	intact_genes <- read.table(snakemake@input[[205]],sep='\t',header=T) %>%
        	filter(trait == traits[i])


	#number of intersecting GWAS hits with INTACT gene regions

	n_intact <- hit_genes %>% filter(V7 %in% intact_genes$hgnc_symbol) %>% select(V1,V2,V3) %>% distinct() %>% nrow()

	out <- data.frame(traits[i],n_intact)

	outdf <- rbind.data.frame(outdf,out)

}


colnames(outdf) <- c('Trait','Number_GWAS_hits_100kb_to_INTACT_gene')


#Find out how many independent GWAS hits are within 100kbp of a core gene

urate_hit_genes <- read.table(snakemake@input[[201]],sep='\t')
igf1_hit_genes <- read.table(snakemake@input[[202]],sep='\t')
test_m_hit_genes <- read.table(snakemake@input[[203]],sep='\t')
test_f_hit_genes <- read.table(snakemake@input[[204]],sep='\t')



urate_core_genes <- read.table(snakemake@input[[198]],sep='\t',header=T)
igf1_core_genes <- read.table(snakemake@input[[199]],sep='\t',header=T)
testosterone_core_genes <- read.table(snakemake@input[[200]],sep='\t',header=T)



urate_ind_hits <- read.table(snakemake@input[[206]],header=T)
igf1_ind_hits <- read.table(snakemake@input[[207]],header=T)
test_m_ind_hits <- read.table(snakemake@input[[208]],header=T)
test_f_ind_hits <- read.table(snakemake@input[[209]],header=T)


#number of intersecting GWAS hits with core gene regions

n_urate_core <- urate_hit_genes %>% filter(V7 %in% urate_core_genes$Gene.symbol) %>% select(V1,V2,V3) %>% distinct() %>% nrow()
n_igf1_core <- igf1_hit_genes %>% filter(V7 %in% igf1_core_genes$Gene.name) %>% select(V1,V2,V3) %>% distinct() %>% nrow()
n_test_m_core <- test_m_hit_genes %>% filter(V7 %in% testosterone_core_genes$Gene.name) %>% select(V1,V2,V3) %>% distinct() %>% nrow()
n_test_f_core <- test_f_hit_genes %>% filter(V7 %in% testosterone_core_genes$Gene.name) %>% select(V1,V2,V3) %>% distinct() %>% nrow()


out_table <- data.frame(Trait = c("Urate_all","IGF_1_all","Testosterone_male","Testosterone_female"),
                        N_hits = c(nrow(urate_ind_hits),nrow(igf1_ind_hits),nrow(test_m_ind_hits),nrow(test_f_ind_hits)),
                        N_core_hits = c(n_urate_core,n_igf1_core,n_test_m_core,n_test_f_core)) %>% mutate(prop = N_core_hits/N_hits)



df_final <- merge(outdf,out_table,by = "Trait") %>% select(Trait,N_hits,Number_GWAS_hits_100kb_to_INTACT_gene,N_core_hits)


write.table(df_final[c(4,1,3,2),],snakemake@output[[1]],col.names = T, row.names = F, quote = F, sep = '\t')
	    
	   
	    
	    
	    
	    
	    
	    
	    
	    
	    
	    
#Table 4 and Supplemental Excel File 2



urate_file <- urate_core_genes 
igf_file <- igf1_core_genes
testosterone_file <- testosterone_core_genes

colnames(urate_file)[2] <- "hgnc_symbol"
colnames(igf_file)[2] <- "hgnc_symbol"
colnames(testosterone_file)[4] <- "hgnc_symbol"

urate_ensID <- left_join(urate_file,annotLookup,by = "hgnc_symbol")
colnames(urate_ensID)[3] <- 'Gene'

igf1_ensID <- left_join(igf_file,annotLookup,by='hgnc_symbol')
colnames(igf1_ensID)[4] <- "Gene"
colnames(igf1_ensID)[1] <- "Pathway"

testosterone_ensID <- left_join(testosterone_file,annotLookup,by='hgnc_symbol')
colnames(testosterone_ensID)[5] <- "Gene"
colnames(testosterone_ensID)[1] <- "Pathway"


urate_pathways <- unique(as.character(urate_ensID$Pathway))
igf1_pathways <- unique(as.character(igf1_ensID$Pathway))
testosterone_pathways <- unique(as.character(testosterone_ensID$Pathway))




core_genes_rst_table <- data.frame()


test_urate <- ev_int_fdr("Urate_all",prior = "linear",type = "molecular")
colnames(test_urate)[1] <- "ensembl_gene_id"
test <- left_join(test_urate,annotLookup,by="ensembl_gene_id")
print(length(unique(test$hgnc_symbol[test$hgnc_symbol %in% urate_ensID$hgnc_symbol & test$sig == 1])))
for (i in 1:length(urate_pathways)){
  n_genes <- length(unique(test$hgnc_symbol[test$hgnc_symbol %in% urate_ensID$hgnc_symbol[urate_ensID$Pathway == urate_pathways[i]] & test$sig == 1]))
  outdf <- data.frame(trait = rep("Urate_all",n_genes),
                      path = rep(urate_pathways[i],n_genes),
                      gene_hits = unique(test$hgnc_symbol[test$hgnc_symbol %in% urate_ensID$hgnc_symbol[urate_ensID$Pathway == urate_pathways[i]] & test$sig == 1]))
  core_genes_rst_table <- rbind.data.frame(core_genes_rst_table,outdf)
}
#unique(test$hgnc_symbol[test$hgnc_symbol %in% urate_ensID$hgnc_symbol & test$sig == 1])



test_igf1 <- ev_int_fdr("IGF_1_all",prior = "linear",type = "molecular")
colnames(test_igf1)[1] <- "ensembl_gene_id"
test <- left_join(test_igf1,annotLookup,by="ensembl_gene_id")
print(length(unique(test$hgnc_symbol[test$hgnc_symbol %in% igf1_ensID$hgnc_symbol & test$sig == 1])))
for (i in 1:length(igf1_pathways)){
  gene_hits <- length(unique(test$hgnc_symbol[test$hgnc_symbol %in% igf1_ensID$hgnc_symbol[igf1_ensID$Pathway == igf1_pathways[i]] & test$sig == 1]))
  outdf <- data.frame(trait = rep("IGF_1_all",gene_hits),
                      path = rep(igf1_pathways[i],gene_hits),
                      gene_hits = unique(test$hgnc_symbol[test$hgnc_symbol %in% igf1_ensID$hgnc_symbol[igf1_ensID$Pathway == igf1_pathways[i]] & test$sig == 1]))
  core_genes_rst_table <- rbind.data.frame(core_genes_rst_table,outdf)
}
#unique(test$hgnc_symbol[test$hgnc_symbol %in% igf1_ensID$hgnc_symbol & test$sig == 1])


test_test_m <- ev_int_fdr("Testosterone_male",prior = "linear",type = "molecular")
colnames(test_test_m)[1] <- "ensembl_gene_id"
test <- left_join(test_test_m,annotLookup,by="ensembl_gene_id")
print(length(unique(test$hgnc_symbol[test$hgnc_symbol %in% testosterone_ensID$hgnc_symbol & test$sig == 1])))
for (i in 1:length(testosterone_pathways)){
  gene_hits <- length(unique(test$hgnc_symbol[test$hgnc_symbol %in% testosterone_ensID$hgnc_symbol[testosterone_ensID$Pathway ==     testosterone_pathways[i]] & test$sig == 1]))
  outdf <- data.frame(trait = rep("Testosterone_male",gene_hits),
                      path = rep(testosterone_pathways[i],gene_hits),
                      gene_hits = unique(test$hgnc_symbol[test$hgnc_symbol %in% testosterone_ensID$hgnc_symbol[testosterone_ensID$Pathway == testosterone_pathways[i]] & test$sig == 1]))
  core_genes_rst_table <- rbind.data.frame(core_genes_rst_table,outdf)
}
unique(test$hgnc_symbol[test$hgnc_symbol %in% testosterone_ensID$hgnc_symbol & test$sig == 1])



test_test_f <- ev_int_fdr("Testosterone_female",prior = "linear",type = "molecular")
colnames(test_test_f)[1] <- "ensembl_gene_id"
test <- left_join(test_test_f,annotLookup,by="ensembl_gene_id")
print(length(unique(test$hgnc_symbol[test$hgnc_symbol %in% testosterone_ensID$hgnc_symbol & test$sig == 1])))
for (i in 1:length(testosterone_pathways)){
  gene_hits <- length(unique(test$hgnc_symbol[test$hgnc_symbol %in% testosterone_ensID$hgnc_symbol[testosterone_ensID$Pathway ==     testosterone_pathways[i]] & test$sig == 1]))
  outdf <- data.frame(trait = rep("Testosterone_female",gene_hits),
                      path = rep(testosterone_pathways[i],gene_hits),
                      gene_hits = unique(test$hgnc_symbol[test$hgnc_symbol %in% testosterone_ensID$hgnc_symbol[testosterone_ensID$Pathway ==     testosterone_pathways[i]] & test$sig == 1]))
  core_genes_rst_table <- rbind.data.frame(core_genes_rst_table,outdf)
}

print(dim(core_genes_rst_table))

#colnames(core_genes_rst_table) <- c("Trait","Pathway","Ge")

#Find intersection between proximity hits and our hits
purine_genes <- c("ADCY5","PDE4C",'PDE1A',"NME2","NME6","PKLR",
                  "AK1","AK7","AMPD3","IMPDH2","NT5C1B","NT5C2",
                  "XDH", "GUCY2D","ADCY5","PKLR","PDE1A","PDE9A", "HPRT1")
purine_overlap <- intersect(unique(purine_genes),core_genes_rst_table$gene_hits[core_genes_rst_table$path=="Purine metabolism"])

solute_genes <- c("SLC2A9","ABCG2","SLC22A11","SLC22A7","PDZK1","SLC17A1"
                  ,"ABCC4","SLC22A12","SLC16A9")

solute_overlap <- intersect(unique(solute_genes),core_genes_rst_table$gene_hits[core_genes_rst_table$path=="Solute transport"])

GH_genes <- c("GHSR","CREB3L2","EP300","SSTR2","SSTR5","GNAI3","ADCY1","ADCY3","ADCY4","GNAS","GHRHR","GH1","GH2","FOS")

GH_overlap <- intersect(unique(GH_genes),core_genes_rst_table$gene_hits[core_genes_rst_table$path=="Growth hormone secretion"])

secretion_genes <- c("GHR","JAK1","JAK2","TYK2","SH2B1","STAT5B","STAT5A","STAT3","SOCS2","EP300")

secretion_overlap <- intersect(unique(secretion_genes),core_genes_rst_table$gene_hits[core_genes_rst_table$path=="IGF-1 secretion"])

SB_genes <- c("INS","INSL4","INSL6","IGF1","IGF2","INS-IGF2","IGFBP1","IGFBP3","IGFALS","PAPPA2")

SB_overlap <- intersect(unique(SB_genes),core_genes_rst_table$gene_hits[core_genes_rst_table$path=="IGF-1 serum balance"])

DS_genes <- c("IRS1","IRS2","SH2B1","PIK3CA","AKT3","FOXO3","IGF1R","IGF2R","INSR")

DS_overlap <- intersect(unique(DS_genes),core_genes_rst_table$gene_hits[core_genes_rst_table$path=="Downstream signaling"])

ras_genes <- c("MRAS","RASA2","RASIP1","RAF1","MAP3K4","MAPK3","MAP4K1","MAPKAP1","MAPK8IP3","RIN2","SHC1","RAB26","RAB3D","RAB7L1","RAB3IL1","RAB5B","RABEP2")


ras_overlap <- intersect(unique(ras_genes),core_genes_rst_table$gene_hits[core_genes_rst_table$path=="Ras signaling"])

SBS_m <- c("UGT2B7","UGT2B11",paste0("UGT1A",seq(1,10)),"UGT2A3",
           "UGT2B10","UGT2B15","UGT2B28","SRD5A2","AKR1C3","AKR1C4",
           "CYP19A1","HSD17B1")

SBS_m_overlap <- intersect(SBS_m, core_genes_rst_table$gene_hits[core_genes_rst_table$path == "Synthesis and metabolism" & core_genes_rst_table$trait == "Testosterone_male"])


SBS_f <- c("SOAT1","CYP17A1","CYP21A2","AKR1C1","AKR1C2","AKR1C3","CYP1A2","CYP11B1","CYP11B2","CYP3A4","CYP3A5","UGT2B7","UGT2B11")

SBS_f_overlap <- intersect(SBS_f, core_genes_rst_table$gene_hits[core_genes_rst_table$path == "Synthesis and metabolism" & core_genes_rst_table$trait == "Testosterone_female"])

#Generate results table with number of core genes for each trait and pathway







#Supplemental Excel File 2

proximity_rst_table <- data.frame(trait = c(rep("Urate_all",length(c(purine_genes,solute_genes))),
                                            rep("IGF_1_all",length(c(GH_genes,secretion_genes,SB_genes,DS_genes,ras_genes))),
                                            rep("Testosterone_male",length(SBS_m)),
                                            rep("Testosterone_female",length(SBS_f))),
                                  path = c(rep("Purine metabolism",length(purine_genes)),
                                            rep("Solute transport",length(solute_genes)),
                                            rep("Growth hormone secretion",length(GH_genes)),
                                            rep("IGF-1 secretion",length(secretion_genes)),
                                            rep("IGF-1 serum balance",length(SB_genes)),
                                            rep("Downstream signaling",length(DS_genes)),
                                            rep("Ras signaling",length(ras_genes)),
                                            rep("Synthesis and metabolism",length(SBS_m)),
                                            rep("Synthesis and metabolism",length(SBS_f))),
                                  gene_hits = c(purine_genes,
                                                solute_genes,
                                                GH_genes,
                                                secretion_genes,
                                                SB_genes,
                                                DS_genes,
                                                ras_genes,
                                                SBS_m,
                                                SBS_f),
                                  proximity_flag = 1)

core_genes_rst_table %>%
  mutate(intact_flag = 1) %>%
  merge(proximity_rst_table, by = c("trait","path","gene_hits"),all = T) %>%
  mutate(implicated_by_intact = case_when(is.na(intact_flag) ~ "NO",
                                 intact_flag == 1 ~ "YES")) %>%
  mutate(implicated_by_proximity = case_when(is.na(proximity_flag)  ~ "NO",
                                    proximity_flag == 1 ~ "YES")) %>%
  dplyr::select(trait,path,gene_hits,implicated_by_intact,implicated_by_proximity) %>%
  dplyr::rename(Trait = trait,
                Core_gene = gene_hits,
                Pathway = path,
                Implicated_by_intact = implicated_by_intact,
                Implicated_by_proximity = implicated_by_proximity) %>%
  unique() %>%
  filter(!(Pathway %in% c("No pathway", "Serum homeostasis","HPG signaling"))) %>%
  write.table(snakemake@output[[2]],row.names = F,col.names = T,quote = F,sep ='\t')



core_genes_rst_table %>%
  mutate(intact_flag = 1) %>% 
  merge(proximity_rst_table, by = c("trait","path","gene_hits"),all = T) %>%
  mutate(implicated_by_intact = case_when(is.na(intact_flag) ~ "NO",
                                 intact_flag == 1 ~ "YES")) %>%
  mutate(implicated_by_proximity = case_when(is.na(proximity_flag)  ~ "NO",
                                    proximity_flag == 1 ~ "YES")) %>%
  dplyr::select(trait,path,gene_hits,implicated_by_intact,implicated_by_proximity) %>%
  dplyr::rename(Trait = trait,
                Core_gene = gene_hits,
                Pathway = path,
                Implicated_by_intact = implicated_by_intact,
                Implicated_by_proximity = implicated_by_proximity) %>%
  unique() %>%
  filter(!(Pathway %in% c("No pathway", "Serum homeostasis","HPG signaling"))) %>%
  group_by(Trait,Pathway) %>%
  summarise(n_proximity = length(which(Implicated_by_proximity == "YES")),
            n_intact = length(which(Implicated_by_intact == "YES")), 
            n_overlap = length(which(Implicated_by_intact == "YES" & Implicated_by_proximity == "YES"))) %>%
  write.table(snakemake@output[[3]],row.names = F, col.names = T, quote = F, sep = '\t')









#Table 5

enrich_fun <- function(pathway_name, trait_df, df){
  
  print(pathway_name)
  
  pathway_genes <- as.character(trait_df$Gene[trait_df$Pathway==pathway_name])
  
  print(pathway_genes)
  #indicator for membership
  
  df$pathway_anno <- 0
  
  df$pathway_anno[df$Gene %in% pathway_genes] <- 1
  
  print(df[df$pathway_anno == 1],)
  print(df[df$Posterior > 0],)
  

  #perform gene set enrichment using posteriors and annotation vector
  
  out <- enrich_res2(sig_lev = 0.05,
            pprobs = df$Posterior, 
            d_vec = df$pathway_anno,
            SE_type = 'NDS')
  
  out$pathway <- pathway_name
  
  out$gene_set_size <- length(pathway_genes)
  
  out$prop_rep <- length(which(pathway_genes %in% df$Gene))/length(pathway_genes)
  
  return(out)
}





enrich_fun_testosterone_f <- function(tissue){
 
  dat_idx <- which(tissues == tissue) + 49
  
  ev_int_df <- fread(snakemake@input[[dat_idx]])

  ev_int_df$Posterior = evid_int(GLCP_vec = ev_int_df$GLCP,
                                 prior_fun = truncated_linear,
                                 z_vec = ev_int_df$TWAS_z,
                                 u = pi1_fun(z_vec = ev_int_df$TWAS_z,lambda = 0.5),
                                 t = 0.05)


  out <- data.frame()

  for (i in 1:length(testosterone_pathways)){
      rst <- enrich_fun(pathway_name = testosterone_pathways[i], trait_df = testosterone_ensID, df = ev_int_df)
      rst$Trait <- "Testosterone_female"
      rst$Tissue <- tissue
      out <- rbind.data.frame(out,rst)
  }

  out <- out %>% arrange(pval,Tissue)

  return(out)

}

#testostone male

enrich_fun_testosterone_m <- function(tissue){

  dat_idx <- which(tissues == tissue) + 49*2

  ev_int_df <- fread(snakemake@input[[dat_idx]])

  ev_int_df$Posterior = evid_int(GLCP_vec = ev_int_df$GLCP,
                                 prior_fun = truncated_linear,
                                 z_vec = ev_int_df$TWAS_z,
                                 u = pi1_fun(z_vec = ev_int_df$TWAS_z,lambda = 0.5),
                                 t = 0.05)

  out <- data.frame()

  for (i in 1:length(testosterone_pathways)){
      rst <- enrich_fun(pathway_name = testosterone_pathways[i], trait_df = testosterone_ensID, df = ev_int_df)
      rst$Trait <- "Testosterone_male"
      rst$Tissue <- tissue
      out <- rbind.data.frame(out,rst)
  }

  out <- out %>% arrange(pval,Tissue)

  return(out)

  }


#urate

enrich_fun_urate <- function(tissue){

  dat_idx <- which(tissues == tissue) + 49*3

  ev_int_df <- fread(snakemake@input[[dat_idx]])

  ev_int_df$Posterior = evid_int(GLCP_vec = ev_int_df$GLCP,
                                 prior_fun = truncated_linear,
                                 z_vec = ev_int_df$TWAS_z,
                                 u = pi1_fun(z_vec = ev_int_df$TWAS_z,lambda = 0.5),
                                 t = 0.05)
  out <- data.frame()

  for (i in 1:length(urate_pathways)){
      rst <- enrich_fun(pathway_name = urate_pathways[i], trait_df = urate_ensID, df = ev_int_df)
      rst$Trait <- "Urate_all"
      rst$Tissue <- tissue
      out <- rbind.data.frame(out,rst)
  }

  out <- out %>% arrange(pval,Tissue)

  return(out)
}

#IGF_1

enrich_fun_igf1 <- function(tissue){
 
  dat_idx <- which(tissues == tissue) 

  ev_int_df <- fread(snakemake@input[[dat_idx]])

  ev_int_df$Posterior = evid_int(GLCP_vec = ev_int_df$GLCP,
                                 prior_fun = truncated_linear,
                                 z_vec = ev_int_df$TWAS_z,
                                 u = pi1_fun(z_vec = ev_int_df$TWAS_z,lambda = 0.5),
                                 t = 0.05)
  out <- data.frame()

  for (i in 1:length(igf1_pathways)){
      rst <- enrich_fun(pathway_name = igf1_pathways[i],trait_df = igf1_ensID, df = ev_int_df)
      rst$Trait <- "IGF_1_all"
      rst$Tissue <- tissue
      out <- rbind.data.frame(out,rst)
  }

  out <- out %>% arrange(pval,Tissue)

  return(out)
}




#compile results for all tissues into one file; sort by p value and estimate

#testosterone f

testosterone_f_gse_rst <- lapply(tissues, FUN = enrich_fun_testosterone_f)

testosterone_f_rst <- data.frame()

for (i in 1:length(testosterone_f_gse_rst)){

  testosterone_f_rst <- rbind.data.frame(testosterone_f_rst,testosterone_f_gse_rst[[i]])

}


#testosterone m

testosterone_m_gse_rst <- lapply(tissues, FUN = enrich_fun_testosterone_m)

testosterone_m_rst <- data.frame()

for (i in 1:length(testosterone_m_gse_rst)){

  testosterone_m_rst <- rbind.data.frame(testosterone_m_rst,testosterone_m_gse_rst[[i]])

}




#urate

urate_gse_rst <- lapply(tissues, FUN = enrich_fun_urate)

urate_rst <- data.frame()

for (i in 1:length(urate_gse_rst)){

  urate_rst <- rbind.data.frame(urate_rst,urate_gse_rst[[i]])

}




#IGF1

igf1_gse_rst <- lapply(tissues, FUN = enrich_fun_igf1)

igf1_rst <- data.frame()

for (i in 1:length(igf1_gse_rst)){

  igf1_rst <- rbind.data.frame(igf1_rst,igf1_gse_rst[[i]])

}


molecular_trait_gse_rst <- rbind.data.frame(urate_rst,igf1_rst,testosterone_m_rst,testosterone_f_rst) %>%
  dplyr::select(Trait,Tissue,pathway,Estimate,SE,pval,CONVERGED) %>%
  filter(pathway != "No pathway") %>%
  arrange(Trait,Tissue,desc(Estimate))

write.table(molecular_trait_gse_rst,file = snakemake@output[[4]],sep = "\t",row.names = F,col.names = T,quote = F)




#Table 5 (filter down to tissues of interest for each trait)



molecular_trait_gse_rst %>% 
	filter((Trait == "Urate_all" & 
		pathway %in% c("Solute transport","Purine metabolism") & 
		Tissue == "Kidney_Cortex") |
	       (Trait == "IGF_1_all" & 
		pathway %in% c("IGF-1 serum balance",
			       "Downstream signaling",
			       "IGF-1 secretion",
			       "Ras signaling",
			       "Growth hormone secretion") &
		Tissue %in% c("Pituitary","Liver")) |
	       (Trait == "Testosterone_male" & 
		pathway %in% c("Synthesis and metabolism",
			       "Serum homeostasis",
			       "HPG signaling") &
		Tissue %in% c("Brain_Hypothalamus",
			      "Pituitary",
			      "Testis",
			      "Liver")) |
	       (Trait == "Testosterone_female" &
		pathway %in% c("Brain_Hypothalamus",
			       "Serum homeostasis",
			       "HPG signaling") & 
		Tissue %in% c("Brain_Hypothalamus",
			      "Pituitary",
			      "Adrenal_Gland",
			      "Ovary",
			      "Liver"))) %>% 
	write.table(snakemake@output[[5]],sep = '\t',row.names = F, col.names = T, quote = F)





#Table S6


molecular_table_list <- lapply(c("Urate_all","IGF_1_all","Testosterone_male","Testosterone_female"),FUN = sig_genes_fun,type = "molecular")

molecular_table <- data.frame()

for (i in 1:length(molecular_table_list)){
  
  molecular_table <- rbind.data.frame(molecular_table, molecular_table_list[[i]])
  
}

molecular_table %>% mutate(post_hoc_col = paste0(post_hoc_n_hits, " (",post_hoc_unique_hits,")"))%>%
  mutate(step_col = paste0(step_n_hits, " (",step_n_unique_hits,")")) %>%
  mutate(linear_col = paste0(linear_n_hits, " (", linear_unique_hits,")")) %>%
  rbind.data.frame(data.frame(trait = "Total",
                              post_hoc_n_hits = NA,
                              post_hoc_unique_hits = NA,
                              step_n_hits = NA,
                              step_n_unique_hits = NA,
                              linear_n_hits = NA,
                              linear_unique_hits = NA,
                              post_hoc_col = sum(molecular_table$post_hoc_n_hits),
                              step_col = sum(molecular_table$step_n_hits),
                              linear_col = sum(molecular_table$linear_n_hits))) %>%
  dplyr::select(trait,post_hoc_col,step_col,linear_col) %>%
  dplyr::rename("Trait" = trait,
         "Hukku et al. 2022" = post_hoc_col,
         "Step Prior" = step_col,
         "Linear Prior" = linear_col) %>%
  mutate(Trait = c("Serum urate","IGF-1","Testosterone (male)","Testosterone (female)","Total")) %>%
  write.table(snakemake@output[[6]],sep = '\t',row.names = F, col.names = T, quote = F)

