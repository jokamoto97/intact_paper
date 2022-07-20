library(reshape2)
library(stringr)
library(data.table)
library(dplyr)
library(tidyr)
library(qvalue)


ttc <- read.table("data/hukku_et_al/trait_tissue_nonmolecular.txt",header=F)$V1


rst_df <- data.frame()

outdf <- data.frame()

for (i in 1:length(snakemake@input[1:196])){

                print(snakemake@input[[i]])

                real_dat_res <- fread(snakemake@input[[i]]) %>% 
			select(Parameter, Estimate, SE, z, pval, CI_Leftlim, CI_Rightlim, CONVERGED, GO_term, term_size, prop_rep)
		 
		real_dat_res$tissue_trait <- ttc[i]

                rst_df <- rbind.data.frame(rst_df,real_dat_res)

                outdf <- rbind.data.frame(outdf,
                                          data.frame(Trait_tissue = ttc[i],
                                                     n_GO_terms_tested = nrow(real_dat_res),
                                                     n_tests_converged = length(which(real_dat_res$CONVERGED == 1)),
                                                     n_sig05 = length(which(real_dat_res$pval < 0.05))))


}

#Supplemental Excel File 1

rst_df %>% select(Estimate,pval,CONVERGED,GO_term,term_size,SE,tissue_trait) %>%
           extract(tissue_trait, into = c("Trait","Tissue"), "([^_]+)_(.*)$") %>%
           select(Trait,Tissue,GO_term,Estimate,SE,pval,CONVERGED) %>%
           arrange(Trait,Tissue,desc(Estimate)) %>%
           write.table(snakemake@output[[2]],sep='\t',row.names = F, col.names = T)



#Table S5



rbind(colSums(outdf[1:49,2:4]),
      colSums(outdf[50:98,2:4]),
      colSums(outdf[99:149,2:4]),
      colSums(outdf[149:196,2:4])) %>%
data.frame() %>% 
mutate(trait = c("CAD","HDL","Heights","LDL")) %>%
                        select(trait,n_GO_terms_tested,n_sig05) %>%
                        write.table(snakemake@output[[1]],sep = '\t',row.names = F,col.names = T)


