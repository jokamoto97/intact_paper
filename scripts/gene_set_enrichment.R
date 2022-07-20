library(reshape2)
library(stringr)
library(data.table)
library(dplyr)
library(tidyr)
library(qvalue)

source("scripts/intact_gse.R")



ensembl_go_bp <- read.table(snakemake@input[[197]],sep = '\t',header = T)

go_term_list <- as.character(unique(ensembl_go_bp$go_id))

for (i in 1:length(snakemake@input[1:196])){
        print(snakemake@input[[i]])

        real_dat <- fread(snakemake@input[[i]])
	
	real_dat$Posterior <- evid_int(GLCP_vec = real_dat$GLCP,
                                                    prior_fun = truncated_linear,
                                                    z_vec = real_dat$TWAS_z,
                                                    u = pi1_fun(z_vec = real_dat$TWAS_z,lambda = 0.5),
                         t = 0.05)

        real_dat_res <- enrich_go_wrapper(go_id_vec = go_term_list,
                  df = real_dat,
                  min_size = 1,
                  max_size = Inf,
                  prop_threshold = 0.7)
        write.table(real_dat_res,snakemake@output[[i]],sep='\t',quote=F,row.names = F, col.names = T)

}

