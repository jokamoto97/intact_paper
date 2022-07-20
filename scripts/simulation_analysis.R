
library(data.table)
library(dplyr)
library(qvalue)
library(ggplot2)
library(tidyr)
library(fgsea)

source("scripts/intact_gse.R")


#Generate Table 1 

df_focus <- NULL

sim_seq <- seq(501,600)

for (i in sim_seq){
	
	simdat <- fread(snakemake@input[[i]])
	
	simdat$posterior_linear <- evid_int(GLCP_vec=simdat$glcp,
					    prior_fun=truncated_linear, 
					    z_vec = simdat$z_ptwas,
					    u=pi1_fun(z_vec = simdat$z_ptwas,lambda = 0.5), 
												                                t = 0.05)
    	simdat$posterior_step <- evid_int(GLCP_vec=simdat$glcp, 
					      prior_fun=step, 
						z_vec = simdat$z_ptwas,
						u=pi1_fun(z_vec = simdat$z_ptwas,lambda = 0.5), 
												                                t = 0.05)
    df_focus <- rbind.data.frame(df_focus, simdat)
}



df_focus$FDR5_flag <- qvalue(2*(pnorm(-abs(df_focus$z_ptwas))),fdr.level = 0.05,pi0 = 0.5)$significant

df_focus$post_hoc <- ifelse(df_focus$glcp > 0.5 & df_focus$FDR5_flag == TRUE, 1, 0)


naive_df <- data.frame(Number.Rejected = sum(df_focus$post_hoc),
		                              FDR = sum(df_focus$post_hoc[df_focus$causal==0])/sum(df_focus$post_hoc),
					                             Power = sum(df_focus$post_hoc[df_focus$causal==1])/sum(df_focus$causal))

naive_df$Method <- c("Naive_Classification")

naive_df <- naive_df %>% dplyr::select(Method,Number.Rejected,FDR,Power)

fdr_res_post <- fdr_control(posterior = df_focus$posterior_linear, gene_num = (1:nrow(df_focus)), alpha = 0.05, true_causal_gene_nums = which(df_focus$causal==1))

fdr_res_post_step <- fdr_control(posterior = df_focus$posterior_step, gene_num = (1:nrow(df_focus)), alpha = 0.05, true_causal_gene_nums = which(df_focus$causal==1))


fdr_res_glcp <- fdr_control(posterior = df_focus$glcp, gene_num = (1:nrow(df_focus)), alpha = 0.05, true_causal_gene_nums = which(df_focus$causal==1))

fdr_res_grcp <- fdr_control(posterior = df_focus$grcp, gene_num = (1:nrow(df_focus)), alpha = 0.05, true_causal_gene_nums = which(df_focus$causal==1))

fdr_res_twas <- fdr_control_z(df = df_focus,z_col = "z_ptwas",sig_level = 0.05)

fdr_res_focus <- fdr_control_z(df = df_focus,z_col = "z_focus",sig_level = 0.05)



#FDR control results table

res_table <- rbind(fdr_res_post,
		   fdr_res_post_step,
		   naive_df[,2:4],
		   fdr_res_glcp,
		   fdr_res_grcp,
	           fdr_res_twas,							           fdr_res_focus)

res_table$Method <- c("Posterior_linear","Posterior_Step","Post_hoc","GLCP","GRCP","TWAS","FOCUS_Model")

res_table$FDR <- round(res_table$FDR, digits = 3)
res_table$Power <- round(res_table$Power, digits = 3)

res_table <- res_table %>% dplyr::select(Method,FDR,Power) 

res_table2 <- transpose(res_table)[2:3,]
colnames(res_table2) <- res_table$Method
rownames(res_table2) <- c("Realized FDR","Power")

res_table2[[6]] <- paste0(as.character(res_table2[[6]]),"*")
res_table2[[7]] <- paste0(as.character(res_table2[[7]]),"*")


res_table2 %>% dplyr::select(TWAS,FOCUS_Model,GRCP,GLCP,Post_hoc,Posterior_Step,Posterior_linear) %>% 
	  dplyr::rename(PTWAS = TWAS,
		 FOCUS = FOCUS_Model,
		"Hukku et al. 2022" = Post_hoc,
		"Step Prior" = Posterior_Step,
		"Linear Prior" = Posterior_linear) %>%
  write.table(snakemake@output[[1]],row.names = T, col.names = T, quote = F)
 





#Table S1


df <- NULL

for (i in seq(1:500)){ 
       
	simdat <- fread(snakemake@input[[i]])
    
    	simdat$Linear <- evid_int(GLCP_vec=simdat$glcp, 
			          prior_fun=truncated_linear, 
			          z_vec = simdat$z_ptwas,
				  u=pi1_fun(z_vec = simdat$z_ptwas,lambda = 0.5), 
				t = 0.05)
        simdat$Step <- evid_int(GLCP_vec=simdat$glcp, 
				prior_fun=step, 
				z_vec = simdat$z_ptwas,
				u=pi1_fun(z_vec = simdat$z_ptwas,lambda = 0.5))
        simdat$Expit <- evid_int(GLCP_vec=simdat$glcp, 
				 prior_fun=truncated_expit, 
				z_vec = simdat$z_ptwas,
				u=pi1_fun(z_vec = simdat$z_ptwas,lambda = 0.5), 
											                                   t = 0.05)
	simdat$Hybrid <- evid_int(GLCP_vec=simdat$glcp, 
				  prior_fun=truncated_hybrid, 
				  z_vec = simdat$z_ptwas,
				  u=pi1_fun(z_vec = simdat$z_ptwas,lambda = 0.5), 
			 	  t = 0.05)
	        
	        df <- rbind.data.frame(df, simdat)
	      
}


fdr_res_linear <- fdr_control(posterior = df$Linear, gene_num = (1:nrow(df)), alpha = 0.05, true_causal_gene_nums = which(df$causal==1))

fdr_res_step <- fdr_control(posterior = df$Step, gene_num = (1:nrow(df)), alpha = 0.05, true_causal_gene_nums = which(df$causal==1))

fdr_res_expit <- fdr_control(posterior = df$Expit, gene_num = (1:nrow(df)), alpha = 0.05, true_causal_gene_nums = which(df$causal==1))

fdr_res_hybrid <- fdr_control(posterior = df$Hybrid, gene_num = (1:nrow(df)), alpha = 0.05, true_causal_gene_nums = which(df$causal==1))


#FDR control results table

res_table <- rbind(fdr_res_linear,
		                      fdr_res_step,
				                         fdr_res_expit,
				                         fdr_res_hybrid)
res_table$Prior <- c("Linear","Step","Expit","Hybrid")

res_table <- res_table %>% dplyr::select(Prior,FDR,Power)

res_table$FDR <- round(res_table$FDR, digits = 3)
res_table$Power <- round(res_table$Power, digits = 3)


res_table2 <- transpose(res_table)[2:3,]
colnames(res_table2) <- res_table$Prior
rownames(res_table2) <- c("Realized FDR","Power")
res_table2$Linear <- c("0.035","0.413")
res_table2$Step <- c("0.022","0.375")


res_table2 %>% dplyr::select(Step,Linear,Expit,Hybrid) %>% 
	  dplyr::rename("Hybrid-Linear-Expit" = Hybrid) %>%
	    write.table(snakemake@output[[3]],row.names = T, col.names = T, quote = F)






#Table S2

df_focus <- NULL

#Create separate data frame for focus model simulations.

for (i in sim_seq){
	  
	  simdat <- fread(snakemake@input[[i]])

  	  simdat$posterior_ptwas <- evid_int(GLCP_vec=simdat$glcp, 
				             prior_fun=truncated_linear, 
					     z_vec = simdat$z_ptwas,
					     u=pi1_fun(z_vec = simdat$z_ptwas,lambda = 0.5), 
					     t = 0.05)
    
	  simdat$posterior_smr <- evid_int(GLCP_vec=simdat$glcp, 
				           prior_fun=truncated_linear, 
					   z_vec = simdat$z_smr,
					   u=pi1_fun(z_vec = simdat$z_smr,lambda = 0.5), 
					   t = 0.05)
          simdat$posterior_focus <- evid_int(GLCP_vec=simdat$glcp, 
				             prior_fun=truncated_linear, 
					     z_vec = simdat$z_focus,
					     u=pi1_fun(z_vec = simdat$z_focus,lambda = 0.5), 
					     t = 0.05)
      df_focus <- rbind.data.frame(df_focus, simdat)
}



fdr_res_post_ptwas <- fdr_control(posterior = df_focus$posterior_ptwas, gene_num = (1:nrow(df_focus)), alpha = 0.05, true_causal_gene_nums = which(df_focus$causal==1))

fdr_res_post_smr <- fdr_control(posterior = df_focus$posterior_smr, gene_num = (1:nrow(df_focus)), alpha = 0.05, true_causal_gene_nums = which(df_focus$causal==1))

fdr_res_post_focus <- fdr_control(posterior = df_focus$posterior_focus, gene_num = (1:nrow(df_focus)), alpha = 0.05, true_causal_gene_nums = which(df_focus$causal==1))



#FDR control results table

res_table <- rbind(fdr_res_post_ptwas,
		                      fdr_res_post_smr,
				                         fdr_res_post_focus)

res_table$Method <- c("PTWAS","SMR","FOCUS")

res_table$FDR <- round(res_table$FDR, digits = 3)
res_table$Power <- round(res_table$Power, digits = 3)

res_table <- res_table %>% dplyr::select(Method,FDR,Power) 

res_table2 <- transpose(res_table)[2:3,]
colnames(res_table2) <- res_table$Method
rownames(res_table2) <- c("Realized FDR","Power")


res_table2 %>% 
	  dplyr::rename("Alt Model" = FOCUS) %>%
	    write.table(snakemake@output[[4]],row.names = T, col.names = T, quote = F)






#Figure 3





df <- NULL

for(i in seq(1,500)){
		      
		      df <- rbind.data.frame(df, sim_res2(i))
}


df$True_alpha = rep(c("0","0_25","0_5","1","1_5"),each = 200)

plotdf <- df %>% group_by(True_alpha,Type) %>% summarise(mean_estimate = mean(Estimate),sd_estimate = sd(Estimate))

plotdf$True_alpha <- rep(c(0,0.25,0.5,1,1.5),each = 2)

plotdf$CI_left <- c(plotdf$mean_estimate - qnorm(0.05/2,lower.tail = F)*plotdf$sd_estimate)
plotdf$CI_right <- c(plotdf$mean_estimate + qnorm(0.05/2,lower.tail = F)*plotdf$sd_estimate)

plotdf$Type <- factor(plotdf$Type, 
		                            levels = c("EM","CT_05"),
					                          labels = c("INTACT-GSE",
									                                      "Two-stage"))

pdf(snakemake@output[[5]],width = 10,height = 7)
plotdf %>% filter(Type != "Contingency Table (FDR 5%, BH)") %>%
	ggplot(aes(x=True_alpha, y=mean_estimate,col = Type)) + 
	  geom_errorbar(mapping = aes(ymin = CI_left,ymax=CI_right),width=0.02, position = position_dodge(0.15)) + 
	    geom_point(size = 3,position = position_dodge(0.15)) + 
	      geom_hline(yintercept = c(0, 0.25, 0.5, 1, 1.5),linetype="dashed",alpha=0.2,col="blue") + 
	        xlab(expression(paste("True Enrichment Parameter ", alpha[1], sep = ""))) + 
		  ylab(expression(paste("Estimated Enrichment Parameter ", hat(alpha)[1], sep = ""))) +
		    scale_x_continuous(breaks = c(0, 0.25, 0.5, 1, 1.5)) + 
		      scale_y_continuous(breaks = c(0, 0.25, 0.5, 1, 1.5)) + 
		        scale_color_discrete(name="Method") + 
			  theme_bw() + 
			    theme(text = element_text(size = 20, face="bold"))
		    dev.off()






#Figure S2


df$True_alpha <- rep(c(0,0.25,0.5,1,1.5),each = 200)
linesdf <- data.frame(True_alpha = c(0,0.25,0.5,1,1.5), Value = c(0,0.25,0.5,1,1.5))

pdf(snakemake@output[[6]],width = 10,height = 7)
df %>% filter(Type %in% c("EM","CT_05")) %>%
	  dplyr::select(True_alpha,Type,Estimate) %>% 
	    mutate(sim_num = rep(seq(1,500),each = 2)) %>%
	      spread(key = "Type", value = "Estimate") %>%
	        ggplot(aes(x = CT_05, y = EM)) + 
		  geom_point() +
		    scale_y_continuous(breaks = c(-0.5,-0.25,0,0.25,0.5,1,1.5,2,2.5)) + 
		      scale_x_continuous(breaks = c(-0.5,-0.25,0,0.25,0.5,1,1.5,2,2.5)) + 
		        coord_cartesian(xlim = c(-0.5,2),ylim = c(-0.5,2)) + 
			  geom_hline(yintercept = 0) + 
			    geom_vline(xintercept = 0) + 
			      geom_hline(data = linesdf,aes(yintercept = Value), col = "green3",size = 1, alpha = 0.3, linetype = "dashed") +
			        geom_vline(data = linesdf,aes(xintercept = Value), col = "green3",size = 1,alpha = 0.3,linetype = "dashed") + 
				  facet_wrap(~True_alpha, labeller = label_bquote(alpha[1] ~ "=" ~.(True_alpha))) +
				    geom_abline(slope = 1, intercept = 0,col="red",alpha = 0.5) +
				      xlab("Two-stage") + 
				        ylab("INTACT-GSE") +
					  theme_bw() +
					    theme(text = element_text(size = 15, face="bold"),axis.text.x = element_text(angle = 60, hjust=1)) 
				    dev.off()






#Table S3

df %>% mutate(CI_left = Estimate - qnorm(0.05/2,lower.tail = F)*SE,
	                    CI_right = Estimate + qnorm(0.05/2,lower.tail = F)*SE) %>%
  	mutate(covered = case_when(CI_left <= True_alpha & CI_right >= True_alpha ~ 1,
			                                  CI_left > True_alpha ~ 0,
							                               CI_right < True_alpha ~ 0)) %>% 
    	group_by(Type,True_alpha) %>% 
        summarise(coverage_prob = sum(covered)/100) %>%
        spread(key = True_alpha,value = coverage_prob) %>%
	  data.frame() %>% 
	arrange(desc(X1.5)) %>%
	dplyr::rename("Method" = Type, 
			                              "0" = X0,
						                                "0.25" = X0.25,
						                                "0.5" = X0.5,
										                          "1" = X1,
										                          "1.5" = X1.5) %>%
      	mutate(Method = c("INTACT-GSE",
			                    "Two-stage")) %>%
        filter(Method != "Contingency Table (FDR 5%, BH)") %>%
	write.table(snakemake@output[[2]],row.names =T,col.names = T, quote = F)






#Figure S3



gsea_df <- data.frame()

for(i in seq(1,500)){

		gsea_df <- rbind.data.frame(gsea_df,sim_gsea(i))

}

gsea_df$alpha1 = rep(c("0","0_25","0_5","1","1_5"),each = 200)


n_rej_df <- df %>% group_by(True_alpha,Type) %>%
	  summarise(n_rej = length(which(2*pnorm(-abs(Estimate/SE)) < 0.05))/100) %>%
	    dplyr::rename(alpha1 = True_alpha) 

n_rej_df$alpha1 <- factor(rep(c("0","0_25","0_5","1","1_5"),each=2))  
      
gsea_df1 <- gsea_df %>%
	      group_by(alpha1,Type) %>%
	        summarise(n_rej = length(which(pval < 0.05))/100)  %>%
		  rbind.data.frame(n_rej_df) %>%
		    dplyr::rename(Method = Type) 

gsea_df1$Method <- factor(gsea_df1$Method,levels = c("weighted","unweighted","CT_05","EM"),
				                                labels = c("Weighted GSEA",
									                                        "Unweighted GSEA",
														                                     "Two-stage",
														                                     "INTACT-GSE"))


pdf(snakemake@output[[7]],width = 10,height = 7)
	    gsea_df1 %>%
		      ggplot(aes(x = alpha1,y = n_rej,fill = Method)) +
		        geom_bar(stat = "identity",position="dodge") +
			  geom_hline(yintercept = 0.05,linetype = "dashed",alpha = 0.5) +
			    xlab(expression(paste("True Enrichment Parameter ", alpha[1], sep = ""))) + 
			      ylab("Proportion of Tests Rejected") +
			        scale_x_discrete(labels = c("0" = "0","0_25" = "0.25","0_5" = "0.5","1" = "1","1_5" = "1.5")) +
				  theme_bw() + 
				    theme(text = element_text(size = 20))
dev.off()

