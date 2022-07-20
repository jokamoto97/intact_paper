library(data.table)
library(numDeriv)
library(bdsmatrix)
library(SQUAREM)


#prior 1 function (takes GLCP)
expit <- function(GLCP,t=NULL,u=1){
  out <- 1/(1 + exp(-1*(GLCP-0.5)/0.1))
  out <- out * u
  return(out)
}

#prior 2 function (takes GLCP and threshold t)
truncated_expit <- function(GLCP,t=NULL,u=1){
  out <- ifelse(GLCP < t, 0,
                1/(1 + exp(-1*(GLCP-0.5)/0.1)))
  out <- out * u
  return(out)
}

#prior 3 function (takes GLCP)
step <- function(GLCP,t=NULL,u=1){
  out <- ifelse(GLCP < 0.5, 0, 1)
  out <- out * u
  return(out)
}

#prior 4 function (takes a GLCP and threshold t)
truncated_linear <- function(GLCP, t=NULL,u=1){
  out <- ifelse(GLCP > t, GLCP, 0)
  out <- out * u
  return(out)
}

#prior 5 function (takes a GLCP and threshold t)
truncated_hybrid <- function(GLCP, t=NULL,u=1){
  out <- ifelse(GLCP < t, 0, 
                ifelse(GLCP >= t & GLCP < 0.5,  1/(1 + exp(-1*(GLCP-0.5)/0.1)),
                       GLCP))
  out <- out * u
  return(out)
}



evid_int <- function(GLCP_vec, prior_fun, z_vec, t=NULL, u=1, K = c(1,2,4,8,16)){

  #Compute log10(BF) grid from TWAS z score for each K value

  bf_grid_log10 <- sapply(K,FUN = function(K,z_vec){
                  out <- 0.5*log10(1/(1+K)) + (0.5*(z_vec^2) * (K/(1+K)))*log10(exp(1))
                  return(out)
                  },
                  z_vec=z_vec)


  #Perform Bayesian model averaging using log sum exp trick

  bf_log10 <- apply(bf_grid_log10,FUN=function(bf_vec){
    x_star <- max(bf_vec)
    logsumexp <- x_star + log10(sum(10^(bf_vec-x_star)))
    out <- logsumexp - log10(length(bf_vec))
    return(out)
    },
    MARGIN = 1)

  #Compute prior from GLCPs

  prior <- prior_fun(GLCP_vec,t,u)

  #compute log10(posterior odds)

  log10odds <- bf_log10 + log10((prior/(1-prior)))

  #compute posterior probability

  post <- 1/ (1 + 10^(-log10odds))

  return(post)
}


logistic.em <- function(d_vec,pprobs,alpha) {

  #new alpha vector

  alpha_new <- rep(NA,2)

  pi <- mean(pprobs)

  #Induced Bayes factor vector (posterior odds/prior odds), on log10 scale to prevent overflow:

  log10BF <- log10(pprobs/(1-pprobs)) -log10(pi/(1-pi))

  ####E-step####
  pi_t <- 1/(1 + exp(-1*(alpha[1] + alpha[2]*d_vec)))

  #compute log10(posterior odds)

  log10odds <- log10BF + log10((pi_t/(1-pi_t)))

  #Update gamma estimate

  Egamma <- 1 / (1 + 10^(-log10odds))

  #M step

  C00 <- sum((1-Egamma)*(1-d_vec))

  C10 <- sum(Egamma*(1-d_vec))

  C01 <- sum((1-Egamma)*d_vec)

  C11 <- sum(Egamma * d_vec)

  #Update alpha

  alpha0 <- log(C10/C00)

  alpha1 <- log(C00*C11/(C10*C01))

  alpha_new <- c(alpha0,alpha1)

  diff <- sqrt(crossprod(alpha_new-alpha))

  #print(alpha_new)

  alpha <- alpha_new

  return(alpha_new)
}


logistic.em2 <- function(d_vec,pprobs,alpha) {

  #new alpha vector

  alpha_new <- rep(NA,2)

  pi <- mean(pprobs)

  #Induced Bayes factor vector (posterior odds/prior odds), on log10 scale to prevent overflow:

  log10BF <- log10(pprobs/(1-pprobs)) -log10(pi/(1-pi))

  ####E-step####
  pi_t <- 1/(1 + exp(-1*(alpha[1] + alpha[2]*d_vec)))

  #compute log10(posterior odds)

  log10odds <- log10BF + log10((pi_t/(1-pi_t)))

  #Update gamma estimate

  Egamma <- 1 / (1 + 10^(-log10odds))

  #M step

  C00 <- sum((1-Egamma)*(1-d_vec)) + 1/length(d_vec) #add pseudocounts == 1/(#genes)

  C10 <- sum(Egamma*(1-d_vec)) + 1/length(d_vec)

  C01 <- sum((1-Egamma)*d_vec) + 1/length(d_vec)

  C11 <- sum(Egamma * d_vec) + 1/length(d_vec)

  #Update alpha

  alpha0 <- log(C10/C00)

  alpha1 <- log(C00*C11/(C10*C01))

  alpha_new <- c(alpha0,alpha1)

  diff <- sqrt(crossprod(alpha_new-alpha))

  #print(alpha_new)

  alpha <- alpha_new

  return(alpha_new)
}

logistic.loglik <- function(alpha,d_vec,pprobs){

  #Empirical Bayes estimate for prior
  pi <- mean(pprobs)

  #logBF <- log(pprobs/(1-pprobs)) -log(pi/(1-pi))

  #account for genes with pprob == 1 so that likelihood does not become Inf
  pprobs <- ifelse(pprobs > 0, pprobs - 1.e-8, pprobs)

  BF <- (pprobs/(1-pprobs))/(pi/(1-pi))

  liklhd <- (1/(1+exp(-1*(alpha[1] + alpha[2]*d_vec))))*BF +
    1/(1+exp((alpha[1] + alpha[2]*d_vec)))

  loglik <- log(liklhd)

  return(sum(loglik))
}



em_est2<-function(pprobs, d_vec){

    alpha_start <- c(log(mean(pprobs)/(1-mean(pprobs))),0)

    CONVERGED <- TRUE


    #Compute MLEs
      square_obj  <- try(SQUAREM::squarem(par = alpha_start,
                       fixptfn = logistic.em,
                      objfn = logistic.loglik,
                      control = list(tol = 1.e-08, minimize=F, maxiter=50),
                      pprobs=pprobs,
                      d_vec=d_vec),silent=T)
      if("try-error" %in% class(square_obj)){
      try(square_obj <- SQUAREM::squarem(par = alpha_start,
                      fixptfn = logistic.em2,
                      objfn = logistic.loglik,
                      control = list(tol = 1.e-08, minimize=F, maxiter=5),
                      pprobs=pprobs,
                      d_vec=d_vec),silent = T)
          CONVERGED <- FALSE
        if("try-error" %in% class(square_obj)){
          return(c(NA,NA,CONVERGED))
        }
        else{
          return(c(square_obj$par,CONVERGED))
          }
      }else{
          return(c(square_obj$par,CONVERGED))
        }

    #}
}

#function that computes bootstrap standard errors of MLE estimates, given posterior probabilities, gene set annotations, and number of bootstrap samples. Throw out bootstrap samples for which both annotated and non-annotated sets do not have at least one gene with a positive posterior.

enrich_bootstrap_se <- function(pprobs, d_vec, reps=100){

      #create a matrix of sampled posterior probabilities with replacement. Rows are bootstrap samples

      BD = t(sapply(1:reps, function(x) sample(pprobs, length(pprobs), replace = T)) )

      #compute alpha estimates for each bootstap sample

      bootstrap_samples = apply(BD, 1, function(x) em_est(x, d_vec))

      #compute sd across bootstrap estimates

      sdv = apply(bootstrap_samples,1,sd)

 return(sdv)
}

enrich_bootstrap_se2 <- function(pprobs, d_vec, reps = 100){

  #generate bootstrap samples indices

  samples <- replicate(n = 100,
                     expr = sample(seq(length(pprobs)),
                                   size = length(pprobs), replace = TRUE))

  #compute gene set enrichment for each bootstrap sample

  res <- apply(samples, 2, FUN = function(x){
    boot_dat <- cbind(pprobs,d_vec)[x,]
    return(em_est(pprobs = boot_dat[,1], d_vec = boot_dat[,2]))
  })

  #compute SE across gene set enrichment estimates

  sdv = apply(res,1,sd,na.rm = T)

 return(sdv)
}



enrich_res2 <- function(sig_lev, pprobs, d_vec, SE_type = "profile_likelihood", boot_rep = NULL){

  #compute a1 and a0 MLEs using the SQUAREM package

  alpha <- em_est2(pprobs = pprobs, d_vec = d_vec)

  CONVERGED <- alpha[3]

  alpha <- alpha[1:2]

  #alpha0 MLE

  a0 <- alpha[1]

  #alpha1 MLE

  a1 <- alpha[2]

  #compute the specified standard error type

  if (SE_type == "profile_likelihood"){

  #compute a0 standard error

  SE_a0 <- sqrt(a0^2 / (2*(logistic.loglik(alpha = c(a0,a1),
                                                        pprobs = pprobs,
                                                        d_vec = d_vec) -
                             logistic.loglik(alpha = c(0,a1),
                                                        pprobs = pprobs,
                                                        d_vec = d_vec))))


  #compute a1 standard error

  SE_a1 <- sqrt(a1^2 / (2*(logistic.loglik(alpha = c(a0,a1),
                                                        pprobs = pprobs,
                                                        d_vec = d_vec) -
                             logistic.loglik(alpha = c(a0,0),
                                                        pprobs = pprobs,
                                                        d_vec = d_vec))))

  }

  if (SE_type == "bootstrap"){

  SE_vec <- enrich_bootstrap_se2(pprobs = pprobs, d_vec = d_vec, reps = boot_rep)

  SE_a0 <- SE_vec[1]

  SE_a1 <- SE_vec[2]

  }

  if (SE_type == "NDS"){

    if (NA %in% alpha){

      SE_a0 <- SE_a1 <- NA

    }else{

    #Numerical differentiation of log likelihood function to get Hessian evaluated at MLEs

    hess <- numDeriv::hessian(func = logistic.loglik,
                              x = alpha,
                              method = "Richardson",
                              d_vec = d_vec,
                              pprobs = pprobs)

    #Generalized Cholesky decomposition

    cholesk <- bdsmatrix::gchol(-1*hess)

    cov_mat <- solve(cholesk)

    SE_a0 <- sqrt(cov_mat[1,1])

    SE_a1 <- sqrt(cov_mat[2,2])

    }
  }
  if (!(SE_type %in% c("profile_likelihood","bootstrap","NDS"))){

    print("Standard error type must be profile_likelihood, bootstrap, or NDS")

  }

  #assume that likelihood \hat\alpha1|\alpha_1 ~ N(\alpha_1, se^2(\hat alpha_1)) and
  #alpha1 ~ N(0,1)
  # so alpha1|\hat\alpha1 ~ N(1/(1+se^2)alpha1, se^2/(se^2 + 1))

  posterior_mean_a1 <- a1*1/(1+SE_a1^2)

  posterior_var_a1 <- SE_a1^2/(SE_a1^2 + 1)

  posterior_se_a1 <- sqrt(posterior_var_a1)

  #If Hessian fails to converge, use a1 prior as the posterior.

  if (posterior_se_a1 %in% c(0,NA)){

    posterior_se_a1 <- 1

    posterior_mean_a1 <- 0

  }

  #If NDS approach fails, use the prior as the posterior.


  out_df <- data.frame(Parameter = c('a1'),
                    Estimate = c(posterior_mean_a1), #shrink a1 estimate
                    SE = c(posterior_se_a1), #report posterior se for a1
                    z = posterior_mean_a1/posterior_se_a1,
                    pval = 2*pnorm(-abs(posterior_mean_a1/posterior_se_a1)),
                    CI_Leftlim = c(posterior_mean_a1 - qnorm(sig_lev/2,lower.tail = F)*posterior_se_a1),
                    CI_Rightlim = c(posterior_mean_a1 + qnorm(sig_lev/2,lower.tail = F)*posterior_se_a1),
                    CONVERGED = CONVERGED)

  return(out_df)

}




pi1_fun <- function(z_vec,lambda){

  p_vec <- 2*pnorm(abs(z_vec),lower.tail = F)

  p_vec <- p_vec[which(p_vec != 1)]

  pi0 <- length(which(p_vec > lambda))/(length(p_vec)*(1-lambda))

  pi0_max <-  0.99

  pi1 <- 1- min(pi0_max,pi0)

  return(pi1)

}



enrich_GO <- function(go_id, df, prior = truncated_linear, t = 0.05, u = pi1_fun(z_vec = df$TWAS_z,lambda = 0.5)){
  
  #get list of member genes for the GO term
  
  #print(go_id)
  
  go_genes <- as.character(ensembl_go_bp$ensembl_gene_id[ensembl_go_bp$go_id==go_id])
  
  #indicator for membership in the GO term
  
  df$go_term_anno <- 0
  
  df$go_term_anno[df$Gene %in% go_genes] <- 1
  
  #perform evidence integration 
  
  df$posteriors <- evid_int(GLCP_vec=df$GLCP, 
                              prior_fun=prior, 
                              z_vec = df$TWAS_z,
                              u=u, 
                              t = t)
  

  #perform gene set enrichment using posteriors and annotation vector
  
  out <- enrich_res2(sig_lev = 0.05,
            pprobs = df$posteriors, 
            d_vec = df$go_term_anno,
            SE_type = 'NDS')
  
  out$GO_term <- go_id
  
  out$term_size <- length(go_genes)
  
  return(out)
}

#Wrapper function that checks the size of the GO term and proportion of investigated genes represented in df before testing for enrichment. GO term must be between min_size and max_size genes in size and exceed prop_threshold for the proportion of investigated genes.
enrich_go_wrapper <- function(go_id_vec,df,prior = truncated_linear, t = 0.05, u = pi1_fun(z_vec = df$TWAS_z,lambda = 0.5),min_size,max_size,prop_threshold){
  
  res_df <- data.frame()
  
  for (i in 1:length(go_id_vec)){
  
  #create a complete list of genes annotated with the GO term
    
  go_genes <- as.character(ensembl_go_bp$ensembl_gene_id[ensembl_go_bp$go_id==go_id_vec[i]])
  
  #compute number of genes represented in both investigated genes and GO term genes
  
  num_rep <- length(which(go_genes %in% df$Gene))

  #compute proportion of genes represented in both investigated genes and GO term genes
  
  prop_rep <- num_rep/length(go_genes)

  if(length(go_genes) < min_size | length(go_genes) > max_size | prop_rep < prop_threshold){
    next
  }else{
    
    go_res <- enrich_GO(go_id = go_id_vec[i], df = df, prior = prior, t = t, u = u)
    
    go_res$prop_rep <- prop_rep
    
    res_df <- rbind(res_df, go_res)
    
  }
  }
  
  res_df <- res_df[order(res_df$pval),]
  
  rownames(res_df) <- NULL
  
  return(res_df)
}



#takes posteriors, gene index, significance threshold, and indices of true causal genes

fdr_control<-function(posterior,gene_num, alpha=0.05, true_causal_gene_nums){

        lfdr_sort = sort(1-posterior)

        FDR = cumsum(lfdr_sort)/(1:length(lfdr_sort))

        thresh = 1 - lfdr_sort[max(which(FDR<=alpha))]

        rej_gene = as.numeric(gene_num[which(posterior>thresh)])

        rej = length(rej_gene)

        FP = length(which(!(rej_gene %in% true_causal_gene_nums)))

        TP = rej - FP

        out <- data.frame("Number Rejected" = rej, "FDR" = FP/rej, "Power" = TP/length(true_causal_gene_nums))

        return (out)

}

#takes the simdat data frame , z score column name, and significance level

fdr_control_z <- function(df, z_col, sig_level = 0.05){

  df$significant <- qvalue(2*pnorm(abs(df[[z_col]]),lower.tail = F),fdr.level = 0.05)$significant

  rej <- length(which(df$significant == TRUE))

  FP <- length(which(df$significant[df$causal == 0] == TRUE))

  TP <- rej - FP

  out <- data.frame("Number Rejected" = rej, "FDR" = FP/rej, "Power" = TP/length(which(df$causal == 1)))

  return(out)

}

fdr_control_alpha <- function(df, z_col, sig_lev){

  df$significant <- qvalue(2*pnorm(abs(df[[z_col]]),lower.tail = F),fdr.level = sig_lev)$significant

  rej <- length(which(df$significant == TRUE))

  FP <- length(which(df$significant[df$True_alpha=='0']))

  TP <- rej - FP

  out <- data.frame("Number Rejected" = rej, "FDR" = FP/rej, "Power" = TP/length(which(df$True_alpha != "0")))

  return(out)
}


fdr_rst<-function(posterior, alpha=0.05){
	  
	        gene_num <- seq(1,length(posterior))

        lfdr_sort = sort(1-posterior)
	        
	        FDR = cumsum(lfdr_sort)/(1:length(lfdr_sort))
	        
	        thresh = 1 - lfdr_sort[max(which(FDR<=alpha))]
		        
		        rej_gene = as.numeric(gene_num[which(posterior>thresh)])
		        
		        out <- rep(FALSE,length(posterior))
			        
			        out[rej_gene] <- TRUE

			        return(out)
}





sim_res2 <- function(sim_num){

	simdat <- fread(snakemake@input[[sim_num]])
  
        #Compute posterior using PTWAS z score

        simdat$pprobs <- evid_int(GLCP_vec=simdat$glcp, 
				                                    prior_fun=truncated_linear,                         
								                                      z_vec = simdat$z_ptwas,
								                                      u=pi1_fun(z_vec = simdat$z_ptwas,lambda = 0.5), 
												                                        t = 0.05)
	    
	simdat$sig_05_q <- qvalue(2*pnorm(-abs(simdat$z_ptwas)),fdr.level = 0.05,pi0 = 0.5)$significant
	    
	est05_q <- as.vector(coef(summary(glm(data = simdat, sig_05_q ~ annot, family = binomial())))[2,1])
	SE05_q <- as.vector(coef(summary(glm(data = simdat, sig_05_q ~ annot, family = binomial())))[2,2])

	ptwas_res <- enrich_res2(pprobs = simdat$pprobs, d_vec = simdat$annot,sig_lev = 0.05, SE_type = 'NDS',boot_rep = NULL)
			        
			        
			        out_EM <- data.frame("Type" = "EM", "Estimate" = ptwas_res$Estimate,"SE" = ptwas_res$SE)
			        out_CT05 <- data.frame("Type" = "CT_05", "Estimate" = est05_q, "SE" = SE05_q)

				        return(rbind.data.frame(out_EM, out_CT05))
}



sim_gsea <- function(sim_num){

	        set.seed(123)

        dat <- fread(snakemake@input[[sim_num]])

	        gs <- list(as.character(dat$gene[which(dat$annot == 1)]))

	        dat1 <- dat %>% mutate(p.value = 2*pnorm(-abs(z_ptwas))) %>% mutate(z_ptwas = rank(z_ptwas,  ties.method = "random"))

		        score_vec <- dat1$z_ptwas

		        names(score_vec) <- dat1$gene

			        score_vec <- score_vec[order(score_vec,decreasing=F)]

			        names(gs) <- "sim_gs"

				        out1 <- fgsea(pathways = gs, stats = score_vec, nperm = 1000,gseaParam = 0)$pval
				        out2 <- fgsea(pathways = gs, stats = score_vec, nperm = 1000,gseaParam = 1)$pval

					        out <- data.frame("sim_num" = rep(sim_num,2),
								                            "Type" = c("unweighted","weighted"),
											                              "pval" = c(out1,out2))

					        return(out)
}








cad_gse <- function(path_name, tissue){

          print(path_name)

          path_genes <- ctd_genes %>% filter(pathway_id == path_name)

          path_genes <- left_join(path_genes, annotLookup,by = "hgnc_symbol")

          dat_idx <- which(tissues == tissue)

          dat <- fread(snakemake@input[[dat_idx]])

          dat$Posterior <- evid_int(GLCP_vec = dat$GLCP,
                                    prior_fun = truncated_linear,
                                    z_vec = dat$TWAS_z,
                                    u = pi1_fun(z_vec = dat$TWAS_z,lambda = 0.5),
                                    t = 0.05)

          #indicator for membership in the gene set

          dat$anno <- 0

          dat$anno[dat$Gene %in% path_genes$ensembl_gene_id] <- 1

          print(head(dat%>%filter(anno == 1)))

          #perform gene set enrichment using posteriors and annotation vector

          out <- enrich_res2(sig_lev = 0.05,
                                         pprobs = dat$Posterior,
                                                     d_vec = dat$anno,
                                                     SE_type = 'NDS')

          out$path <- path_name

          out$path_def <- as.character(unique(path_genes$pathway_name))

          out$gene_set_size <- length(path_genes$ensembl_gene_id)

          out$prop_rep <- length(which(path_genes$ensembl_gene_id %in% dat$Gene))/length(path_genes$ensembl_gene_id)

          out$tissue <- tissue

          return(out)
}










ev_int_fdr <- function(trait,type,prior){

        if (type == "molecular"){

                res_df <- data.frame()

                for (i in 1:length(tissues)){

                        if(prior == "linear"){

                          trait_idx <- ifelse(trait == 'IGF_1_all',0,
                                             ifelse(trait == "Testosterone_female",49,
                                                    ifelse(trait == "Testosterone_male",49*2,49*3)))

                          dat <- fread(snakemake@input[[i+trait_idx]])

                          dat$Posterior <- evid_int(GLCP_vec = dat$GLCP,
                                                    prior_fun = truncated_linear,
                                                    z_vec = dat$TWAS_z,
                                                    u = pi1_fun(z_vec = dat$TWAS_z,lambda = 0.5),
                                                    t = 0.05)
                        }
                        else{

                          trait_idx <- ifelse(trait == "IGF_1_all",0,
                                             ifelse(trait == "Testosterone_female",49,
                                                 ifelse(trait == "Testosterone_male",49*2,49*3)))

                          dat <- fread(snakemake@input[[i+trait_idx]])

                          dat$Posterior <- evid_int(GLCP_vec = dat$GLCP,
                                                    prior_fun = step,
                                                    z_vec = dat$TWAS_z,
                                                    u = pi1_fun(z_vec = dat$TWAS_z,lambda = 0.5))
                        }

                dat$Tissue <- tissues[i]

                res_df <- rbind.data.frame(res_df,dat)

                }
        }
        else{

              res_df <- data.frame()

              for (i in 1:length(tissues)){

                  if(prior == "linear"){

                          trait_idx <- ifelse(trait == "CAD",0,
                                                ifelse(trait == "HDL",49,
                                                       ifelse(trait == "Heights",2*49,3*49)))
                          dat <- fread(snakemake@input[[i + trait_idx]])

                          dat$Posterior <- evid_int(GLCP_vec = dat$GLCP,
                                                    prior_fun = truncated_linear,
                                                    z_vec = dat$TWAS_z,
                                                    u = pi1_fun(z_vec = dat$TWAS_z,lambda = 0.5),
                         t = 0.05)
                  }

                  else{

                          trait_idx <- ifelse(trait == "CAD",0,
                                           ifelse(trait == "HDL",49,
                                                ifelse(trait == "Heights",2*49,3*49)))

                          dat <- fread(snakemake@input[[i + trait_idx]])

                          dat$Posterior <- evid_int(GLCP_vec = dat$GLCP,
                                                    prior_fun = step,
                                                    z_vec = dat$TWAS_z,
                                                    u = pi1_fun(z_vec = dat$TWAS_z,lambda = 0.5))

              }

                dat$Tissue <- tissues[i]

                res_df <- rbind.data.frame(res_df,dat)

              }
        }

  res_df1 <- res_df %>% group_by(Tissue) %>% mutate(sig = fdr_rst(posterior = Posterior)) %>%   ungroup()
  res_df1$LCP_flag <- ifelse(res_df1$GLCP >=0.5, 1, 0)
  res_df1 <- res_df1 %>% group_by(Tissue) %>% mutate(twas_flag = qvalue(2*pnorm(-abs(TWAS_z)),fdr.level = 0.05,pi0 = 0.5)$significant) %>% ungroup()
  res_df1$post_hoc <- ifelse(res_df1$LCP_flag == 1 & res_df1$twas_flag == TRUE, 1, 0)
  out <- res_df1

    return(out)
}







sig_genes_fun <- function(trait,type){

        dat_linear <- ev_int_fdr(trait,type,prior = "linear")
        dat_step <- ev_int_fdr(trait,type,prior = "step")
        out <- data.frame(trait = trait,
                          post_hoc_n_hits = sum(dat_linear$post_hoc),
                          post_hoc_unique_hits = length(unique(dat_linear$Gene[dat_linear$post_hoc == 1])),
                          step_n_hits = sum(dat_step$sig),
                          step_n_unique_hits = length(unique(dat_step$Gene[dat_step$sig == 1])),
                          linear_n_hits = sum(dat_linear$sig),
                          linear_unique_hits = length(unique(dat_linear$Gene[dat_linear$sig == 1])))

        return(out)
}











sig_genes_list <- function(trait,type){

  dat_linear <- ev_int_fdr(trait,type,prior = "linear")
  out <- unique(dat_linear$Gene[dat_linear$sig == 1])
  out <- data.frame(trait = rep(trait,length(out)), gene = out)

  return(out)
}

