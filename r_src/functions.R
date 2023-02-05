#' @param morantest an object from moran.test()
#' @param var character vector of variable
#' @return data.frame with 1 row
tidyMI <- function(morantest, var){
  tidy(morantest) %>% 
    rename(Moran_I_statistic = estimate1, 
           Expectation = estimate2, 
           Variance = estimate3,
           Moran_I_statistic_standard_deviate = statistic) %>% 
    dplyr::select(-method, -alternative) %>% 
    mutate(name = var) %>% 
    relocate(name)
}

collect_3chains_fromCARBayes <- function(get_metrics = TRUE){
  
  out_list <- list()
  
  chain1 <- S.CARleroux(sum_y ~ offset(log(sum_E)) + x + log_sum_E, W = W, data = df_agg, 
                        family = "poisson", thin = 10, prior.var.beta = c(2,2,2),
                        burnin = 10000, n.sample = 50000, verbose = FALSE)
  chain2 <- S.CARleroux(sum_y ~ offset(log(sum_E)) + x + log_sum_E, W = W, data = df_agg, 
                        family = "poisson", thin = 10, prior.var.beta = c(2,2,2),
                        burnin = 10000, n.sample = 50000, verbose = FALSE)
  chain3 <- S.CARleroux(sum_y ~ offset(log(sum_E)) + x + log_sum_E, W = W, data = df_agg, 
                        family = "poisson", thin = 10, prior.var.beta = c(2,2,2),
                        burnin = 10000, n.sample = 50000, verbose = FALSE)
						
  out_list$fitted.values <- rowMeans(matrix(c(chain1$fitted.values,
											  chain2$fitted.values,
											  chain3$fitted.values), 
                  byrow = F, nrow = length(chain1$fitted.values), ncol = 3))
  
  tot_its <- nrow(chain1$samples$beta)
  
  # beta1
  beta1 <- matrix(c(chain1$samples$beta[,1], chain2$samples$beta[,1],
                    chain3$samples$beta[,1]), 
                  byrow = F, nrow = tot_its, ncol = 3)
  
  # beta2
  beta2 <- matrix(c(chain1$samples$beta[,2], chain2$samples$beta[,2],
                    chain3$samples$beta[,2]), 
                  byrow = F, nrow = tot_its, ncol = 3)
				  
  # beta3
  beta3 <- matrix(c(chain1$samples$beta[,3], chain2$samples$beta[,3],
                    chain3$samples$beta[,3]), 
                  byrow = F, nrow = tot_its, ncol = 3)
  
  # rho
  rho <- matrix(c(chain1$samples$rho, chain2$samples$rho,
                  chain3$samples$rho), 
                byrow = F, nrow = tot_its, ncol = 3)
  
  # tau2
  tau2 <- matrix(c(chain1$samples$tau2, chain2$samples$tau2, 
                   chain3$samples$tau2), 
                 byrow = F, nrow = tot_its, ncol = 3)
  
  # phi 
  phi <- rbind(chain1$samples$phi,
               chain2$samples$phi,
               chain3$samples$phi)
  
  ## Summarize all parameters
  out_list$results <- 
    list(
    (dplyr::mutate(as.data.frame(beta1), 
                   it = 1:tot_its) %>% 
       pivot_longer(-it, values_to = "(Intercept)")),
    (dplyr::mutate(as.data.frame(beta2), 
                   it = 1:tot_its) %>% 
       pivot_longer(-it, values_to = "x")),
	(dplyr::mutate(as.data.frame(beta3), 
                   it = 1:tot_its) %>% 
       pivot_longer(-it, values_to = "log_sum_E")),
    (dplyr::mutate(as.data.frame(rho), 
                   it = 1:tot_its) %>% 
       pivot_longer(-it, values_to = "rho")),
    (dplyr::mutate(as.data.frame(tau2), 
                   it = 1:tot_its) %>% 
       pivot_longer(-it, values_to = "tau2"))
  ) %>% 
    reduce(inner_join, by = c("it", "name")) %>% 
    dplyr::select(-name) %>% 
    pivot_longer(-it) %>% 
    dplyr::select(-it) %>% 
    group_by(name) %>% 
    summarise_all(list(mean = mean,
                       median = median,
                       sd = sd,
                       var = var,
                       lower = ~quantile(., probs = 0.025),
                       upper = ~quantile(., probs = 0.975))) %>% 
    mutate(cisize = upper - lower,
           Rhat = c(Rhat(beta1), Rhat(beta2), Rhat(beta3), Rhat(rho), Rhat(tau2)),
           ess_bulk = c(ess_bulk(beta1), ess_bulk(beta2), ess_bulk(beta3), ess_bulk(rho), ess_bulk(tau2)),
           ess_tail = c(ess_tail(beta1), ess_tail(beta2), ess_tail(beta3), ess_tail(rho), ess_tail(tau2))) %>% 
    rename(term = name) %>% 
    setNames(c("term", paste0("lerouxcb_", names(.)[-1])))
  
  if(get_metrics == TRUE){
  
  # Calculate loglikelihood for each iteration
  samples.fitted <- exp(rbind(chain1$samples$beta,
                           chain2$samples$beta,
                           chain3$samples$beta) %*% t(model.matrix(~df_agg$x + df_agg$log_sum_E)) + phi + matrix(log(df_agg$sum_E), byrow = T,
                                                                                              nrow = tot_its * 3, 
                                                                                              ncol = nrow(df_agg)))
  full_its <- tot_its * 3
  samples.loglike <- matrix(NA, byrow = T,
                            nrow = full_its, 
                            ncol = nrow(df_agg))
  for(i in 1:full_its){
    samples.loglike[i,] <- dpois(x=df_agg$sum_y, lambda=samples.fitted[i,], log=TRUE)
  }
  
  # calculate p.w
  p.w <- sum(apply(samples.loglike,2, var), na.rm=TRUE)
    
  # CARBayes calculation of loglikelihood
  regression <- mean(beta1) + mean(beta2) * df_agg$x + apply(phi, 2, mean) + log(df_agg$sum_E) + mean(beta3) * df_agg$log_sum_E
  loglike <- sum( dpois(x=df_agg$sum_y, lambda=exp(regression), log=TRUE) )
  deviance.fitted <- -2 * loglike
  mean.deviance <- -2 * sum(samples.loglike, na.rm=TRUE) /   nrow(samples.loglike)
  p.d <- mean.deviance - deviance.fitted
  
  # Calculate normalized AIC for three chains
      ## AIC which adjusts for sample size
      ## (-2 * log p(y | \hat{\theta} ) + 2 * k)/n
  out_list$leroux_normAIC_pw <- (-2*loglike + 2 * p.w)/nrow(df_agg)
  out_list$leroux_normAIC_pd <- (-2*loglike + 2 * p.d)/nrow(df_agg)
  
  }
  
  # return list
  return(out_list)
  
}