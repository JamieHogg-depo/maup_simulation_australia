##
# Produce result plots
##

rm(list=ls())
library(tidyverse)
library(patchwork)
library(psych)
library(sf)
library(maptools)
library(gridExtra)
library(latex2exp)
options(scipen = 999)

#cur_date <- "20220329"
cur_date <- "20220606"
cur_loc <- paste0(getwd(), "/Zonations/Results/", cur_date)
plot_loc <- paste0(cur_loc, "/plots")
rds_loc <- paste0(cur_loc, "/rds")

## Import rds files ## ---------------------------------------------------------
rds_files <- gsub(".rds", "", list.files(rds_loc))
for(i in 1:length(rds_files)){
  object_name <- gsub("simlist_", "", rds_files[i])
  object_name <- gsub("storeddatalist_", "", object_name)
  assign(object_name, readRDS(paste0(rds_loc, "/", rds_files[i], ".rds")))
}

## Condense ABS model results ## -----------------------------------------------

abs_fits <- list((
  sa1models %>% rename(alpha = y_av, gamma = rho) %>% 
  mutate(t = 0)
  ),
  (
    sa2models %>% rename(alpha = y_av, gamma = rho) %>% 
      mutate(t = 2)
    ),
  (
    sa3models %>% rename(alpha = y_av, gamma = rho) %>% 
      mutate(t = 4)
  ),
  (
    sa4models %>% rename(alpha = y_av, gamma = rho) %>% 
      mutate(t = 6)
  )) %>% 
  reduce(bind_rows) %>% 
  mutate(t = as.factor(t)) %>% 
  rename(Leroux = lerouxcb_median,
         Poisson = nonspatial_estimate,
         MPoisson = rand_estimate)
  


## Setup objects for plotting ## -----------------------------------------------

# convert target pop to factor variables
concor <- data.frame(target_pop = unique(models$target_pop),
                     t = as.factor(1:7))
models <- left_join(models, concor, by = "target_pop") %>% 
  rename(alpha = y_av,
         gamma = rho)
MI <- left_join(MI, concor, by = "target_pop") %>% 
  rename(alpha = y_av,
         gamma = rho)
df_agg_nogeo <- left_join(df_agg_nogeo, concor, by = "target_pop") %>% 
  rename(alpha = y_av,
         gamma = rho)
df_agg <- left_join(df_agg, concor, by = "target_pop") %>% 
  rename(alpha = y_av,
         gamma = rho)
data_zone_metrics <- left_join(data_zone_metrics, concor, by = "target_pop") %>% 
  rename(alpha = y_av,
         gamma = rho)

# explore the median and IQR of zone sizes
df_agg_nogeo %>% 
  group_by(target_pop, t) %>% 
  summarise(n = n(),
            Median = median(sum_pop),
            q1 = quantile(sum_pop, prob = 0.25),
            q3 = quantile(sum_pop, prob = 0.75)) 
sa1data %>% 
  filter(rho == 0.5, y_av == 5) %>% 
  st_drop_geometry() %>% 
  mutate(sa2 = str_sub(SA1_MAIN16, 1, 3)) %>% 
  group_by(sa2) %>% 
  summarise(pop = sum(pop)) %>% 
  summarise(Median = median(pop),
            q1 = quantile(pop, prob = 0.25),
            q3 = quantile(pop, prob = 0.75))

####
## PLOTS
####

## Convergence of leroux model ----
# average ess and rhat
models %>% 
  filter(term == "x") %>%
  summarise_each(funs(min,median,max), lerouxcb_Rhat, lerouxcb_ess_bulk) %>% 
  pivot_longer(everything())

# average Rhat for beta
models %>% 
  filter(term == "x") %>%
  group_by(t) %>% 
  summarise_each(funs(min,median,max), lerouxcb_ess_bulk, lerouxcb_Rhat)

# plots of convergence
models %>% 
  ggplot(aes(y = lerouxcb_ess_bulk, fill = term)) +
  theme_bw()+
  geom_boxplot()+
  labs(title = "Bulk ESS of parameters from Leroux model",
       y = "Bulk ESS",
       fill = "Parameter")
  ggsave(file = paste0(plot_loc, "/ESS.png"),
         width = 11.3, height = 6.76)
  
models %>% 
  group_by(term) %>% 
  summarize(min = min(lerouxcb_ess_bulk),
            q1 = quantile(lerouxcb_ess_bulk, 0.25),
            mean = mean(lerouxcb_ess_bulk),
            median = median(lerouxcb_ess_bulk),
            q3 = quantile(lerouxcb_ess_bulk, 0.75),
            max = max(lerouxcb_ess_bulk))

models %>% 
  ggplot(aes(y = lerouxcb_Rhat, fill = term)) +
  theme_bw()+
  geom_boxplot()+
  labs(title = "Rhat of parameters from Leroux model",
       y = "Rhat",
       fill = "Parameter")
ggsave(file = paste0(plot_loc, "/Rhat.png"),
       width = 11.3, height = 6.76)

models %>% 
  group_by(term) %>% 
  summarize(min = min(lerouxcb_Rhat),
            q1 = quantile(lerouxcb_Rhat, 0.25),
            mean = mean(lerouxcb_Rhat),
            median = median(lerouxcb_Rhat),
            q3 = quantile(lerouxcb_Rhat, 0.75),
            max = max(lerouxcb_Rhat))

## Leroux vs poisson - point estimates ----
models %>% 
  filter(term == "x") %>% 
  mutate(poisson_bigger = nonspatial_estimate > lerouxcb_median) %>% 
  ggplot(aes(y = nonspatial_estimate,
             x = lerouxcb_median, 
             col = as.factor(poisson_bigger)))+
  geom_point()+theme_bw()+
  geom_abline(slope = 1, intercept = 0)+
  facet_grid(t~gamma, labeller = purrr::partial(label_both, sep = " = "))+
  labs(title = "Poisson point estimates vs Leroux posterior medians",
       col = "Poisson > Leroux",
       x = "Leroux model (Median of posterior)",
       y = "Poisson model (MLE point estimate)")
ggsave(file = paste0(plot_loc, "/leroux_vs_poisson_x.png"),
       width = 11.3, height = 6.76)

## Boxplot: zoning distributions for both models ----- covariate x

grrr <- expand.grid(gamma = c(0.05, 0.5, 0.95),
                    alpha =c(5, 25, 50)) %>% 
  mutate(Scenario = paste0("S", 1:9))
test1 <- models %>%
  filter(term == "x") %>% 
  dplyr::select(lerouxcb_median, nonspatial_estimate, rand_estimate, t, alpha, gamma) %>% 
  pivot_longer(-c(t, alpha, gamma)) %>% 
  mutate(name = factor(name, labels = c("Leroux", "Baseline", "MPoisson"))) %>% 
  filter(name != "MPoisson")
test2 <- abs_fits %>% 
  filter(term == "x") %>% 
  dplyr::select(Leroux, Poisson, MPoisson, t, alpha, gamma) %>% 
  pivot_longer(-c(t, alpha, gamma)) %>% 
  rename(abs = value) %>% 
  mutate(name = ifelse(name == "Poisson", "Baseline", name)) %>% 
  filter(name != "MPoisson")
ggplot()+
  theme_bw()+
  geom_boxplot(data = test1, aes(y = value, x = t, col = as.factor(name)),
               alpha = 0.4) +
  geom_point(data = test2, aes(y = abs, x = t, col = as.factor(name)), 
             shape = 15, size = 3) + 
  geom_hline(data = filter(test2, t == "0"), 
             aes(yintercept = abs, col = as.factor(name)), 
             linetype = "dotted")+
  facet_grid(alpha~gamma, labeller = purrr::partial(label_both, sep = " = "))+
  labs(col = "Model",
       title = "Zoning distributions for covariate effect",
       x = "Target Population",
       y = "Covariate point estimate")+
  geom_text(data = grrr, 
            mapping = aes(x = 2, y = 1.3, label = Scenario),
            hjust = -0.1,
            vjust = -1)
ggsave(file = paste0(plot_loc, "/boxplot_x.png"),
       width = 11.3, height = 6.76)

## Compare points estimates ## -------------------------------------------------
zz <- 6
(models %>% 
  filter(t == zz) %>% 
  ggplot(aes(y = nonspatial_estimate, x = lerouxcb_mean))+
  theme_bw()+
  geom_point()+
  facet_grid(alpha~gamma, labeller = purrr::partial(label_both, sep = " = "))+
  geom_abline(slope = 1, intercept = 0)+
    labs(title = "Poisson vs Leroux"))/

(models %>% 
  filter(t == zz) %>% 
  ggplot(aes(y = rand_estimate, x = lerouxcb_mean))+
  theme_bw()+
  geom_point()+
  facet_grid(alpha~gamma, labeller = purrr::partial(label_both, sep = " = "))+
  geom_abline(slope = 1, intercept = 0)+
   labs(title = "Mixed poisson vs Leroux"))

models %>% 
  group_by(t, alpha, gamma) %>% 
  summarise(cor_pl = cor(nonspatial_estimate, lerouxcb_mean, use = "complete.obs"),
            cor_rl = cor(rand_estimate, lerouxcb_mean, use = "complete.obs")) %>% view

## Boxplot: zoning distributions for both models -- intercept
test1 <- models %>%
  filter(term == "(Intercept)") %>% 
  dplyr::select(lerouxcb_median, nonspatial_estimate, rand_estimate, t, alpha, gamma) %>% 
  pivot_longer(-c(t, alpha, gamma)) %>% 
  mutate(name = factor(name, labels = c("Leroux", "Baseline", "MPoisson"))) %>% 
  filter(name != "MPoisson")
test2 <- abs_fits %>% 
  filter(term == "(Intercept)") %>% 
  dplyr::select(Leroux, Poisson, MPoisson, t, alpha, gamma) %>% 
  pivot_longer(-c(t, alpha, gamma)) %>% 
  rename(abs = value) %>% 
  mutate(name = ifelse(name == "Poisson", "Baseline", name)) %>% 
  filter(name != "MPoisson")
ggplot()+
  theme_bw()+
  geom_boxplot(data = test1, aes(y = value, x = t, col = as.factor(name)),
               alpha = 0.4) +
  geom_point(data = test2, aes(y = abs, x = t, col = as.factor(name)), 
             shape = 15, size = 3) + 
  geom_hline(data = filter(test2, t == "0"), 
             aes(yintercept = abs, col = as.factor(name)), 
             linetype = "dotted")+
  facet_grid(alpha~gamma, labeller = purrr::partial(label_both, sep = " = "))+
  labs(col = "Model",
       title = "Zoning distributions for intercept",
       x = "Target Population",
       y = "Intercept point estimate")+
  geom_text(data = grrr, 
            mapping = aes(x = 2, y = 1.3, label = Scenario),
            hjust = -0.1,
            vjust = -1)
ggsave(file = paste0(plot_loc, "/boxplot_intercept.png"),
       width = 11.3, height = 6.76)

## Boxplot: zoning distributions for both models -- log_sum_E
test1 <- models %>%
  filter(term == "log_sum_E") %>% 
  dplyr::select(lerouxcb_median, nonspatial_estimate, rand_estimate, t, alpha, gamma) %>% 
  pivot_longer(-c(t, alpha, gamma)) %>% 
  mutate(name = factor(name, labels = c("Leroux", "Poisson", "MPoisson")))
test2 <- abs_fits %>% 
  filter(term == "log_sum_E") %>% 
  dplyr::select(Leroux, Poisson, MPoisson, t, alpha, gamma) %>% 
  pivot_longer(-c(t, alpha, gamma)) %>% 
  rename(abs = value)
ggplot()+
  theme_bw()+
  geom_boxplot(data = test1, aes(y = value, x = t, col = as.factor(name)),
               alpha = 0.4) +
  geom_point(data = test2, aes(y = abs, x = t, col = as.factor(name)), 
             shape = 15, size = 3) + 
  geom_hline(data = filter(test2, t == "0"), 
             aes(yintercept = abs, col = as.factor(name)), 
             linetype = "dotted")+
  facet_grid(alpha~gamma, labeller = purrr::partial(label_both, sep = " = "))+
  labs(col = "Model",
       title = "Zoning distributions for coefficient of log_sum_E",
       x = "Target Population",
       y = "Coefficient point estimate")
ggsave(file = paste0(plot_loc, "/boxplot_log_sum_E.png"),
       width = 11.3, height = 6.76)

## Line plots -----
temp <- models %>% 
  filter(term == "x") %>% 
  group_by(alpha, gamma, t) %>% 
  summarise(n = n(),
            Poisson_mean = mean(nonspatial_estimate),
            Poisson_sd = sd(nonspatial_estimate),
            Leroux_mean = mean(lerouxcb_median),
            Leroux_sd = sd(lerouxcb_median),
            .groups = "drop") %>% 
  mutate(Poisson_u = Poisson_mean + 1.96 * Poisson_sd,
         Poisson_l = Poisson_mean - 1.96 * Poisson_sd,
         Leroux_u = Leroux_mean + 1.96 * Leroux_sd,
         Leroux_l = Leroux_mean - 1.96 * Leroux_sd) %>% 
  arrange(t)
just_t1 <- filter(temp, t == 1) %>% 
  dplyr::select(-t) %>% dplyr::select(alpha, gamma, Leroux_sd, Poisson_sd) %>% 
  rename(lsd = Leroux_sd,
         psd = Poisson_sd)
# add ASGS aggregation results
abs_points <- abs_fits %>% 
  filter(term == "x") %>% 
  rename(absLeroux_mean = Leroux,
         absPoisson_mean = Poisson) %>% 
  dplyr::select(t, alpha, gamma, absPoisson_mean, absLeroux_mean)

# pivot data to correct format
temp <- temp %>% 
  left_join(.,just_t1, by = c("alpha", "gamma")) %>% 
  mutate(Leroux_ratio1 = Leroux_sd/lsd,
         Poisson_ratio1 = Poisson_sd/psd) %>% 
  dplyr::select(-lsd, -psd) %>% 
  group_by(alpha, gamma) %>% 
  mutate(Leroux_ratiosd = Poisson_sd/lag(Poisson_sd),
         Poisson_ratiosd = Leroux_sd/lag(Leroux_sd)) %>% 
  ungroup() %>% 
  full_join(.,abs_points, by = c("t", "alpha", "gamma")) %>% 
  pivot_longer(-c(gamma, alpha, t, n), 
               names_to = c("model", "metric"),
               names_pattern = "(.*)_(.*)") %>% 
  pivot_wider(values_from = value,
              names_from = metric) %>% 
  mutate(t = fct_relevel(t, "0"),
         model = as.character(model))

## both models - line plots -----
temp %>% 
ggplot(aes(y = mean, 
                    ymin = l,
                    ymax = u, 
                    x = as.numeric(t) - 1,
                    col = model))+
  theme_bw()+ geom_line()+geom_point()+ geom_errorbar()+
  facet_grid(alpha~gamma, labeller = purrr::partial(label_both, sep = " = ")) + 
  labs(x = "Target population",
       y = "Point estimates for covariate",
       title = "Zoning distribution with 95% confidence intervals", 
       col = "Model")
ggsave(file = paste0(plot_loc, "/lineplot_x.png"),
       width = 11.3, height = 6.76)

## both models - lagged relative increase in sd -----
temp %>% 
  filter(model != "absLeroux") %>% 
  filter(model != "absPoisson") %>% 
ggplot(aes(y = ratiosd, 
                 x = as.numeric(t),
                 col = as.factor(model)))+
  geom_line()+geom_point()+theme_bw()+
  facet_grid(alpha~gamma, labeller = purrr::partial(label_both, sep = " = ")) + 
  labs(x = "Target population",
       y = "Ratio of zoning standard deviations",
       title = "Relative increase in zoning standard deviations", 
       col = "Model", 
       subtitle = "Each point is relative to the previous target population")
ggsave(file = paste0(plot_loc, "/lineplot_ratiosd.png"),
       width = 11.3, height = 6.76)

## both models - lagged relative increase in sd to target pop 1 -----
ggplot(temp, aes(y = ratio1, 
                 x = as.numeric(t),
                 col = as.factor(model)))+
  geom_line()+geom_point()+theme_bw()+
  facet_grid(alpha~gamma, labeller = purrr::partial(label_both, sep = " = ")) + 
  labs(x = "Target population",
       y = "Ratio of zoning standard deviations",
       title = "Relative increase in zoning standard deviations", 
       col = "Model", 
       subtitle = "Each point is relative to t = 1")
ggsave(file = paste0(plot_loc, "/lineplot_ratio1.png"),
       width = 11.3, height = 6.76)


## Boxplot: standard deviate in Morans I - raw SIR ----
MI %>%
  filter(name == "raw_SIR") %>% 
  ggplot(aes(y = Moran_I_statistic_standard_deviate,
             x = t, 
             #fill = t, 
             group = t))+
  theme_bw()+
  geom_boxplot() +
  facet_grid(alpha~gamma, labeller = purrr::partial(label_both, sep = " = "))+
  labs(x = "Target Population",
       title = "Standard Deviate of Moran's I for raw incidence ratio",
       y = "Moran's I Standard Deviate")+
  geom_hline(yintercept = 0, col = "black", linetype = "dashed")+
  geom_hline(yintercept = 1.96, col = "black", linetype = "dotted")+
  geom_hline(yintercept = -1.96, col = "black", linetype = "dotted")
ggsave(file = paste0(plot_loc, "/MI_standard_deviate_rawSIR.png"),
       width = 11.3, height = 6.76)

## Boxplot: standard deviate in Morans I - y ----
MI %>%
  filter(name == "y") %>% 
  ggplot(aes(x = Moran_I_statistic_standard_deviate, 
             fill = t, 
             group = t))+
  theme_bw()+
  geom_boxplot() +
  facet_grid(alpha~gamma, labeller = purrr::partial(label_both, sep = " = "))+
  labs(fill = "Target Population",
       title = "Standard Deviate of Moran's I for y",
       x = "Moran's I Standard Deviate")+
  geom_vline(xintercept = 0, col = "black", linetype = "dashed")+
  geom_vline(xintercept = 1.96, col = "black", linetype = "dotted")+
  geom_vline(xintercept = -1.96, col = "black", linetype = "dotted")+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
ggsave(file = paste0(plot_loc, "/MI_standard_deviate_y.png"),
         width = 11.3, height = 6.76)

## Compare zoning distribution for rho for target_pop = 1 ----

models %>% 
  filter(t == 1, term == "rho" | term == "tau2") %>% 
  ggplot(aes(y = lerouxcb_median, 
             x = as.factor(alpha),
             fill = as.factor(gamma)))+
  theme_bw()+
  geom_boxplot()+
  facet_grid(~term, labeller = purrr::partial(label_both, sep = " = "))+
  labs(title = "Zoning distribution for rho in Leroux model",
       subtitle = "Target population of 5,000",
       fill = "gamma",
       x = "alpha",
       y = "Median of posterior for rho")
ggsave(file = paste0(plot_loc, "/rhotau2_postmed_targetpop1.png"),
         width = 11.3, height = 6.76)

models %>% 
  filter(t == 1, term == "rho" | term == "tau2") %>% 
  ggplot(aes(y = lerouxcb_cisize, 
             x = as.factor(alpha),
             fill = as.factor(gamma)))+
  theme_bw()+
  geom_boxplot()+
  facet_grid(~term, labeller = purrr::partial(label_both, sep = " = "))+
  labs(title = "Zoning distribution for rho in Leroux model",
       subtitle = "Target population of 5,000",
       fill = "gamma",
       x = "alpha",
       y = "Size of credible interval for rho")
ggsave(file = paste0(plot_loc, "/rhotau2_cisize_targetpop1.png"),
       width = 11.3, height = 6.76)

models %>% 
  filter(t == 2, term == "rho" | term == "tau2") %>% 
  ggplot(aes(y = lerouxcb_median, 
             x = as.factor(alpha),
             fill = as.factor(gamma)))+
  theme_bw()+
  geom_boxplot()+
  facet_grid(~term, labeller = purrr::partial(label_both, sep = " = "))+
  labs(title = "Zoning distribution for rho in Leroux model",
       subtitle = "Target population of 10,000",
       fill = "gamma",
       x = "alpha",
       y = "Median of posterior for rho")
ggsave(file = paste0(plot_loc, "/rhotau2_postmed_targetpop2.png"),
       width = 11.3, height = 6.76)

models %>% 
  filter(t == 2, term == "rho" | term == "tau2") %>% 
  ggplot(aes(y = lerouxcb_cisize, 
             x = as.factor(alpha),
             fill = as.factor(gamma)))+
  theme_bw()+
  geom_boxplot()+
  facet_grid(~term, labeller = purrr::partial(label_both, sep = " = "))+
  labs(title = "Zoning distribution for rho in Leroux model",
       subtitle = "Target population of 10,000",
       fill = "gamma",
       x = "alpha",
       y = "Size of credible interval for rho")
ggsave(file = paste0(plot_loc, "/rhotau2_cisize_targetpop2.png"),
       width = 11.3, height = 6.76)

## Median mean error across all t, alpha and gamma ## -----

df_agg_nogeo %>% 
  mutate(rel_mean_abs_error = mean_abs_error/mean_y) %>% 
  group_by(alpha, gamma, t, zone_seed) %>% 
  summarise(MRMAE = median(rel_mean_abs_error),
            .groups = "drop") %>% 
  mutate(alpha = as.factor(alpha),
         gamma = as.factor(gamma)) %>% 
  ggplot(aes(y = MRMAE, x = t)) +
  geom_boxplot(alpha = 0.4)+
    facet_grid(alpha~gamma, labeller = purrr::partial(label_both, sep = " = ")) + #scales = "free_y")+
  theme_bw()+
  labs(title = "Relative mean absolute error for different zonations",
       y = "Median of RMAE",
       x = "Target Population")
ggsave(file = paste0(plot_loc, "/MRMAE.png"),
       width = 11.3, height = 6.76)

## Median Variance of y's in each aggregated region ----
## Across all t, alpha, gamma 

df_agg_nogeo %>% 
  group_by(alpha, gamma, t, zone_seed) %>% 
  summarise(M_var_y = median(var_y),
            n = n(),
            .groups = "drop") %>% 
  mutate(alpha = as.factor(alpha),
         gamma = as.factor(gamma)) %>% 
  ggplot(aes(y = M_var_y, x = t)) +
  geom_boxplot(alpha = 0.4)+
  facet_grid(alpha~gamma, labeller = purrr::partial(label_both, sep = " = "),
             scales = "free_y")+
  theme_bw()+
  labs(title = "Variance of SA1s for new zones",
       y = "Median of variance of y",
       x = "Target Population")
ggsave(file = paste0(plot_loc, "/vary.png"),
       width = 11.3, height = 6.76)

## Mean y's in each zonation ----
## Across all t, alpha, gamma 

df_agg_nogeo %>% 
  group_by(alpha, gamma, t, zone_seed) %>% 
  summarise(mean_y = median(mean_y),
            n = n(),
            .groups = "drop") %>% 
  mutate(alpha = as.factor(alpha),
         gamma = as.factor(gamma)) %>% 
  ggplot(aes(y = mean_y, x = t)) +
  geom_boxplot(alpha = 0.4)+
  facet_grid(alpha~gamma, labeller = purrr::partial(label_both, sep = " = "),
             scales = "free_y")+
  theme_bw()+
  labs(title = "Mean of SA1s for new zones",
       y = "Median of mean of y",
       x = "Target Population")
ggsave(file = paste0(plot_loc, "/mean.png"),
       width = 11.3, height = 6.76)

## Compare model fits ## -------------------------------------------------------

data_zone_metrics %>% 
  ggplot(aes(y = poisson_normAIC, x = t))+
  geom_boxplot(alpha = 0.4)+
  facet_grid(alpha~gamma, labeller = purrr::partial(label_both, sep = " = "),
             scales = "free_y")+
  theme_bw()+
  labs(y = "Normalized AIC - Poisson model", 
       x = "Target Population",
       col = "", 
       title = "Normalized AIC for all Poisson model fits")
ggsave(file = paste0(plot_loc, "/poissonAIC.png"),
       width = 11.3, height = 6.76)

data_zone_metrics %>% 
  ggplot(aes(y = leroux_normAIC_pw, x = t))+
  geom_boxplot(alpha = 0.4)+
  facet_grid(alpha~gamma, labeller = purrr::partial(label_both, sep = " = "))+
  theme_bw()+
  labs(y = "Normalized AIC - Leroux model", 
       x = "Target Population",
       col = "", 
       title = "Normalized AIC for all leroux model fits")
ggsave(file = paste0(plot_loc, "/lerouxAIC_pw.png"),
       width = 11.3, height = 6.76)

## Create subset maps for simulated data ## ------------------------------------

# choose one map combination
map <- sa1data %>% 
  filter(rho == 0.05 & y_av == 50) %>% 
  st_transform(crs = st_crs(3112))

# map border
map.border <- st_union(map)

# we must use crs 3112
# https://epsg.io/map#srs=3112&x=1579959.78&y=-3939741.29&z=10&layer=streets
lims <- list(Brisbane = list(xlim = c(1832429.78, 1861917.08),
                             ylim = c(-3236788.81, -3269337.58)),
             Sydney = list(xlim = c(1554778.93, 1584318.50),
                           ylim = c(-3898068.10, -3940799.69))
)

# Trim extreme values
cut.offs <- c(1/1.5, 1.5)
# SIR.clip <- apply(SIR, 2, median) # Finds the median SIR from MCMC
map$raw_SIR[which(map$raw_SIR < cut.offs[1], arr.ind = TRUE)] <- cut.offs[1]
map$raw_SIR[which(map$raw_SIR > cut.offs[2], arr.ind = TRUE)] <- cut.offs[2]

# Fill colours
Fill.colours <- c("#2C7BB6", "#2C7BB6", "#ABD9E9", "#FFFFBF", "#FDAE61", "#D7191C", "#D7191C")
End <- log2(1.6)
Breaks.fill <- c(1/1.5, 1/1.25, 1, 1.25, 1.5)
Fill.values <- c(-End, log2(Breaks.fill), End)
Fill.values.r <- scales::rescale(x = Fill.values, from = range(Fill.values, na.rm = T))

# Plots without transparency
gg.base <-  ggplot(data = map) +
  geom_sf(aes(fill = log2(raw_SIR)), colour = NA) +
  geom_sf(data = map.border, colour = "black", fill = NA)+
  scale_fill_gradientn("", colours = Fill.colours, values = Fill.values.r ,
                       labels = as.character(round(Breaks.fill, 3)),
                       breaks = log2(Breaks.fill), limits = range(Fill.values)) +
  theme_void()+
  theme(
    legend.position = "bottom",
    plot.title = element_text(margin = margin(0,0,2,0)),
    plot.margin = unit(c(1,1,1,1), "mm")
  ) +
  guides(fill = guide_colourbar(barwidth = 20)
  ) +
  coord_sf()
gg.base

# Brisane inset
brisbane_inset <- gg.base +
  guides(fill = "none")+
  coord_sf(xlim = lims$Brisbane$xlim,
           ylim = lims$Brisbane$ylim)

# Sydney inset
sydney_inset <- gg.base +
  guides(fill = "none")+
  coord_sf(xlim = lims$Sydney$xlim,
           ylim = lims$Sydney$ylim)

all <- grid.arrange(grobs = list(gg.base, brisbane_inset, sydney_inset), nrow = 1)

## Create a plot of a single grid of simulated areas ## ------------------------

df_agg <- readRDS(paste0(rds_loc, "/simlist_df_agg.rds")) %>% 
  rename(alpha = y_av, gamma = rho) %>% 
  left_join(., concor, by = "target_pop")
df_agg_MIN <- st_drop_geometry(df_agg)

lims <- list(Brisbane = list(xlim = c(1832429.78, 1861917.08),
                             ylim = c(-3236788.81, -3269337.58)),
             Sydney = list(xlim = c(1554778.93, 1584318.50),
                           ylim = c(-3898068.10, -3940799.69))
)

sa_codes <- data.frame(t = as.factor(1:7),
                       sa = paste0("SA", seq(1.5,4.5,0.5)))

# Trim extreme values
cut.offs <- c(1/4, 4)
# SIR.clip <- apply(SIR, 2, median) # Finds the median SIR from MCMC
df_agg_MIN$raw_SIR[which(df_agg_MIN$raw_SIR < cut.offs[1], arr.ind = TRUE)] <- cut.offs[1]
df_agg_MIN$raw_SIR[which(df_agg_MIN$raw_SIR > cut.offs[2], arr.ind = TRUE)] <- cut.offs[2]

# Fill colours
Fill.colours <- c("#2C7BB6", "#2C7BB6", "#ABD9E9", "#FFFFBF", "#FDAE61", "#D7191C", "#D7191C")
End <- log2(4.1)
Breaks.fill <- c(1/4, 1/2, 1/1.5, 1/1.25, 1, 1.25, 1.5, 2, 4)
Fill.values <- c(-End, log2(Breaks.fill), End)
Fill.values.r <- scales::rescale(x = Fill.values, from = range(Fill.values, na.rm = T))

# Create plot
df_agg %>% 
  filter(gamma == 0.95 & 
           alpha ==  25 & 
           zone_seed %in% c(1,2)) %>%
  mutate(zone_seed = paste0("Seed ", zone_seed)) %>% 
  left_join(.,sa_codes, by = "t") %>% 
  st_transform(crs = st_crs(3112)) %>% 
  ggplot() +
  geom_sf(aes(fill = log2(raw_SIR)))+
  coord_sf(xlim = lims$Brisbane$xlim,
           ylim = lims$Brisbane$ylim)+
  theme_bw()+
  scale_fill_gradientn("", colours = Fill.colours, values = Fill.values.r ,
                       labels = as.character(round(Breaks.fill, 3)),
                       breaks = log2(Breaks.fill), limits = range(Fill.values))+ 
  theme(
    legend.position = "bottom",
    plot.title = element_text(margin = margin(0,0,2,0)),
    plot.margin = unit(c(1,1,1,1), "mm"),
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank()#,
    #strip.text.x = element_blank(),
    #strip.text.y = element_blank()
  ) +
  guides(fill = guide_colourbar(barwidth = 20)
  ) +
  labs(fill = "",
       title = "Raw incidence ratio for simulated data",
       subtitle = TeX(r'($\gamma = 0.95, \alpha = 25$)'))+
  facet_grid(zone_seed ~ sa)

ggsave(file = paste0(plot_loc, "/brisbane_RIR.png"),
       width = 7.3, height = 3.76)


## Create a series of Brisbane maps for a particular scenario ## ---------------

pos_grid <- expand.grid(gamma = c(0.05, 0.5, 0.95),
                        alpha =c(5, 25, 50))
lims <- list(Brisbane = list(xlim = c(1832429.78, 1861917.08),
                             ylim = c(-3236788.81, -3269337.58)),
             Sydney = list(xlim = c(1554778.93, 1584318.50),
                           ylim = c(-3898068.10, -3940799.69))
)

for(i in 1:nrow(pos_grid)){
  
  plot_name <- paste0("gamma", pos_grid$gamma[i],
                      "alpha",  pos_grid$alpha[i])
  
  df_agg %>% 
    filter(gamma ==  pos_grid$gamma[i] & 
             alpha ==  pos_grid$alpha[i] & 
             zone_seed %in% c(1, 2)) %>%
    st_transform(crs = st_crs(3112)) %>% 
    ggplot() +
    geom_sf(aes(fill = raw_SIR))+
    coord_sf(xlim = lims$Brisbane$xlim,
             ylim = lims$Brisbane$ylim)+
    #theme_void()+
    scale_fill_gradient2(low = "blue", mid = "white",
                         high = "red", midpoint = 1) + 
    theme(
      legend.position = "bottom",
      plot.title = element_text(margin = margin(0,0,2,0)),
      plot.margin = unit(c(1,1,1,1), "mm")
    ) +
    guides(fill = guide_colourbar(barwidth = 20)
    ) +
    labs(fill = "",
         title = "Raw SIR for Brisbane",
         subtitle = paste0("gamma: ", pos_grid$gamma[i],
                        ", alpha: ",  pos_grid$alpha[i]))+
    facet_grid(zone_seed ~ t, labeller = purrr::partial(label_both, sep = " = "))
  
  ggsave(file = paste0(plot_loc, "/brisbane_", plot_name, ".png"),
         width = 7.3, height = 6.76)
  
  message(paste0("Finished loop ", i, " out of ", nrow(pos_grid), "."))
  
}



################################################################################
## DEPRECIATED ## --------------------------------------------------------------

## Boxplot: zoning distributions for leroux model ----
models %>%
  filter(term == "x") %>% 
  ggplot(aes(x = lerouxcb_median, fill = t, group = t))+
  theme_bw()+
  geom_boxplot() +
  facet_grid(alpha~gamma, labeller = purrr::partial(label_both, sep = " = "))+
  labs(fill = "Target Population",
       title = "Zoning distributions for covariate effect",
       subtitle = "Posterior median from Leroux model",
       x = "Posterior median for covariate")
ggsave(file = paste0(plot_loc, "/boxplot_leroux_x.png"),
       width = 11.3, height = 6.76)

## Boxplot: zoning distributions for poisson model ----
models %>%
  filter(term == "x") %>% 
  ggplot(aes(x = nonspatial_estimate, fill = t, group = t))+
  theme_bw()+
  geom_boxplot() +
  facet_grid(alpha~gamma, labeller = purrr::partial(label_both, sep = " = "))+
  labs(fill = "Target Population",
       title = "Zoning distributions for covariate effect",
       subtitle = "MLE estimate from Poisson model",
       x = "MLE estimate")
ggsave(file = paste0(plot_loc, "/boxplot_poisson_x.png"),
       width = 11.3, height = 6.76)