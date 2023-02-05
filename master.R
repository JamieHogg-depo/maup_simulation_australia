##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                                                            ~~
##                   EXTRACT COEFFICIENTS AFTER AGGREGATING                 ----
##                                                                            ~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Load packages ## ------------------------------------------------------------
library(tidyverse)
library(spdep)
library(sf)
library(Matrix)
library(MASS)
library(tidymodels)
library(CARBayes)
library(spatialreg)
library(rstan)
library(lme4)
library(broom.mixed)
library(AER)

location_src <- paste0(getwd(), "/CDS_MAUP2/r_src")

## ---- Load functions ---- ##
source(paste0(location_src, "/funGENERATEdata.R"))
source(paste0(location_src, "/functions.R"))

## Create grid of rho and y_av params
 param_grid <- expand.grid(rho = c(0.05, 0.5, 0.95),
                           y_av = c(5, 25, 50))

# number of repeated datasets
nu_datasets <- 1 # DO NOT TOUCH

# target populations
target_pop <- c(5000, 10000, 50000, 90000, 200000, 300000, 700000)
# target_pop <- c(5000, 10000)

# number of different seeds
# nu_zone_seeds <- 2

# create loop counter
c <- 0
c_max <- c(length(target_pop),
           nu_zone_seeds) %>% 
		   reduce(`*`)

# Collect all zonation file names
zones_filenames <- list.files(paste0(location_src, "/zones"))
zones_pop_targets <- as.numeric(regmatches(zones_filenames, regexpr("([0-9]+)", zones_filenames)))

## ---- Create list to store results ---- ##

sim_list <- list(name = list(),
                 target_pop = list(),
                 model = list(),
                 MI = list(),
				 OD = list(),
                 data_zone_metrics = list(),
                 df_agg = list())

stored_data_list <- list(data = list(), 
                         MI = list(), 
                         sa1model = list(),
                         sa2model = list(),
                         sa3model = list(),
                         sa4model = list())

## ---- LOOP: rho and y_av ---- ##
#for(QaS in 1:nrow(param_grid)){

# NOTE: setupSIMS is not stochastic
set <- setupSIMS(rho = param_grid$rho[[QaS]], 
                 y_av = param_grid$y_av[[QaS]], 
                 dir = location_src)

## ---- LOOP: new datasets under same scenario

for(BaS in 1:nu_datasets){

# NOTE: Each run of simDATA() will generate a different data set  
sa1 <- simDATA(set)

# update list with dataset
stored_data_list$data[[QaS]] <- sa1 %>%
  mutate(rho = param_grid$rho[[QaS]], 
         y_av = param_grid$y_av[[QaS]])

# update list with moran's I
nb <- poly2nb(sa1, queen = F)
listw <- nb2listw(nb, style="B", zero.policy = TRUE)
stored_data_list$MI[[QaS]] <- rbind(
  tidyMI(moran.test(sa1$y, listw = listw, zero.policy = T), "y"),
  tidyMI(moran.test(sa1$raw_SIR, listw = listw, zero.policy = T), "raw_SIR")) %>% 
  mutate(rho = param_grid$rho[[QaS]], 
         y_av = param_grid$y_av[[QaS]])

# fit the leroux model to the sa1 data
		df_agg <- sa1 %>% 
		rename(sum_y = y, sum_E = E) %>%
		mutate(log_sum_E = log(sum_E),
			   zoneID = row_number())
		
if(fit_sa1_data == TRUE){
		 
		W <- nb2mat(nb, style = "B")

		# leroux
		start_time_model <- Sys.time()
		tsil <- collect_3chains_fromCARBayes(get_metrics = FALSE)
		duration_model <- Sys.time() - start_time_model
		duration_model <- round(as.numeric(duration_model, units = "mins"), 1)
		message(paste0("CARBayes fit to SA1 data took ", duration_model, " minutes."))
		
		# poisson
		fit_nonspatial <- glm(sum_y ~ offset(log(sum_E)) + x + log_sum_E, data = df_agg, family = poisson)
		R_nonspatial <- tidy(fit_nonspatial) %>% 
		  setNames(c("term", paste0("nonspatial_", names(.)[-1])))
		  
		# Generic poisson model: with generic random effect
		poisson_rand <- glmer(sum_y ~ offset(log(sum_E)) + x + log_sum_E + (1 | zoneID), data = df_agg,
					  family = "poisson")
		R_rand <- tidy(poisson_rand) %>% 
		dplyr::select(-c(effect, group)) %>% 
		setNames(c("term", paste0("rand_", names(.)[-1])))
		
		# append model results to output list
		stored_data_list$sa1model[[QaS]] <- list(R_nonspatial, tsil$results, R_rand) %>% 
		  reduce(full_join, by = "term") %>% 
		  mutate(rho = param_grid$rho[[QaS]], 
		         y_av = param_grid$y_av[[QaS]],
				 poisson_normAIC = as.numeric((-2*logLik(fit_nonspatial) + 2*2)/(nrow(df_agg))))
		
		# remove objects
		rm(nb, listw, W, df_agg, 
		   start_time_model, duration_model,
		   tsil, fit_nonspatial,
		   R_nonspatial)
		
} else{		
  
		# poisson
		fit_nonspatial <- glm(sum_y ~ offset(log(sum_E)) + x + log_sum_E, data = df_agg, family = poisson)
		R_nonspatial <- tidy(fit_nonspatial) %>% 
		  setNames(c("term", paste0("nonspatial_", names(.)[-1])))
		  
		# Generic poisson model: with generic random effect
		poisson_rand <- glmer(sum_y ~ offset(log(sum_E)) + x + log_sum_E + (1 | zoneID), data = df_agg,
					  family = "poisson")
		R_rand <- tidy(poisson_rand) %>% 
		dplyr::select(-c(effect, group)) %>% 
		setNames(c("term", paste0("rand_", names(.)[-1])))
		
		# append model results to output list
		stored_data_list$sa1model[[QaS]] <- list(R_nonspatial, R_rand) %>% 
		  reduce(full_join, by = "term") %>% 
		  mutate(rho = param_grid$rho[[QaS]], 
		         y_av = param_grid$y_av[[QaS]],
				 poisson_normAIC = as.numeric((-2*logLik(fit_nonspatial) + 2*2)/(nrow(df_agg))))

		# remove objects
		rm(df_agg, listw, fit_nonspatial)
		
}
		
# fit the leroux model to the sa2 data
		df_agg <- sa1 %>% 
		  mutate(sa2 = str_sub(SA1_MAIN16, 1, 9)) %>% 
		  group_by(sa2) %>% 
		  summarise(x = mean(x),
					sum_y = sum(y),
					sum_E = sum(E), 
					log_sum_E = log(sum_E),
					geometry = st_union(geometry)) %>%
		  mutate(zoneID = row_number())
		
		nb <- poly2nb(df_agg, queen = F)
		W <- nb2mat(nb, style = "B")

		# leroux
		start_time_model <- Sys.time()
		tsil <- collect_3chains_fromCARBayes(get_metrics = FALSE)
		duration_model <- Sys.time() - start_time_model
		duration_model <- round(as.numeric(duration_model, units = "mins"), 1)
		message(paste0("CARBayes fit to SA2 data took ", duration_model, " minutes."))
		
		# poisson
		fit_nonspatial <- glm(sum_y ~ offset(log(sum_E)) + x + log_sum_E, data = df_agg, family = poisson)
		R_nonspatial <- tidy(fit_nonspatial) %>% 
		  setNames(c("term", paste0("nonspatial_", names(.)[-1])))
		  
		# Generic poisson model: with generic random effect
		poisson_rand <- glmer(sum_y ~ offset(log(sum_E)) + x + log_sum_E + (1 | zoneID), data = df_agg,
					  family = "poisson")
		R_rand <- tidy(poisson_rand) %>% 
		dplyr::select(-c(effect, group)) %>% 
		setNames(c("term", paste0("rand_", names(.)[-1])))
		
		# append model results to output list
		stored_data_list$sa2model[[QaS]] <- list(R_nonspatial, tsil$results, R_rand) %>% 
		  reduce(full_join, by = "term") %>% 
		  mutate(rho = param_grid$rho[[QaS]], 
		         y_av = param_grid$y_av[[QaS]],
				 poisson_normAIC = as.numeric((-2*logLik(fit_nonspatial) + 2*2)/(nrow(df_agg))))
		
		# remove objects
		rm(nb, W, df_agg, 
		   start_time_model, duration_model,
		   tsil, fit_nonspatial,
		   R_nonspatial)
		
# fit the leroux model to the sa3 data
		df_agg <- sa1 %>% 
		  mutate(sa3 = str_sub(SA1_MAIN16, 1, 5)) %>% 
		  group_by(sa3) %>% 
		  summarise(x = mean(x),
					sum_y = sum(y),
					sum_E = sum(E), 
					log_sum_E = log(sum_E),
					geometry = st_union(geometry)) %>%
		  mutate(zoneID = row_number())
		
		nb <- poly2nb(df_agg, queen = F)
		W <- nb2mat(nb, style = "B")

		# leroux
		start_time_model <- Sys.time()
		tsil <- collect_3chains_fromCARBayes(get_metrics = FALSE)
		duration_model <- Sys.time() - start_time_model
		duration_model <- round(as.numeric(duration_model, units = "mins"), 1)
		message(paste0("CARBayes fit to SA3 data took ", duration_model, " minutes."))
		
		# poisson
		fit_nonspatial <- glm(sum_y ~ offset(log(sum_E)) + x + log_sum_E, data = df_agg, family = poisson)
		R_nonspatial <- tidy(fit_nonspatial) %>% 
		  setNames(c("term", paste0("nonspatial_", names(.)[-1])))
		  
		# Generic poisson model: with generic random effect
		poisson_rand <- glmer(sum_y ~ offset(log(sum_E)) + x + log_sum_E + (1 | zoneID), data = df_agg,
					  family = "poisson")
		R_rand <- tidy(poisson_rand) %>% 
		dplyr::select(-c(effect, group)) %>% 
		setNames(c("term", paste0("rand_", names(.)[-1])))
		
		# append model results to output list
		stored_data_list$sa3model[[QaS]] <- list(R_nonspatial, tsil$results, R_rand) %>% 
		  reduce(full_join, by = "term") %>% 
		  mutate(rho = param_grid$rho[[QaS]], 
		         y_av = param_grid$y_av[[QaS]],
				 poisson_normAIC = as.numeric((-2*logLik(fit_nonspatial) + 2*2)/(nrow(df_agg))))
		
		# remove objects
		rm(nb, W, df_agg, 
		   start_time_model, duration_model,
		   tsil, fit_nonspatial,
		   R_nonspatial)
		
# fit the leroux model to the sa4 data
		df_agg <- sa1 %>% 
		  mutate(sa4 = str_sub(SA1_MAIN16, 1, 3)) %>% 
		  group_by(sa4) %>% 
		  summarise(x = mean(x),
					sum_y = sum(y),
					sum_E = sum(E), 
					log_sum_E = log(sum_E),
					geometry = st_union(geometry)) %>%
		  mutate(zoneID = row_number())
		
		nb <- poly2nb(df_agg, queen = F)
		W <- nb2mat(nb, style = "B")

		# leroux
		start_time_model <- Sys.time()
		tsil <- collect_3chains_fromCARBayes(get_metrics = FALSE)
		duration_model <- Sys.time() - start_time_model
		duration_model <- round(as.numeric(duration_model, units = "mins"), 1)
		message(paste0("CARBayes fit to SA4 data took ", duration_model, " minutes."))
		
		# poisson
		fit_nonspatial <- glm(sum_y ~ offset(log(sum_E)) + x + log_sum_E, data = df_agg, family = poisson)
		R_nonspatial <- tidy(fit_nonspatial) %>% 
		  setNames(c("term", paste0("nonspatial_", names(.)[-1])))
		  
		# Generic poisson model: with generic random effect
		poisson_rand <- glmer(sum_y ~ offset(log(sum_E)) + x + log_sum_E + (1 | zoneID), data = df_agg,
					  family = "poisson")
		R_rand <- tidy(poisson_rand) %>% 
		dplyr::select(-c(effect, group)) %>% 
		setNames(c("term", paste0("rand_", names(.)[-1])))
		
		# append model results to output list
		stored_data_list$sa4model[[QaS]] <- list(R_nonspatial, tsil$results, R_rand) %>% 
		  reduce(full_join, by = "term") %>% 
		  mutate(rho = param_grid$rho[[QaS]], 
		         y_av = param_grid$y_av[[QaS]],
				 poisson_normAIC = as.numeric((-2*logLik(fit_nonspatial) + 2*2)/(nrow(df_agg))))
		
		# remove objects
		rm(nb, W, df_agg, 
		   start_time_model, duration_model,
		   tsil, fit_nonspatial,
		   R_nonspatial)

## ---- LOOP: Different zonations ---- ##

start <- Sys.time()

for(LaS in 1:length(target_pop)){
  
cur_target_pop <- target_pop[LaS]
cur_file_index <- which(zones_pop_targets == target_pop[LaS])
  
## ---- LOOP: Different seeds for same zonation ---- ##

for(JaS in 1:nu_zone_seeds){
  
c <- c + 1  # update counter

cur_file_zonation <- zones_filenames[cur_file_index[JaS]]

# extract the target population from file name
sim_list$target_pop[[c]] <- as.character(cur_target_pop)
sim_list$name[[c]] <- str_remove(cur_file_zonation, ".csv")
# assign a temp data
temp_data <- sa1

# load zones
df_zones <- suppressMessages(read_csv(paste0(location_src, "/zones/", cur_file_zonation)) %>% 
  rename(id = BldBlID) %>% 
  mutate(id = id + 1))

# relabel the zones
cur_zone_labs <- data.frame(TractID = unique((df_zones %>% arrange(id))$TractID),
                            zoneID = 1:length(unique((df_zones %>% arrange(id))$TractID)))
df_zones <- df_zones %>% 
  left_join(.,cur_zone_labs, by = "TractID") %>% 
  dplyr::select(-TractID)

# join to dataset
temp_data <- temp_data %>%  
  inner_join(.,df_zones, by = "id")

# aggregate
df_agg <- temp_data %>% 
  group_by(zoneID) %>% 
  mutate(mean_y = mean(y, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(abs_error = abs(y - mean_y)) %>% 
  group_by(zoneID) %>% 
  summarise(Nu_sa1s = n(),
            sum_pop = sum(pop),
            mean_pop = mean(pop),
            var_pop = var(pop, na.rm = T),
            x = mean(x),
            sum_y = sum(y),
            mean_y = mean(y, na.rm = T),
            var_y = var(y, na.rm = T),
            mean_abs_error = mean(abs_error, na.rm = T), 
            sum_E = sum(E), 
			log_sum_E = log(sum_E),
            geometry = st_union(geometry)) %>% 
  mutate(raw_SIR = sum_y/sum_E,
         zone_seed = JaS,
         target_pop = cur_target_pop, 
         rho = param_grid$rho[[QaS]], 
         y_av = param_grid$y_av[[QaS]])

# add df_agg to sim_list
sim_list$df_agg[[c]] <- df_agg

# get the spatial objects from our collapsed dataset
temp_nb <- poly2nb(df_agg, queen = F)
temp_listw <- nb2listw(temp_nb, style="B", zero.policy = TRUE)
W <- nb2mat(temp_nb, style = "B")

# Fit the models
# generic poisson model: NOT ACCOUNTING FOR SPATIAL AUTOCORRELATION
fit_nonspatial <- glm(sum_y ~ offset(log(sum_E)) + x + log_sum_E, data = df_agg, family = poisson)
R_nonspatial <- tidy(fit_nonspatial) %>% 
  setNames(c("term", paste0("nonspatial_", names(.)[-1])))
sim_list$OD[[c]] <- tidy(dispersiontest(fit_nonspatial))

# Generic poisson model: with generic random effect
poisson_rand <- glmer(sum_y ~ offset(log(sum_E)) + x + log_sum_E + (1 | zoneID), data = df_agg,
                      family = "poisson")
R_rand <- tidy(poisson_rand) %>% 
  dplyr::select(-c(effect, group)) %>% 
  setNames(c("term", paste0("rand_", names(.)[-1])))

# CARBAYES: Leroux
tsil <- collect_3chains_fromCARBayes()
R_leroux_cb <- tsil$results

# Add some summary statistics for the aggregation
sim_list$data_zone_metrics[[c]] <- data.frame(median_sa1s_zones = median(df_agg$Nu_sa1s),
											  poisson_normAIC = as.numeric((-2*logLik(fit_nonspatial) + 2*2)/(nrow(df_agg))),
											  leroux_normAIC_pw = tsil$leroux_normAIC_pw,
											  leroux_normAIC_pd = tsil$leroux_normAIC_pd,
											  cor_leroux_rand = cor(tsil$fitted.values, predict(poisson_rand, type = "response")),
											  cor_leroux_nonspatial = cor(tsil$fitted.values, predict(fit_nonspatial, type = "response")),
                                              Nbar = median(df_agg$sum_pop),
                                              N = nrow(df_agg),
                                              target_pop = cur_target_pop, 
                                              zone_seed = JaS,
                                              rho = param_grid$rho[[QaS]], 
                                              y_av = param_grid$y_av[[QaS]], 
                                              data_seed = BaS)
											  
# append model results to output list
sim_list$model[[c]] <- list(R_nonspatial, 
							R_leroux_cb,
							R_rand) %>% 
  reduce(full_join, by = "term") %>% 
  mutate(target_pop = cur_target_pop, 
         zone_seed = JaS,
         rho = param_grid$rho[[QaS]], 
         y_av = param_grid$y_av[[QaS]], 
         data_seed = BaS)

## Moran's I ## ----------------------------------------------------------------
sim_list$MI[[c]] <- rbind(
tidyMI(moran.test(df_agg$sum_y, listw = temp_listw, zero.policy = T), "y"),
tidyMI(moran.test(df_agg$raw_SIR, listw = temp_listw, zero.policy = T), "raw_SIR"),
tidyMI(moran.test(df_agg$sum_E, listw = temp_listw, zero.policy = T), "E"),
tidyMI(moran.test(df_agg$sum_pop, listw = temp_listw, zero.policy = T), "pop")
) %>% 
  mutate(target_pop = cur_target_pop, 
         zone_seed = JaS,
         rho = param_grid$rho[[QaS]], 
         y_av = param_grid$y_av[[QaS]], 
         data_seed = BaS)

# drop some objects
rm(temp_nb, temp_listw, df_agg, temp_data, df_zones, cur_zone_labs, W, 
   fit_nonspatial, R_nonspatial, R_leroux_cb)

# Print some helpful messages to the console
  message(paste0("Loop: ", c, "(", c_max, "), ",
                 "Rho: ", param_grid$rho[[QaS]], ", ",
                 "y_av: ", param_grid$y_av[[QaS]], ", ",
                 "Data_seed: ", BaS, "(", nu_datasets, "), ",
                 "Target pop: ", cur_target_pop, ", ",
                 "Zone_seed: ", JaS, "(", nu_zone_seeds, ")" ))
  
} # for JaS - number of zone seeds

} # for LaS - different zonations

} # for BaS - number of datasets

#} # for QaS - rho and y_av

end <- Sys.time()
duration <- end - start
duration <- round(as.numeric(duration, units = "mins"), 1)
message(paste0("With ", c_max, " loops, this simulation took ", duration, 
               " minutes. The average time per loop is ", round(duration/c_max, 1), " minutes."))

# Save output
saveRDS(sim_list, file = paste0(getwd(), "/CDS_MAUP2/outputs/", cur_date, "/r/Qas", QaS, "_simlist.rds"))
saveRDS(stored_data_list, file = paste0(getwd(), "/CDS_MAUP2/outputs/", cur_date, "/r/Qas", QaS, "_storeddatalist.rds"))

## ---- Results ---- ##

# average time per loop is 0.3 minutes
# With 96 loops, this simulation took 30.9 minutes. The average time per loop is 0.3 minutes.
