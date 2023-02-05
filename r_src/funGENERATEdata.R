
setupSIMS <- function(rho, y_av, dir){

dir <- paste0(dir, "/")
sa1 <- readRDS(paste0(dir, "sa1_shape/sa1data.rds"))
L <- readRDS(paste0(dir, "phi/L2_rho", rho, ".rds"))

## Create the weight matrix ## -------------------------------------------------

nb <- poly2nb(sa1, queen = F)
listw <- nb2listw(nb, style="B", zero.policy = TRUE)

## Append to out_list ## -------------------------------------------------------
out_list <- list()
out_list$sa1 <- sa1
out_list$rho <- rho
out_list$y_av <- y_av
out_list$nb <- nb
out_list$listw <- listw
out_list$L <- L
return(out_list)

}

simDATA <- function(input_list){
  
sa1 <- input_list$sa1
L <- input_list$L
N <- nrow(input_list$sa1)
y_av <- input_list$y_av

## generate spatial random effect ## -------------------------------------------

Z <- rnorm(nrow(sa1))
phi <- backsolve(L, Z)

## Form the data ## ------------------------------------------------------------

# add the spatial random effects to the dataset
sa1$phi <- phi

# generate continuous covariate
sa1$x <- runif(N, -3, 3)

# create the linear predictor
beta <- 1
phi_sd <- 2
eta = log(sa1$pop) + y_av + (beta * sa1$x) + sa1$phi
mu <- exp(eta)

# draw from the poisson distribution
y_initial <- rpois(N, mu)

# add to dataset
sa1$mu <- mu
sa1$y_initial <- y_initial

# Rescale y so it matches desired level of counts
sa1 <- sa1 %>% 
  mutate(y = round(y_initial / sum(y_initial) * y_av * N))

# get expected counts
# get the overall rate using the simulated y's
grand_rate <- sum(sa1$y)/sum(sa1$pop)
E_initial <- sa1$pop*grand_rate

# Ensure that sum(y) = sum(E)
sa1$E <- (E_initial / sum(E_initial)) * sum(sa1$y)

# get the raw SIR
sa1 <- sa1 %>% 
  mutate(raw_SIR = y/E)

## output data ## -----------------------------------------------------------------

return(sa1)

}
