## ----setup, include=FALSE-----------------------------------------------------
library(EpiPvr)
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(
  comment = "#>",
  tidy = FALSE,
  width = 70
)

## ----include=FALSE------------------------------------------------------------
message("Now running: Chunk 1 setups...")
flush.console()
#options(mc.cores = parallel::detectCores()) 
#Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')
#if (requireNamespace("rstan", quietly = TRUE)) {
#    rstan::rstan_options(auto_write = TRUE)
#    options(mc.cores = 1)
#}

## ----presets, cache=FALSE-----------------------------------------------------

virus_type='SPT' 
lsEst_in=40; # upper bound on per insect survival in the experiment in days


## ----dataset, cache=FALSE-----------------------------------------------------

data("ap_data_sim_SPT", package = "EpiPvr")
print(ap_data_sim_SPT)


## ----sample, cache=FALSE, message=FALSE, warning=FALSE------------------------
if (identical(Sys.getenv("NOT_CRAN"), "true")) {
EVSPT_sim=estimate_virus_parameters_SPT(ap_data_sim_SPT,lsEst_in,D_numPtsPdin=1,
                                        mcmcOptions=c(4500,6000),numChainsIn=4,mc.parallel=0)
} else {
  # Load precomputed object for CRAN
  EVSPT_sim <- readRDS(system.file("extdata", "EVSPT_sim.rds", package = "EpiPvr"))
}

target=EVSPT_sim$array1 # fixing the first output array, contains the sample Markov chains


## ----print, cache=FALSE-------------------------------------------------------

# viewing output content
print(dim(EVSPT_sim$array0))
print(head(EVSPT_sim$array0))
print(dim(EVSPT_sim$array1))
print(head(EVSPT_sim$array1))
print(dim(EVSPT_sim$array2))
print(head(EVSPT_sim$array2))
print(length(EVSPT_sim$array3))
print(head(EVSPT_sim$array3))
print(length(EVSPT_sim$array4))
print(head(EVSPT_sim$array4))
print(length(EVSPT_sim$array5))
print(head(EVSPT_sim$array5))
print(length(EVSPT_sim$converge_results))
print(head(EVSPT_sim$converge_results))

runs=dim(target)[1]
nChains=dim(target)[2]


## ----SPTprintw, cache=FALSE---------------------------------------------------
if (length(warnings()) > 0) {
  message("Warnings occurred during fitting; check convergence diagnostics.")
}

## ----plotSPT, fig.width=5, fig.height=5, echo=FALSE---------------------------

#bayesplot::mcmc_pairs(EVSPT_sim$array0,diag_fun = "dens")



## ----plotSPT2, fig.width=5, fig.height=5, echo=FALSE--------------------------

#bayesplot::mcmc_dens(EVSPT_sim$array1,facet_args = list(ncol = 3))+ ggplot2::theme(aspect.ratio = 0.5)


## ----plotSPT3, echo=FALSE-----------------------------------------------------

plot(EVSPT_sim$array3$lenA,EVSPT_sim$array3$propA, xlab = "AAP duration, hours", ylab = "Prop. test plants SPT infected",ylim=c(0,1))
lines(EVSPT_sim$array3$lenA,EVSPT_sim$array3$simulLA)
lines(EVSPT_sim$array3$lenA,EVSPT_sim$array3$simulUA)
legend(EVSPT_sim$array3$lenA[1], 0.5, legend=c("experiment", "95% Cr.I. estimates"), pch = c(1, NA), lty = c(NA, 1), col = "black", cex=0.8, bty = "n")
plot(EVSPT_sim$array4$lenI,EVSPT_sim$array4$propI, xlab = "IAP duration, hours", ylab = "Prop. test plants SPT #infected",ylim=c(0,1))
lines(EVSPT_sim$array4$lenI,EVSPT_sim$array4$simulLI)
lines(EVSPT_sim$array4$lenI,EVSPT_sim$array4$simulUI)
legend(EVSPT_sim$array4$lenI[1], 0.5, legend=c("experiment", "95% Cr.I. estimates"), pch = c(1, NA), lty = c(NA, 1), col = "black", cex=0.8, bty = "n")



## ----unpack, cache=FALSE------------------------------------------------------

# P EPIDEMIC inferences
# accessing the chains from estimate_virus_parameters
# VIRUS PARAMETERS

# convert from per hours to per day
al_fits=target[,,dimnames(target)$variable=="al[1]", drop = TRUE]*24
be_fits=target[,,dimnames(target)$variable=="be[1]", drop = TRUE]*24
mu_fits=target[,,dimnames(target)$variable=="mu[1]", drop = TRUE]*24

#confirm there are no negative values
print(sum((al_fits<0)))
print(sum((be_fits<0)))
print(sum((mu_fits<0)))

# LOCAL PARAMETERS
# set the local parameters rates (per day)
thet_external <- 0.45 # dispersal
r_external  <- 1/28 # roguing
bf_external <- 1/14  # vector mortality
h_external <- 1/365  # harvesting rate
nu_pl_external <- 1/14  # plant latent progression rate (1/plant latent period)
localParams=c(thet_external, r_external, h_external, bf_external, nu_pl_external)
# generate p Epdemic results for given insect burden (per plant)
numInsects=3

## ----infer, cache=FALSE-------------------------------------------------------

runsPrag=10 

# initial guess for expected time in hours until next event - the function will use this 
# to work out efficient time step (optional, default = 0.1)
interval_ind=1/10 #

# initialising vectors for storing chain samples
result_vec_byPl <- array(rep(0, runsPrag*nChains), dim=c(runsPrag,nChains))
result_vec_byIns <- result_vec_byPl


if (identical(Sys.getenv("NOT_CRAN"), "true")) {

# loop through MCMC chains to sample posterior distributions of virus parameters and 
# calculate epidemic probability for a given vector abundance
for (ggg in 1:nChains){
  print(paste0('nChains',ggg))
  for (ppp in 1:runsPrag) {
    al_estim=al_fits[nrow(al_fits)+1-ppp,ggg]
    be_estim=be_fits[nrow(al_fits)+1-ppp,ggg]
    mu_estim=mu_fits[nrow(al_fits)+1-ppp,ggg]
    virusParams=c(al_estim,be_estim,mu_estim)
    # returns vector of epidemic probabilities for different inoculum states (see ms)
    numVars <- ((numInsects+1)*3)-1  # insects per plant number no. inoculum states
    qm_out=calculate_epidemic_probability(numInsects,interval_ind,localParams,virusParams)
    result_vec_byPl[ppp,ggg]=qm_out[1]  # storing epi prob for state: inoculum=single infected plant
    result_vec_byIns[ppp,ggg]=qm_out[(numVars-(numInsects-1))]  # storing for state: single infected insect
    
  }
}

} else {
  # Load precomputed object for CRAN
  result_vec_byPl <- readRDS(system.file("extdata", "result_vec_byPl.rds", package = "EpiPvr"))
  result_vec_byIns <- readRDS(system.file("extdata", "result_vec_byIns.rds", package = "EpiPvr"))
}



## ----distributions, cache=FALSE-----------------------------------------------

# data wrangling the set of epidemic probabilities from the chains
data_table_Pl=array(rep(0,3), dim=c(3, 1))
data_table_Ins=array(rep(0,3), dim=c(3, 1))
target_A1A2=array(rep(0,2*runsPrag*nChains), dim=c(2,runsPrag*nChains))

# for a given level of insect use 'summarise_draws' to calculate credible intervals
target_A1=array(rep(0,runsPrag*nChains), dim=c(runsPrag, nChains))
target_A2=target_A1
for (gg in 1:nChains){
  target_A1[,gg]=t(result_vec_byPl[,gg])
  target_A2[,gg]=t(result_vec_byIns[,gg])
}
target_A1A2=rbind(as.vector(target_A1),as.vector(target_A2))

# first the first few elements of epidemic probability from plant inocula
print(head(target_A1A2))
# and from insect inocula
print(head(target_A1A2))


## ----tables1, cache=FALSE-----------------------------------------------------


# ensure at least ~20 data points for each of 4 chains (ie cant summarise distribution if v small)
ts=posterior::summarize_draws(posterior::as_draws(t(target_A1A2))) 
data_table_Pl=c(ts[1,]$median,ts[1,]$q5,ts[1,]$q95)
data_table_Ins=c(ts[2,]$median,ts[2,]$q5,ts[2,]$q95)

## ----tables2, cache=FALSE-----------------------------------------------------
names(data_table_Pl)=c('50','5','95')
print(data_table_Pl)

names(data_table_Ins)=c('50','5','95')
print(data_table_Ins)


## ----PTpreandfit, cache=FALSE-------------------------------------------------

virus_type='PT'
data("ap_data_sim_PT", package = "EpiPvr")
print(ap_data_sim_PT)


## ----PTsample, cache=TRUE-----------------------------------------------------
if (identical(Sys.getenv("NOT_CRAN"), "true")) {
EVPT_sim=estimate_virus_parameters_PT(ap_data_sim_PT,lsEst_in,D_numPtsPdin=1,mcmcOptions=c(1000,2000),numChainsIn=4,mc.parallel=0)
} else {
  # Load precomputed object for CRAN
  EVPT_sim <- readRDS(system.file("extdata", "EVPT_sim.rds", package = "EpiPvr"))
}

target=EVPT_sim$array1

runs=dim(target)[1]
nChains=dim(target)[2]


## ----PTprint, cache=FALSE-----------------------------------------------------

print(dim(EVPT_sim$array0))
head(EVPT_sim$array0)
print(dim(EVPT_sim$array1))
head(EVPT_sim$array1)
print(dim(EVPT_sim$array2))
head(EVPT_sim$array2)
print(length(EVPT_sim$array3))
head(EVPT_sim$array3)
print(length(EVPT_sim$array4))
head(EVPT_sim$array4)
print(length(EVPT_sim$array5))
head(EVPT_sim$array5)
print(length(EVPT_sim$array6))
head(EVPT_sim$array6)
print(length(EVPT_sim$converge_results))
head(EVPT_sim$converge_results)


## ----PTprintw, cache=FALSE----------------------------------------------------
if (length(warnings()) > 0) {
  message("Warnings occurred during fitting; check convergence diagnostics.")
}

## ----PTplot, fig.width=5, fig.height=5, echo=FALSE----------------------------

# bayesplot::mcmc_pairs(EVPT_sim$array0,diag_fun = "dens")


## ----PTplot1, fig.width=5, fig.height=5, echo=FALSE---------------------------

#bayesplot::mcmc_dens(EVPT_sim$array1,facet_args = list(ncol = 4))+ ggplot2::theme(aspect.ratio = 0.5)


## ----PTplot2, echo=FALSE------------------------------------------------------

plot(EVPT_sim$array3$lenA,EVPT_sim$array3$propA, xlab = "AAP duration, hours", ylab = "Prop. test plants PT infected",ylim=c(0,1))
lines(EVPT_sim$array3$lenA,EVPT_sim$array3$simulLA)
lines(EVPT_sim$array3$lenA,EVPT_sim$array3$simulUA)
legend(EVPT_sim$array3$lenA[1], 0.5, legend=c("experiment", "95% Cr.I. estimates"), pch = c(1, NA), lty = c(NA, 1), col = "black", cex=0.8, bty = "n")
plot(EVPT_sim$array4$lenL,EVPT_sim$array4$propL, xlab = "LAP duration, hours", ylab = "Prop. test plants PT #infected",ylim=c(0,1))
lines(EVPT_sim$array4$lenL,EVPT_sim$array4$simulLL)
lines(EVPT_sim$array4$lenL,EVPT_sim$array4$simulUL)
legend(EVPT_sim$array4$lenL[1], 0.5, legend=c("experiment", "95% Cr.I. estimates"), pch = c(1, NA), lty = c(NA, 1), col = "black", cex=0.8, bty = "n")
plot(EVPT_sim$array5$lenI,EVPT_sim$array5$propI, xlab = "IAP duration, hours", ylab = "Prop. test plants PT #infected",ylim=c(0,1))
lines(EVPT_sim$array5$lenI,EVPT_sim$array5$simulLI)
lines(EVPT_sim$array5$lenI,EVPT_sim$array5$simulUI)
legend(EVPT_sim$array5$lenI[1], 0.5, legend=c("experiment", "95% Cr.I. estimates"), pch = c(1, NA), lty = c(NA, 1), col = "black", cex=0.8, bty = "n")


## ----PTinfer, cache=FALSE-----------------------------------------------------
# P EPIDEMIC inferences
# accessing the chains from estimate_virus_parameters
# VIRUS PARAMETERS

# convert from per min to per day
al_fits=target[,,dimnames(target)$variable=="al[1]", drop = TRUE]*24
be_fits=target[,,dimnames(target)$variable=="be[1]", drop = TRUE]*24
mu_fits=target[,,dimnames(target)$variable=="mu[1]", drop = TRUE]*24
print(sum((al_fits<0)))
print(sum((be_fits<0)))
print(sum((mu_fits<0)))

# LOCAL PARAMETERS
# set the local parameters (per day)
thet_external <- 0.45 # dispersal
r_external <- 1/28 # roguing
bf_external <- 1/14 # vector mortality
h_external <- 1/365 # harvesting rate
nu_pl_external <- 1/14 # plant latent progression rate (1/plant latent period)
localParams=c(thet_external, r_external, h_external, bf_external, nu_pl_external)

# generate p Epdemic results for various insect burdens (per plant)
interval_ind=(1/10) # this is an optional start setting relating to efficiency of time-step
result_vec_byPl_PT <- array(rep(0, runsPrag*nChains), dim=c(runsPrag,nChains))
result_vec_byIns_PT <- result_vec_byPl

if (identical(Sys.getenv("NOT_CRAN"), "true")) {

# loop through MCMC chains to sample posterior distributions of virus parameters
# and calc epi prob for given vector abundance
for (ggg in 1:nChains){
 print(paste0('nChains',ggg))
 for (ppp in 1:runsPrag) {

  al_estim=al_fits[nrow(al_fits)+1-ppp,ggg]
  be_estim=be_fits[nrow(al_fits)+1-ppp,ggg]
  mu_estim=mu_fits[nrow(al_fits)+1-ppp,ggg]
  virusParams=c(al_estim,be_estim,mu_estim)
  
  # returns vector of epidemic probabilities for different inoculum states (see ms)
  numVars <- ((numInsects+1)*3)-1 # insects per plant number no. inoculum states
  qm_out=calculate_epidemic_probability(numInsects,interval_ind,localParams,virusParams)
  
  result_vec_byPl_PT[ppp,ggg]=qm_out[1] # storing epi prob for state: inoculum=single infected plant
  result_vec_byIns_PT[ppp,ggg]=qm_out[(numVars-(numInsects-1))] # storingfor state: single infected insect
  
 }
}

} else {
  # Load precomputed object for CRAN
  result_vec_byPl_PT <- readRDS(system.file("extdata", "result_vec_byPl_PT.rds", package = "EpiPvr"))
  result_vec_byIns_PT <- readRDS(system.file("extdata", "result_vec_byIns_PT.rds", package = "EpiPvr"))
}




## ----wrangleandout, cache=FALSE-----------------------------------------------

# now run the epidemic probability calculator
# runtime for a single parameter set is fast
# runtime for a MCMC chains (potentially very many parameter sets) may take substantial time

# data wrangling the set of epidemic probabilities from the chains
data_table_Pl=array(rep(0,3), dim=c(3, 1))
data_table_Ins=array(rep(0,3), dim=c(3, 1))
target_A1A2=array(rep(0,2*runsPrag*nChains), dim=c(2,runsPrag*nChains))

target_A1=array(rep(0,runsPrag*nChains), dim=c(runsPrag, nChains))
target_A2=target_A1
for (gg in 1:nChains){
 target_A1[,gg]=t(result_vec_byPl_PT[,gg])
 target_A2[,gg]=t(result_vec_byIns_PT[,gg])
}
target_A1A2=rbind(as.vector(target_A1),as.vector(target_A2))
# ensure at least ~20 data points for each of 4 chains
ts=posterior::summarize_draws(posterior::as_draws(t(target_A1A2))) 
data_table_Pl=c(ts[1,]$median,ts[1,]$q5,ts[1,]$q95)
data_table_Ins=c(ts[2,]$median,ts[2,]$q5,ts[2,]$q95)

print(head(target_A1A2[1,]))
print(head(target_A1A2[2,]))

# nb default summarize_draws 90% interval - please adjust for 95%
names(data_table_Pl)=c('50','5','95')
print(data_table_Pl)
names(data_table_Ins)=c('50','5','95')
print(data_table_Ins)



