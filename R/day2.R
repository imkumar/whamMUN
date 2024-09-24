library(wham)
list.files("data-raw")
load("data-raw/OutputforOM2023_chain1.RData")
wstock <- wstock
mat <- apply(bayesianogives,c(1,2),median)
dim(mat)


load("data-raw/data3Mcod2023.rda")
names(data)


years <- as.integer(1988:2022)
ages <- 1:8
Y <- data$Y
A <- data$A;A
waa <- array(0,dim=c(2,Y,A))
waa[1,,] <- wstock
waa[2,,] <- data$wcatch
mat <- array(mat,dim=c(1,Y,A))
Maa <- array(.2,dim=c(1,1,Y,A))
waa_pointer_indices <- 1#index for stock weight
waa_pointer_fleets <- 2 #index for catch weight
waa_pointer_ssb <- 1 #using stock weight for ssb
fracyr_ssb <-  matrix(0,c(Y,1)) #fraction of year spawing occurs. we assume jan=0

cAA <- exp(data$logC)
# cAA[is.na(cAA)] <- 0
cpa <- cAA/rowSums(cAA,na.rm=TRUE)
cpa[is.na(cpa)] <- 0




catch <- as.matrix(data$logCton)|>exp()
catch_cv <- matrix(.1,Y,1) #asssumed .1
catch_Neff <- matrix(100,Y,1) #asssumed 
catch_paa <- array(cpa,dim=c(1,Y,A))
use_catch_paa <- rowSums(cpa)|>as.matrix()
selblock_pointer_fleets <- matrix(1,Y,1) #assumed 1 for fleet selectivity

iAA <- exp(data$logCPUE.EU)
ipa <- iAA/rowSums(iAA,na.rm=TRUE)
ipa[is.na(ipa)] <- 0


index <- as.matrix(rowSums(exp(data$logCPUE.EU),na.rm=TRUE))#age-aggregate index
index_cv <- matrix(.3,Y,1) #assumed .3
index_Neff <- matrix(50,Y,1)
index_fracyr <- matrix(.5,Y,1)#assuming survey in june=6/12
units_indices <- 2 ; #2=index in number, 1=index in biomass
units_index_paa <- 2;# 2=proportion at age for index in number, 1=prop at age for index in weight
use_indices <-matrix(1,Y,1);#to use aggregated index or not. 1=Yes, 0=No
use_index_paa <- rowSums(cpa)|>as.matrix(); #same as above for proportion
index_paa <- array(ipa,dim=c(1,Y,A)); # index proportion at age
selblock_pointer_indices <- matrix(2,Y,1) #assumed 2 for index selectivity 



basic_info <- list(
  n_stocks = 1L,
  ages = ages,
  n_seasons = 1L,
  n_fleets = 1L,
  fracyr_SSB = fracyr_ssb,
  maturity = mat,
  years = years,
  waa = waa,
  waa_pointer_ssb = waa_pointer_ssb
)


catch_info <- list(
  n_fleets = NCOL(catch),
  agg_catch = catch,
  agg_catch_cv = catch_cv,
  catch_paa = catch_paa,
  use_catch_paa = use_catch_paa,
  catch_Neff = catch_Neff,
  selblock_pointer_fleets = selblock_pointer_fleets,
  waa_pointer_fleets = waa_pointer_fleets
)


index_info <- list(
  n_indices = NCOL(index),
  agg_indices = index,
  units_indices = units_indices,
  units_index_paa = units_index_paa,
  agg_index_cv = index_cv,
  fracyr_indices = index_fracyr,
  use_indices = use_indices,
  use_index_paa = use_index_paa,
  index_paa = index_paa,
  index_Neff = index_Neff,
  selblock_pointer_indices = selblock_pointer_indices,
  waa_pointer_indices = waa_pointer_indices
)


M_in <- list(
  initial_MAA = array(Maa, dim = c(1,1,length(basic_info$years), length(basic_info$ages)))
)
#
# fit 3 logistic selevtity: 1: 1 and 3 for 1 fleet(98-2005 and 2006-2022) and 2 for index selectivity. If fix_pars=NULL then slectivity is estimated.
#intital pars for the 3 blocks.



selectivity <- list(
  model = rep("age-specific",2), n_selblocks = 2, #age specific mean selectiviy for fleet and index. constant over years
  fix_pars = list(NULL,NULL) #abouve will be estimated
)

# alternative mean models
# selectivity <- list(model = c("double-logistic", "decreasing-logistic", "logistic"))
# selectivity$initial_pars <- list(
#   c(2,0.2,5,0.2), #ascending a50, 1/slope, and descending a50, 1/slope
#   c(3, 0.2), #descending a50 and 1/slope
#   c(3, 0.2)) #ascending a50, 1/slope
# selectivity$fix_pars <- list(NULL, NULL, NULL) #all estimated
# input_0 <- set_selectivity(input, selectivity = selectivity)
# nofit_0 <- fit_wham(input_0, do.fit = FALSE)
# plot_wham_output(nofit_0, dir.main = tmp.dir)

F_opts <- list(
  F = cbind(rep(5,Y)),
  map_F = cbind(rep(NA, Y)))




# selectivity$re <- c("ar1_y", "none")
selectivity$re <- c("2dar1", "none")
selectivity$fix_pars <- NULL
selectivity$initial_pars <- list(
  c(rep(0.5,4),0.5,0.5,0.5,0.5), #fixed age 5=1
  c(0.5,0.5,0.5,0.5,0.5, 0.5,0.5,0.5))
selectivity$map_pars <- list(c(1:5,6,6,6),c(9:16))

input_all <- prepare_wham_input(basic_info = basic_info, selectivity = selectivity, catch_info = catch_info, index_info = index_info, M = M_in)
input_all <- set_NAA(input_all,NAA_re = list(N1_model="equilibrium"))
input_all <- set_F(input_all, F_opts)

nofit_all <- fit_wham(input_all, do.fit = FALSE, do.brps = FALSE)
nofit_all$gr()

fit_all <- fit_wham(input_all, do.retro = FALSE, do.osa = FALSE, do.sdrep = FALSE)
fit_all$opt
fit_all$gr(fit_all$opt$par)|>round(3)

plot_wham_output(fit_all)
fit_all <- fit_wham(input_all, do.retro = FALSE, do.osa = FALSE, do.sdrep = TRUE)
fit_all$opt
fit_all$gr(fit_all$opt$par)|>round(3)


