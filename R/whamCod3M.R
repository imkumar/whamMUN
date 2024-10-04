library(wham)

load(file='data-raw/data3M.rda')
list2env(data3M, envir=.GlobalEnv)

waa <- array(0,dim=c(2,Y,A))
waa[1,,] <- wstock
waa[2,,] <- wcatch
mat <- array(mat,dim=c(1,Y,A))
Maa <- array(Msca,dim=c(1,1,Y,A))
waa_pointer_indices <- 1#index for stock weight
waa_pointer_fleets <- 2 #index for catch weight
waa_pointer_ssb <- 1 #using stock weight for ssb
fracyr_ssb <-  matrix(0,c(Y,1)) #fraction of year spawing occurs. we assume jan=0

#prepare catch data
cAA <- exp(logC) #log Catch numbers at age
cpa <- cAA/rowSums(cAA,na.rm=TRUE)
cpa[is.na(cpa)] <- 0

catch <- as.matrix(logCton)|>exp() #log total catch in tonnes
catch_cv <- matrix(.1,Y,1) #asssumed .1
catch_Neff <- matrix(100,Y,1) #asssumed 
catch_paa <- array(cpa,dim=c(1,Y,A))
use_catch_paa <- rowSums(cpa)|>as.matrix() #if equal to 1 only then fit paa
selblock_pointer_fleets <- matrix(1,Y,1) #assumed 1 for fleet selectivity

#prepare index data
iAA <- exp(logCPUE.EU) #log EU survey index
ipa <- iAA/rowSums(iAA,na.rm=TRUE)
ipa[is.na(ipa)] <- 0

index <- as.matrix(rowSums(exp(logCPUE.EU),na.rm=TRUE))#age-aggregate index
index_cv <- matrix(.3,Y,1) #assumed .3
index_Neff <- matrix(50,Y,1)
index_fracyr <- matrix(.5,Y,1)#assuming survey in june=6/12
units_indices <- 2 ; #2=index in number, 1=index in biomass
units_index_paa <- 2;# 2=proportion at age for index in number, 1=prop at age for index in weight
use_indices <-matrix(1,Y,1);#to use aggregated index or not. 1=Yes, 0=No
use_index_paa <- rowSums(cpa)|>as.matrix(); #same as above for proportion
index_paa <- array(ipa,dim=c(1,Y,A)); # index proportion at age
selblock_pointer_indices <- matrix(2,Y,1) #assumed 2 for index selectivity 

# prepare env data
comp3M <- read.csv('data-raw/composite_3M_1988.csv')

#prepare basic, catch, and index info
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

#Various model options with different selectivity, M, NAA, and recruit options
#M options
M1 <- list(initial_MAA = array(0.2, dim = c(1,1,Y,A))) #M constant 0.2

M2 <- list(initial_MAA = array(Maa,dim = c(1,1,Y,A))) #M constant from SCAA

M3 <- list(mean_model="fixed-M",re_model=matrix("2dar1"),
           initial_means=array(colMeans(Maa[1,1,,]),dim=c(1,1,A))) #estimate ar1Md 

M4 <- list(mean_model="weight-at-age",b_prior=TRUE) #Lorenzen M

F_opts <- list(
  F = cbind(rep(5,Y)),
  map_F = cbind(rep(NA, Y)))


#Sel options
sel1 <- list(
  model = rep("age-specific",2), n_selblocks = 2, #age specific mean selectiviy for fleet and index. constant over years
  initial_pars = list(
    c(0,rep(0.5,3),0.5,0.5,0.5,0.5), #fixed age 5=1
    c(0.5,0.5,0.5,0.5,0.5, 1,0.5,0.5)),
  map_pars = list(c(NA,2:5,6,6,6),c(9:13,NA,15,16)),
  re = c("2dar1", "none"),
  fix_pars = NULL)


sel2 <- list(
  model = c("age-specific","logistic"), n_selblocks = 2, #age specific mean selectiviy for fleet and index. constant over years
  fix_pars = list(NULL,NULL), #above will be estimated
  initial_pars = list(
    c(0,rep(0.5,3),0.5,0.5,0.5,0.5), #fixed age 5=1
    c(2,0.2)),
  map_pars <- list(c(NA,2:5,6,6,6),c(9:10)),
  re = c("2dar1", "none"),
  fix_pars = NULL
)

sel3 <- list(
  model = c("age-specific","logistic"), n_selblocks = 2, #age specific mean selectiviy for fleet and index. constant over years
  initial_pars = list(
    c(0,rep(0.5,3),1,0.5,0.5,0.5), #fixed age 5=1
    c(2,0.2)),
  map_pars = list(c(NA,NA,3:4,NA,6,6,6),c(9:10)),
  re = c("2dar1", "none"),
  fix_pars = NULL
)

sel4 <- list(
  model = c("age-specific","logistic"), n_selblocks = 2, #age specific mean selectiviy for fleet and index. constant over years
  initial_pars = list(
    c(0,rep(0.5,3),1,0.5,0.5,0.5), #fixed age 5=1
    c(2,0.2)),
  map_pars = list(c(NA,NA,3:4,NA,6,6,6),c(9:10)),
  re = c("iid", "none"),
  fix_pars = NULL)

NAA1 <- list(N1_model="equilibrium",
             recruit_model=2,
             sigma="rec+1", #1 sigma for rec and 1 sigma for NAA
             cor="iid") # cor is for NAA only


NAA2 <- list(N1_model="equilibrium",
             recruit_model=2,
             sigma="rec+1", #1 sigma for rec and 1 sigma for NAA
             cor="2dar1") # cor is for NAA only

NAA3 <- list(N1_model="equilibrium",
             recruit_model=2,
             sigma="rec+1", #1 sigma for rec and 1 sigma for NAA
             cor="2dar1",
             cor_map=array(c(NA,NA,1),dim=c(1,1,3))) # cor is for NAA only

NAA4 <- list(N1_model="equilibrium",
             recruit_model=4,
             sigma="rec+1", #1 sigma for rec and 1 sigma for NAA
             cor="2dar1",
             cor_map=array(c(NA,NA,1),dim=c(1,1,3))) # cor is for NAA only

NAA5 <- list(N1_model="equilibrium",
             recruit_model=4,
             sigma="rec+1", #1 sigma for rec and 1 sigma for NAA
             cor="iid",
             decouple_recruitment=T) # cor is for NAA only

input <- prepare_wham_input(basic_info = basic_info,
                            catch_info = catch_info,
                            index_info = index_info)

inputM1 <- set_F(input, F_opts) # this is for sel2. comment out if using sel3
inputM1 <- set_selectivity(inputM1,sel1)
inputM1 <- set_NAA(inputM1,NAA_re = NAA1)
inputM1 <- set_M(inputM1, M1) #M1S1N1

inputM2 <- set_M(inputM1, M2) #M2S1N1


inputM3 <- set_selectivity(inputM2,sel2) #M2S2N1

inputM4 <- set_selectivity(input,sel3) #M2S3N1
inputM4 <- set_NAA(inputM4,NAA_re = NAA1)
inputM4<- set_M(inputM4, M2)

############################################################
inputM5<- set_M(inputM3, M3)
inputM6 <- set_M(inputM3,M4)
inputM7 <- set_NAA(inputM6,NAA_re = NAA2)
inputM8 <- set_NAA(inputM6,NAA_re = NAA3)

inputM11 <- set_selectivity(input,sel4) #
inputM11 <- set_M(inputM11, M4)
inputM11 <- set_NAA(inputM11,NAA_re = NAA5)

inputM9 <- set_NAA(inputM8,NAA_re = NAA4)
inputM10 <- set_NAA(inputM9,NAA_re = NAA5)

ecov <- list(
  label = "comp3M",
  mean = as.matrix(comp3M$value),
  logsigma = 'est_1', # estimate obs sigma, 1 value shared across years
  year = comp3M$year,
  use_obs = matrix(1, ncol=1, nrow=dim(comp3M)[1]), # use all obs (=1)
  process_model = 'ar1', # "rw" or "ar1"
  recruitment_how = matrix("controlling-lag-0-linear")) # n_Ecov x n_stocks

input12 <- set_ecov(inputM9, ecov)
  
fitM1 <- fit_wham(inputM1, do.retro = TRUE, do.osa = FALSE, do.sdrep = TRUE, do.brps = FALSE)
fitM1$gr(fitM1$opt$par)
fitM1$final_gradient

fitM2 <- update(fitM1,input=inputM2)
fitM3 <- update(fitM1,input=inputM3)
fitM4 <- update(fitM1,input=inputM4)
res <- compare_wham_models(list(fitM1,fitM2), fdir = "result1-2")

fitM5 <- update(fitM1,input=inputM5)
fitM6 <- update(fitM1,input=inputM6)
fitM7 <- update(fitM1,input=inputM7)
fitM8 <- update(fitM1,input=inputM8)
fitM9 <- update(fitM1,input=inputM9)
fitM10 <- update(fitM1,input=inputM10)
fitM11 <- update(fitM1,input=inputM11)
check_convergence(fitM12)
fitM12 <- update(fitM1,input=input12)


rep11 <- fitM11$rep
nlls <- grep("nll", names(rep11), value = TRUE)
sapply(nlls, function(x) sum(rep11[[x]]))



res <- compare_wham_models(list(fitM3,fitM5,fitM6), fdir = "result3-5")
res <- compare_wham_models(list(fitM7,fitM8,fitM9,fitM10), fdir = "result1-11")

res <- compare_wham_models(list(fitM1,fitM2,fitM3,fitM4,fitM5,fitM6), fdir = "result")
plot_wham_output(fitM1);file.rename("wham_figures_tables.html","m1.html")
plot_wham_output(fitM2);file.rename("wham_figures_tables.html","m2.html")
plot_wham_output(fitM3);file.rename("wham_figures_tables.html","m3.html")
plot_wham_output(fitM4);file.rename("wham_figures_tables.html","m4.html")

plot_wham_output(fitM5);file.rename("wham_figures_tables.html","m5.html")
plot_wham_output(fitM6);file.rename("wham_figures_tables.html","m6.html")
plot_wham_output(fitM7);file.rename("wham_figures_tables.html","m7.html")
plot_wham_output(fitM8);file.rename("wham_figures_tables.html","m8.html")
plot_wham_output(fitM9);file.rename("wham_figures_tables.html","m9.html")
plot_wham_output(fitM10);file.rename("wham_figures_tables.html","m10.html")
plot_wham_output(fitM11);file.rename("wham_figures_tables.html","m11.html")
plot_wham_output(fitM12);file.rename("wham_figures_tables.html","m12.html")
