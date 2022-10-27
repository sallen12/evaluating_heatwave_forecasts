
setwd("~/R/Heatwaves")
source('./functions/ConfigureEnv.R')
ConfigureEnv()


################################################################################
## load toy data

load("./input_data/ToyData.RData")


################################################################################
## evaluate overall forecast performance

source('./functions/verification_uv.R')
source('./functions/verification_mv.R')

##### accuracy

## CRPS
plot_crps(obs_dat, ens_raw, pp_mean, pp_sd, clim_mean, clim_sd)


## ES
print_mv_s(obs_dat, ens_raw, ens_pp, ens_clim, score = "ES")


## VS
print_mv_s(obs_dat, ens_raw, ens_pp, ens_clim, score = "VS")


##### calibration

## PIT and rank histograms
plot_pit(obs_dat, ens_raw, pp_mean, pp_sd, clim_mean, clim_sd, lead = 3)


## reliability index
plot_rel_index(obs_dat, ens_raw, pp_mean, pp_sd, clim_mean, clim_sd)



################################################################################
## evaluate forecasts for heatwave severity

##### accuracy

## twCRPS
plot_twcrpss(obs_dat, ens_raw, pp_mean, pp_sd, clim_mean, clim_sd, t_vec = seq(-5, 30), lead = 3)


## owCRPS
plot_owcrpss(obs_dat, ens_raw, pp_mean, pp_sd, clim_mean, clim_sd, t_vec = seq(-5, 30), lead = 3)


## vrCRPS
plot_vrcrpss(obs_dat, ens_raw, pp_mean, pp_sd, clim_mean, clim_sd, t_vec = seq(-5, 30), lead = 3)


## twES
plot_tw_mv_s(obs_dat, ens_raw, ens_pp, ens_clim, score = "twES")


## twVS
plot_tw_mv_s(obs_dat, ens_raw, ens_pp, ens_clim, score = "twVS")


##### calibration

## cPIT histograms
plot_cpit_hist(obs_dat, ens_raw, pp_mean, pp_sd, clim_mean, clim_sd, t = 20, n_bins = 10, lead = 3)


## cPIT reliability diagrams
plot_cpit_reldiag(obs_dat, ens_raw, pp_mean, pp_sd, clim_mean, clim_sd, t = 20, lead = 3)


## reliability diagrams
plot_reldiag(obs_dat, ens_raw, pp_mean, pp_sd, clim_mean, clim_sd, t = 20, lead = 3)
