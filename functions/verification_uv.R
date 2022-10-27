################################################################################
## helper functions to evaluate univariate forecasts (at each lead time)

################################################################################
## CRPS


# function to compute sample CRPS in presence of NAs 
crps_sample_na <- function(y, dat){
  s <- rep(NA, length(y))
  na_ind <- is.na(y)
  s[!na_ind] <- crps_sample(y[!na_ind], dat[!na_ind, ])
  return(s)
}


# wrapper to get CRPS for the three forecast methods
get_crps <- function(obs, ens_raw, pp_mean, pp_sd, clim_mean, clim_sd){
  n_lt <- nrow(obs)
  
  # raw
  crps_raw <- sapply(1:n_lt, function(lt) crps_sample_na(obs[lt, ], ens_raw[lt, , ]))
  
  # pp
  crps_pp <- sapply(1:n_lt, function(lt) crps_norm(obs[lt, ], pp_mean[lt, ], pp_sd[lt, ]))
  
  # clim
  crps_clim <- sapply(1:n_lt, function(lt) crps_norm(obs[lt, ], clim_mean[lt, ], clim_sd[lt, ]))
  
  return(list(raw = crps_raw, pp = crps_pp, clim = crps_clim))
}


# function to plot CRPS as a function of lead time
plot_crps_vs_lt <- function(crps_raw, crps_pp, crps_clim){
  n_lt <- ncol(crps_raw)
  n_na <- apply(crps_raw, 2, function(x){sum(!is.na(x))})
  
  mean_score <- c(colMeans(crps_clim, na.rm = T), 
                  colMeans(crps_raw, na.rm = T),
                  colMeans(crps_pp, na.rm = T))
  
  sd_score <- c(apply(crps_clim, 2, sd, na.rm = T)/sqrt(n_na),
                apply(crps_raw, 2, sd, na.rm=T)/sqrt(n_na),
                apply(crps_pp, 2, sd, na.rm=T)/sqrt(n_na))
  
  # print scores
  score_mat <- matrix(mean_score, nrow = 3, byrow = T)
  colnames(score_mat) <- 1:3; rownames(score_mat) <- c("Clim.", "COSMO", "PP")
  print(score_mat)
  
  # plot scores
  df <- data.frame(s = mean_score,
                   se = sd_score,
                   lead = rep(1:n_lt, 3), mth = rep(c("Clim.", "COSMO", "PP"), each = n_lt))
  ggplot(df) + geom_line(aes(x = lead, y = s, col = mth), size = 0.8) +
    geom_ribbon(aes(x = lead, ymin = s - 1.96*se, ymax = s + 1.96*se, fill = mth), alpha = 0.3, col = NA) +
    scale_y_continuous(name = "CRPS", limits = c(0.75, 4)) + 
    scale_x_continuous(name = "Lead time (days)", breaks = 1:n_lt) +
    theme_bw() +
    theme(legend.justification = c(1, 0.5), legend.position = c(0.99, 0.5),
          legend.title = element_blank())
}


# wrapper to plot CRPS for the three forecast methods
plot_crps <- function(obs, ens_raw, pp_mean, pp_sd, clim_mean, clim_sd){
  crps_output <- get_crps(obs, ens_raw, pp_mean, pp_sd, clim_mean, clim_sd)
  plot_crps_vs_lt(crps_output$raw, crps_output$pp, crps_output$clim)
}



################################################################################
## PIT and rank histograms


# function to get rank of the observation in an ensemble forecast
get_rank <- function(obs, ens){
  r <- sapply(seq_along(obs), function(i) rank(c(obs[i], ens[i, ]))[1])
  return(r)
}


# wrapper to get PIT values for the three forecast methods
get_pit <- function(obs, ens_raw, pp_mean, pp_sd, clim_mean, clim_sd){
  
  n_lt <- nrow(obs)
  n_bins <- dim(ens_raw)[3] + 1
  
  # raw
  pit_raw <- sapply(1:n_lt, function(lt) get_rank(obs[lt, ], ens_raw[lt, , ]))
  
  # pp
  pit_pp <- sapply(1:n_lt, function(lt) pnorm(obs[lt, ], pp_mean[lt, ], pp_sd[lt, ]))
  pit_pp <- floor(pit_pp*n_bins) + 1; pit_pp[pit_pp > n_bins] <- n_bins # convert to ranks
  
  # clim
  pit_clim <- sapply(1:n_lt, function(lt) pnorm(obs[lt, ], clim_mean[lt, ], clim_sd[lt, ]))
  pit_clim <- floor(pit_clim*n_bins) + 1; pit_clim[pit_clim > n_bins] <- n_bins # convert to ranks
  
  return(list(raw = pit_raw, pp = pit_pp, clim = pit_clim))
}


# function to plot PIT histograms at a given lead time
plot_pit_lt <- function(pit, lead, n_bins, title = NULL, ymax = 0.3){

  rank_freq <- sapply(1:n_bins, function(i) mean(pit[, lead] == i, na.rm = T))
  df <- data.frame(freq = rank_freq, rank = as.factor(1:n_bins))
  ggplot(df, aes(x = rank, y = freq)) + geom_bar(stat = "identity") +
    geom_hline(aes(yintercept = 1/n_bins), col="red", lty = "dashed") + 
    scale_x_discrete(name = "Rank") +
    scale_y_continuous(name = "Rel. Freq.", limits = c(0, ymax), expand = c(0, 0)) +
    theme_bw() +
    theme(legend.title = element_blank(), panel.grid = element_blank()) +
    ggtitle(title)
}


# wrapper to plot PIT histograms for the three forecast methods
plot_pit <- function(obs, ens_raw, pp_mean, pp_sd, clim_mean, clim_sd, lead){
  
  n_bins <- dim(ens_raw)[3] + 1
  
  pit_output <- get_pit(obs, ens_raw, pp_mean, pp_sd, clim_mean, clim_sd)

  plot_raw <- plot_pit_lt(pit_output$raw, lead, n_bins, title = "COSMO")
  plot_pp <- plot_pit_lt(pit_output$pp, lead, n_bins, title = "Post-processed")
  plot_clim <- plot_pit_lt(pit_output$clim, lead, n_bins, title = "Climatology")
  
  gridExtra::grid.arrange(plot_clim, plot_raw, plot_pp, nrow = 3)
}



################################################################################
## reliability index

# function to calculate reliability index from ranks
rel_index <- function(pit, n_bins){
  rank_freq <- sapply(1:n_bins, function(i) mean(pit == i, na.rm = T))
  rel_index <- sum(abs(rank_freq - (1/n_bins)))
  return(rel_index)
}


# wrapper to extract reliability indices for the three forecast methods
get_rel_index <- function(obs, ens_raw, pp_mean, pp_sd, clim_mean, clim_sd, n_bins){
  n_lt <- nrow(obs)
  pit_output <- get_pit(obs, ens_raw, pp_mean, pp_sd, clim_mean, clim_sd)
  ri_raw <- sapply(1:n_lt, function(lt) rel_index(pit_output$raw[, lt], n_bins))
  ri_pp <- sapply(1:n_lt, function(lt) rel_index(pit_output$pp[, lt], n_bins))
  ri_clim <- sapply(1:n_lt, function(lt) rel_index(pit_output$clim[, lt], n_bins))
  return(list(raw = ri_raw, pp = ri_pp, clim = ri_clim))
}


# function to plot the reliability index as a function of lead time
plot_rel_index_vs_lt <- function(ri_raw, ri_pp, ri_clim){
  n_lt <- length(ri_raw)
  mean_score <- c(ri_clim, ri_raw, ri_pp)
  
  # print scores
  score_mat <- matrix(mean_score, nrow = 3, byrow = T)
  colnames(score_mat) <- 1:3; rownames(score_mat) <- c("Clim.", "COSMO", "PP")
  print(score_mat)
  
  # plot scores
  df <- data.frame(s = mean_score, lead = rep(1:n_lt, 3), mth = rep(c("Clim.", "COSMO", "PP"), each = n_lt))
  ggplot(df) + geom_line(aes(x = lead, y = s, col = mth), size = 1.2) +
    scale_y_continuous(name = "Rel. Index", limits = c(0, 1)) + 
    scale_x_continuous(name = "Lead time (days)", breaks = 1:n_lt) +
    theme_bw() +
    theme(legend.justification = c(1, 1), legend.position = c(0.99, 0.99),
          legend.title = element_blank())
}


# wrapper to plot reliability index for the three forecast methods
plot_rel_index <- function(obs, ens_raw, pp_mean, pp_sd, clim_mean, clim_sd){
  n_bins <- dim(ens_raw)[3] + 1
  ri_output <- get_rel_index(obs, ens_raw, pp_mean, pp_sd, clim_mean, clim_sd, n_bins)
  plot_rel_index_vs_lt(ri_output$raw, ri_output$pp, ri_output$clim)
}



################################################################################
## twCRPS

# function to calculate the twCRPS for a normal distribution
twcrps_norm <- function(y, mean, sd, t){
  crps_cnorm(y = pmax(y, t), mean, sd, lower = t)
}


# wrapper to get twCRPS for the three forecast methods at a given threshold
get_twcrps <- function(obs, ens_raw, pp_mean, pp_sd, clim_mean, clim_sd, t = -Inf){
  n_lt <- nrow(obs)
  
  # raw
  twcrps_raw <- sapply(1:n_lt, function(lt) crps_sample_na(pmax(obs[lt, ], t), pmax(ens_raw[lt, , ], t)))
  
  # pp
  twcrps_pp <- sapply(1:n_lt, function(lt) twcrps_norm(obs[lt, ], pp_mean[lt, ], pp_sd[lt, ], t))
  
  # clim
  twcrps_clim <- sapply(1:n_lt, function(lt) twcrps_norm(obs[lt, ], clim_mean[lt, ], clim_sd[lt, ], t))
  
  return(list(raw = twcrps_raw, pp = twcrps_pp, clim = twcrps_clim))
}


# function to plot weighted CRPSS at a given lead time
plot_wcrpss_vs_lt <- function(wcrps_raw, wcrps_pp, wcrps_clim, t_vec, lead, ylab = "wCRPSS"){
  
  df <- data.frame(s = c(1 - wcrps_clim[lead, ]/wcrps_raw[lead, ], 1 - wcrps_pp[lead, ]/wcrps_raw[lead, ]),
                   t = t_vec, mth = rep(c("Clim.", "Post-proc."), each = length(t_vec)))
  ggplot(df) + geom_line(aes(x = t, y = s, col = mth), size = 1) + 
    geom_hline(aes(yintercept = 0), lty = "dotted") + 
    geom_vline(aes(xintercept = 25), col = "grey", lty = "dashed") +
    geom_vline(aes(xintercept = 27), col = "grey", lty = "dashed") +
    scale_x_continuous(name = "Threshold (\u00B0C)") +
    scale_y_continuous(name = ylab, limits = c(-2.3, 1)) + 
    theme_bw() +
    theme(legend.justification = c(0, 1), legend.position = c(0.01, 0.99), 
          panel.grid = element_blank(), legend.title = element_blank())
}


# wrapper to get weighted CRPS for the three forecast methods at a range of thresholds
get_wcrps <- function(obs, ens_raw, pp_mean, pp_sd, clim_mean, clim_sd, t_vec, score = "twCRPS"){
  
  if(score == "owCRPS"){
    score_func <- get_owcrps
  }else if(score == "vrCRPS"){
    score_func <- get_vrcrps
  }else{
    score_func <- get_twcrps
  }
  
  output <- lapply(t_vec, function(t) score_func(obs, ens_raw, pp_mean, pp_sd, clim_mean, clim_sd, t))
  
  # raw
  wcrps_raw <- sapply(seq_along(output), function(i) colMeans(output[[i]]$raw, na.rm = T))
  
  # pp
  wcrps_pp <- sapply(seq_along(output), function(i) colMeans(output[[i]]$pp, na.rm = T))
  
  # clim
  wcrps_clim <- sapply(seq_along(output), function(i) colMeans(output[[i]]$clim, na.rm = T))
  
  return(list(raw = wcrps_raw, pp = wcrps_pp, clim = wcrps_clim))
}


# wrapper to plot twCRPSS for the three forecast methods
plot_twcrpss <- function(obs, ens_raw, pp_mean, pp_sd, clim_mean, clim_sd, t_vec = seq(-5, 30), lead){
  wcrps_output <- get_wcrps(obs, ens_raw, pp_mean, pp_sd, clim_mean, clim_sd, t_vec, score = "twCRPS")
  plot_wcrpss_vs_lt(wcrps_output$raw, wcrps_output$pp, wcrps_output$clim, t_vec, lead, ylab = "twCRPSS")
}



################################################################################
## owCRPS

# function to calculate the owCRPS for a normal distribution
owcrps_norm <- function(y, mean, sd, t, BS = T){
  w_y <- as.numeric(y > t)
  obj <- crps_tnorm(y, mean, sd, lower = t)*w_y
  
  if(BS){
    p <- 1 - pnorm(t, mean, sd)
    obj <- obj + (p - w_y)^2
  }
  
  return(obj)
}


# wrapper to get owCRPS for the three forecast methods at a given threshold
get_owcrps <- function(obs, ens_raw, pp_mean, pp_sd, clim_mean, clim_sd, t = -Inf){
  n_lt <- nrow(obs)
  
  # raw
  owcrps_raw <- sapply(1:n_lt, function(lt){
    m <- rowMeans(ens_raw[lt, ,])
    s <- apply(ens_raw[lt, ,], 1, sd)
    owcrps <- owcrps_norm(obs[lt, ], m, s, t)
    return(owcrps)
  })
  
  # pp
  owcrps_pp <- sapply(1:n_lt, function(lt) owcrps_norm(obs[lt, ], pp_mean[lt, ], pp_sd[lt, ], t))
  
  # clim
  owcrps_clim <- sapply(1:n_lt, function(lt) owcrps_norm(obs[lt, ], clim_mean[lt, ], clim_sd[lt, ], t))
  
  return(list(raw = owcrps_raw, pp = owcrps_pp, clim = owcrps_clim))
}


# wrapper to plot owCRPSS for the three forecast methods
plot_owcrpss <- function(obs, ens_raw, pp_mean, pp_sd, clim_mean, clim_sd, t_vec = seq(-5, 30), lead){
  wcrps_output <- get_wcrps(obs, ens_raw, pp_mean, pp_sd, clim_mean, clim_sd, t_vec, score = "owCRPS")
  plot_wcrpss_vs_lt(wcrps_output$raw, wcrps_output$pp, wcrps_output$clim, t_vec, lead, ylab = "owCRPSS")
}



################################################################################
## vrCRPS

# function to calculate the owCRPS for a normal distribution
vrcrps_norm <- function(y, mean, sd, t){
  y_star <- pmax(y, t)
  
  w_y <- as.numeric(y > t)
  
  Ft <- pnorm(t, mean = mean, sd = sd)
  Fys <- pnorm(y_star, mean = mean, sd = sd)
  F0 <- pnorm(0, mean = mean, sd = sd)
  Ft_sq2 <- pnorm(t, mean = mean, sd = sd/sqrt(2))
  dt <- dnorm(t, mean = mean, sd = sd)
  dys <- dnorm(y_star, mean = mean, sd = sd)
  d0 <- dnorm(0, mean = mean, sd = sd)
  
  o1 <- (y - mean)*(2*Fys - 1 - Ft)
  o2 <- (sd^2)*(2*dys - dt)
  o3 <- sd*(1 - Ft_sq2)/sqrt(pi)
  o4 <- (sd^2)*dt*(1 - Ft)
  
  obj1 <- w_y*(o1 + o2) - o3 + o4
  if(t <= 0){
    obj2_1 <- mean*(1 - Ft - 2*F0) + (sd^2)*(2*d0 - dt)
  }else{
    obj2_1 <- mean*(1 - Ft) + (sd^2)*dt
  }
  
  obj2 <- (obj2_1 - abs(y)*w_y)*((1 - Ft) - w_y)
  
  obj <- obj1 + obj2
  
  return(obj)
}
verify_vrcrps <- function(){
  
  mean <- runif(1, -5, 5)
  sd <- runif(1, 1, 3)
  y <- runif(1, -3, 3)
  t <- runif(1, -2, 2)
  
  y_star <- pmax(y, t)
  Ft <- pnorm(t, mean = mean, sd = sd)
  Fys <- pnorm(y_star, mean = mean, sd = sd)
  F0 <- pnorm(0, mean = mean, sd = sd)
  Ft_sq2 <- pnorm(t, mean = mean, sd = sd/sqrt(2))
  dt <- dnorm(t, mean = mean, sd = sd)
  dys <- dnorm(y_star, mean = mean, sd = sd)
  d0 <- dnorm(0, mean = mean, sd = sd)
  
  o1 <- (y - mean)*(2*Fys - 1 - Ft)
  o2 <- (sd^2)*(2*dys - dt)
  obj1 <- o1 + o2
  
  obj1_func <- function(x, y, t, mean, sd) abs(x - y)*(x > t)*dnorm(x, mean, sd)
  print(c(obj1, integrate(obj1_func, lower = -Inf, upper = Inf, y = y, t = t, mean = mean, sd = sd)$val))
  
  
  o3 <- sd*(1 - Ft_sq2)/sqrt(pi)
  o4 <- (sd^2)*dt*(1 - Ft)
  obj2 <- o3 - o4
  
  
  obj2_inner_func <- function(x1, mean, sd) x1*dnorm(x1, mean, sd)
  obj2_inner_int = function(x2, mean, sd){
    sapply(x2, function(z){integrate(obj2_inner_func, lower = z, upper = Inf, mean = mean, sd = sd)$value*dnorm(z, mean, sd)})}
  obj2_1 <- integrate(obj2_inner_int, lower = t, upper = Inf, mean = mean, sd = sd)$val
  
  obj2_inner_func <- function(x1, mean, sd) dnorm(x1, mean, sd)
  obj2_inner_int = function(x2, mean, sd){
    sapply(x2, function(z){integrate(obj2_inner_func, lower = z, upper = Inf, mean = mean, sd = sd)$value*z*dnorm(z, mean, sd)})}
  obj2_2 <- integrate(obj2_inner_int, lower = t, upper = Inf, mean = mean, sd = sd)$val
  
  print(c(obj2, obj2_1 - obj2_2))
  
  y <- runif(10000, -3, 3)
  ens <- t(replicate(length(y), qnorm(seq(0.0001, 0.9999, 0.0001), mean, sd)))
  vrs_samp <- crps_sample(y*(y > t), ens*(ens > t))
  vrs <- vrcrps_norm(y, mean, sd, t)
  print(c(mean(vrs_samp), mean(vrs)))
}


# wrapper to get vrCRPS for the three forecast methods at a given threshold
get_vrcrps <- function(obs, ens_raw, pp_mean, pp_sd, clim_mean, clim_sd, t = -Inf){
  n_lt <- nrow(obs)
  
  # raw
  vrcrps_raw <- sapply(1:n_lt, function(lt) crps_sample_na(obs[lt, ]*(obs[lt, ] > t), ens_raw[lt, , ]*(ens_raw[lt, , ] > t)))
  
  # pp
  vrcrps_pp <- sapply(1:n_lt, function(lt) vrcrps_norm(obs[lt, ], pp_mean[lt, ], pp_sd[lt, ], t))
  
  # clim
  vrcrps_clim <- sapply(1:n_lt, function(lt) vrcrps_norm(obs[lt, ], clim_mean[lt, ], clim_sd[lt, ], t))
  
  return(list(raw = vrcrps_raw, pp = vrcrps_pp, clim = vrcrps_clim))
}


# wrapper to plot vrCRPSS for the three forecast methods
plot_vrcrpss <- function(obs, ens_raw, pp_mean, pp_sd, clim_mean, clim_sd, t_vec = seq(-5, 30), lead){
  wcrps_output <- get_wcrps(obs, ens_raw, pp_mean, pp_sd, clim_mean, clim_sd, t_vec, score = "vrCRPS")
  plot_wcrpss_vs_lt(wcrps_output$raw, wcrps_output$pp, wcrps_output$clim, t_vec, lead, ylab = "vrCRPSS")
}



################################################################################
## conditional PIT histograms and PIT reliability diagrams

# function to calculate conditional PIT value for a normal distribution
norm_cpit <- function(o, m, s, t){
  
  pit <- (pnorm(o, mean = m, sd = s) - pnorm(t, mean = m, sd = s))/
    (1 - pnorm(t, mean = m, sd = s))
  pit[o < t] <- NA
  
  return(pit)
}


# wrapper to get conditional PIT values for the three forecast methods
get_cpit <- function(obs, ens_raw, pp_mean, pp_sd, clim_mean, clim_sd, t = -Inf){
  n_lt <- nrow(obs)
  
  # raw
  cpit_raw <- sapply(1:n_lt, function(lt){
    m <- rowMeans(ens_raw[lt, ,])
    s <- apply(ens_raw[lt, ,], 1, sd)
    cpit <- norm_cpit(obs[lt, ], m, s, t)
    return(cpit)
  })
  
  # pp
  cpit_pp <- sapply(1:n_lt, function(lt) norm_cpit(obs[lt, ], pp_mean[lt, ], pp_sd[lt, ], t))
  
  # clim
  cpit_clim <- sapply(1:n_lt, function(lt) norm_cpit(obs[lt, ], clim_mean[lt, ], clim_sd[lt, ], t))
  
  return(list(raw = cpit_raw, pp = cpit_pp, clim = cpit_clim))
}


# function to convert conditional PIT values to ranks
get_cpit_ranks <- function(cpit, n_bins, t = -Inf){
  cpit <- floor(cpit*n_bins) + 1
  cpit[cpit > n_bins] <- n_bins
  return(cpit)
}


# wrapper to plot conditional PIT histograms for the three forecast methods
plot_cpit_hist <- function(obs, ens_raw, pp_mean, pp_sd, clim_mean, clim_sd, t, n_bins, lead){

  cpit_output <- get_cpit(obs, ens_raw, pp_mean, pp_sd, clim_mean, clim_sd, t)
  
  cpit_raw <- get_cpit_ranks(cpit_output$raw, n_bins)
  cpit_pp <- get_cpit_ranks(cpit_output$pp, n_bins)
  cpit_clim <- get_cpit_ranks(cpit_output$clim, n_bins)
  
  plot_raw <- plot_pit_lt(cpit_raw, lead, n_bins, title = paste("COSMO: T >", t), ymax = 0.4)
  plot_pp <- plot_pit_lt(cpit_pp, lead, n_bins, title = paste("Post-processed: T >", t), ymax = 0.4)
  plot_clim <- plot_pit_lt(cpit_clim, lead, n_bins, title = paste("Climatology: T >", t), ymax = 0.4)
  
  gridExtra::grid.arrange(plot_clim, plot_raw, plot_pp, nrow = 3)
}


# function to plot PIT reliability diagrams (from https://github.com/resinj/replication_GR21)
PITreldiag <- function(pit, resampling = TRUE, n_resamples = 1000, region_level = 0.9){
  # auxiliary function
  color_step_poly <- function(lim,x_up,x_low,y_up,y_low,col = "lightblue1"){
    polygon(c(lim[1],rbind(x_up,c(x_up[-1],lim[2])),lim[2],rbind(rev(x_low),rev(x_low)),lim[1]),
            c(lim[1],c(rbind(y_up,y_up)),lim[2],rbind(rev(y_low),c(rev(y_low)[-1],lim[1])),lim[1]),
            border = NA,col = "lightblue1")
  }
  
  dist_PIT = ecdf(pit)
  
  par(mgp = c(-1.5, 1, 0))
  plot(NULL,xlim = c(0,1),ylim = c(0,1),main = "", 
       xlab = expression(z),ylab = expression("fraction of PIT-values" <= z))
  
  if(resampling){
    low = floor(n_resamples * (1-region_level)/2)
    up = n_resamples - low
    
    resamples = sapply(1:n_resamples,function(i) runif(length(pit)))
    
    t = seq(0,1,0.001)
    dist_resamples_t = apply(resamples, 2, function(s) ecdf(s)(t))
    dist_resamples_t_sorted = apply(dist_resamples_t,1,sort)
    
    color_step_poly(c(0,1),t,t,dist_resamples_t_sorted[up,],dist_resamples_t_sorted[low,])
    # polygon(c(0,t,1,rev(t),0),
    #         c(0,dist_resamples_t_sorted[low,],1,rev(dist_resamples_t_sorted[up,]),0),border = NA,col = "lightblue1")
    points(t, dist_resamples_t_sorted[low,],type = "l",lty = 1,col = "lightblue2")
    points(t, dist_resamples_t_sorted[up,],type = "l",lty = 1,col = "lightblue2")
  }
  
  points(c(0,1),c(0,1),type = "l",lty = 2,col = "grey")
  plot(dist_PIT,do.points = FALSE,xlim = c(0,1),col.01line = NULL,verticals = TRUE,add=TRUE)
}


# function to plot conditional reliability diagram at a given lead time
plot_cpit_reldiag_lt <- function(cpit, lead, level = 0.99){
  PITreldiag(na.omit(as.vector(cpit[, lead])), region_level = level)
}


# wrapper to plot conditional reliability diagrams for the three forecast methods
plot_cpit_reldiag <- function(obs, ens_raw, pp_mean, pp_sd, clim_mean, clim_sd, t, lead){
  cpit_output <- get_cpit(obs, ens_raw, pp_mean, pp_sd, clim_mean, clim_sd, t)
  
  par(mfrow = c(3, 1))  
  plot_cpit_reldiag_lt(cpit_output$clim, lead, level = 0.99)
  plot_cpit_reldiag_lt(cpit_output$raw, lead, level = 0.99)
  plot_cpit_reldiag_lt(cpit_output$pp, lead, level = 0.99)
}



################################################################################
## reliability diagrams

# wrapper to get probability of exceeding a threshold t for the three forecast methods
get_exc_p <- function(obs, ens_raw, pp_mean, pp_sd, clim_mean, clim_sd, t, lead){
  n_lt <- nrow(obs)
  
  # raw
  p_raw <- rowMeans(ens_raw[lead, , ] > t)
  
  # pp
  p_pp <- 1 - pnorm(t, pp_mean[lead, ], pp_sd[lead, ])
  
  # clim
  p_clim <- 1 - pnorm(t, clim_mean[lead, ], clim_sd[lead, ])
  
  return(list(raw = p_raw, pp = p_pp, clim = p_clim))
}


# wrapper to plot reliability diagrams for the three forecast methods
plot_reldiag <- function(obs, ens_raw, pp_mean, pp_sd, clim_mean, clim_sd, t, lead){
  p_output <- get_exc_p(obs, ens_raw, pp_mean, pp_sd, clim_mean, clim_sd, t, lead)
  exc_obs <- as.numeric(obs[lead, ] > t)
  ind <- !is.na(exc_obs)
  reliabilitydiag(as.data.frame(p_output)[ind, ], y = exc_obs[ind], xtype = "continuous")
}


