################################################################################
## helper functions to evaluate multivariate forecasts (across all lead times)

################################################################################
## ES and VS

# wrapper to get ES or VS for the three forecast methods
get_mv_s <- function(obs, ens_raw, ens_pp, ens_clim, score = "ES"){
  
  if(score == "VS"){
    score_sample <- vs_sample
  }else{
    score_sample <- es_sample
  }
  
  # raw
  s_raw <- sapply(1:ncol(obs), function(i) score_sample(obs[, i], ens_raw[, i, ]))
  
  # pp
  s_pp <- sapply(1:ncol(obs), function(i) score_sample(obs[, i], ens_pp[, i, ]))
  
  # clim
  s_clim <- sapply(1:ncol(obs), function(i) score_sample(obs[, i], ens_clim[, i, ]))
  
  return(list(raw = s_raw, pp = s_pp, clim = s_clim))
}


# print ES or VS for the three forecast methods
print_mv_s <- function(obs, ens_raw, ens_pp, ens_clim, score =  "ES"){
  s_output <- get_mv_s(obs, ens_raw, ens_pp, ens_clim, score)
  score_mean <- c(mean(s_output$clim, na.rm = T), mean(s_output$raw, na.rm = T), mean(s_output$pp, na.rm = T))
  names(score_mean) <- c("Clim.", "COSMO", "PP")
  print(score_mean)
}


################################################################################
## twES and twVS

# function to calculate twES for an ensemble forecast
twes_sample <- function(y, x, w_y, w_x, x_0){
  # x is a numerical array of dimension (d, M)
  # y is a numerical vector of length d
  # w_y is a logical
  # w_x is a logical vector of length d
  # x_0 is a numerical vector of length d
  
  v_y <- w_y*y + (1 - w_y)*x_0
  v_x <- sapply(seq_along(w_x), function(i) w_x[i]*x[, i] + (1 - w_x[i])*x_0)
  
  obj <- es_sample(v_y, v_x)
  
  return(obj)
}


# function to calculate twVS for an ensemble forecast
twvs_sample <- function(y, x, w_y, w_x, x_0){
  # x is a numerical array of dimension (d, M)
  # y is a numerical vector of length d
  # w_y is a logical
  # w_x is a logical vector of length d
  # x_0 is a numerical vector of length d
  
  v_y <- w_y*y + (1 - w_y)*x_0
  v_x <- sapply(seq_along(w_x), function(i) w_x[i]*x[, i] + (1 - w_x[i])*x_0)
  
  obj <- vs_sample(v_y, v_x)
  
  return(obj)
}


# wrapper to calculate twES for a sequence of ensemble forecasts
twes_sample_n <- function(y, x, w_y, w_x, x_0){
  obj <- sapply(1:nrow(y), function(i) twes_sample(y[i, ], x[i, , ], w_y[i], w_x[i, ], x_0))
  return(obj)
}


# wrapper to calculate twVS for a sequence of ensemble forecasts
twvs_sample_n <- function(y, x, w_y, w_x, x_0){
  obj <- sapply(1:nrow(y), function(i) twvs_sample(y[i, ], x[i, , ], w_y[i], w_x[i, ], x_0))
  return(obj)
}


# wrapper to get twES or twVS for the three forecast methods at a given threshold
get_tw_mv_s <- function(obs, ens_raw, ens_pp, ens_clim, weight_func = function(x) T, x_0 = c(0, 0, 0), score = "twES"){
  
  n_lt <- nrow(obs)
    
  if(score == "twVS"){
    score_sample <- twvs_sample_n
  }else{
    score_sample <- twes_sample_n
  }
  
  obs <- t(obs)
  w_o <- apply(obs, 1, weight_func)
  
  # raw
  w_ens <- apply(ens_raw, c(2, 3), weight_func)
  tws_raw <- score_sample(obs, aperm(ens_raw, c(2, 1, 3)), w_o, w_ens, x_0)
  
  # pp
  w_ens <- apply(ens_pp, c(2, 3), weight_func)
  tws_pp <- score_sample(obs, aperm(ens_pp, c(2, 1, 3)), w_o, w_ens, x_0)
  
  # clim
  w_ens <- apply(ens_clim, c(2, 3), weight_func)
  tws_clim <- score_sample(obs, aperm(ens_clim, c(2, 1, 3)), w_o, w_ens, x_0)
  
  return(list(raw = tws_raw, pp = tws_pp, clim = tws_clim))
}


# wrapper to plot twES or twVS for the three forecast methods and four weight functions
plot_tw_mv_s <- function(obs, ens_raw, ens_pp, ens_clim, score = "twES"){
  
  score_means <- matrix(NA, nrow = 5, ncol = 3)
  colnames(score_means) <- c("Clim.", "COSMO", "PP")
  rownames(score_means) <- 0:4
  
  # all events (same as unweighted score)
  s_output <- get_tw_mv_s(obs, ens_raw, ens_pp, ens_clim, score = score)
  score_means[1, ] <- c(mean(s_output$clim, na.rm = T), mean(s_output$raw, na.rm = T), mean(s_output$pp, na.rm = T))
  
  # category 1 heat event
  weight_func <- function(x) all(x < 25)                
  x_0 <- c(25, 25, 25)
  s_output <- get_tw_mv_s(obs, ens_raw, ens_pp, ens_clim, weight_func, x_0, score = score)
  score_means[2, ] <- c(mean(s_output$clim, na.rm = T), mean(s_output$raw, na.rm = T), mean(s_output$pp, na.rm = T))
  
  
  # category 2 heat event
  weight_func <- function(x) sum(x >= 25) %in% c(1, 2)        
  x_0 <- c(25, 25, 25)
  s_output <- get_tw_mv_s(obs, ens_raw, ens_pp, ens_clim, weight_func, x_0, score = score)
  score_means[3, ] <- c(mean(s_output$clim, na.rm = T), mean(s_output$raw, na.rm = T), mean(s_output$pp, na.rm = T))*10
  
  # category 3 heat event
  weight_func <- function(x) all(x >= 25) & any(x < 27)       
  x_0 <- c(25, 25, 25)
  s_output <- get_tw_mv_s(obs, ens_raw, ens_pp, ens_clim, weight_func, x_0, score = score)
  score_means[4, ] <- c(mean(s_output$clim, na.rm = T), mean(s_output$raw, na.rm = T), mean(s_output$pp, na.rm = T))*10
  
  # category 4 heat event
  weight_func <- function(x) all(x >= 27)
  x_0 <- c(27, 27, 27)
  s_output <- get_tw_mv_s(obs, ens_raw, ens_pp, ens_clim, weight_func, x_0, score = score)
  score_means[5, ] <- c(mean(s_output$clim, na.rm = T), mean(s_output$raw, na.rm = T), mean(s_output$pp, na.rm = T))*100
  
  # print scores
  print(score_means)
  
  # plot scores
  df <- data.frame(s = as.vector(score_means),
                   mth = rep(c("Clim.", "COSMO", "PP"), each = 5), 
                   cat = c("All", "Cat1", "Cat2", "Cat3", "Cat4"))
  ggplot(df) + geom_point(aes(x = cat, y = s, col = mth, shape = mth)) + 
    scale_x_discrete(name = NULL) +
    scale_y_continuous(name = score) +
    theme_bw() +
    theme(legend.justification = c(1, 1), legend.position = c(0.99, 0.99),
          legend.title = element_blank())
  
}

