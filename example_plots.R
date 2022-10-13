################################################################################
# code to plot the examples of weighted verification tools in Figures 1-4


require(ggplot2)
require(reliabilitydiag)
require(scoringRules)


# function to plot PIT histograms
PIThist <- function(y, bins = 10, ymax = 0.2, title = NULL){
  
  y <- na.omit(y)
  
  intervals <- seq(0, 1, length = bins + 1)
  mids <- (intervals[1:bins] + intervals[2:(bins + 1)])/2
  counts <- sapply(1:bins, function(i) sum(y >= intervals[i] & y < intervals[i + 1]))
  counts[bins] <- counts[bins] + sum(y == 1)
  counts <- counts/sum(counts)
  
  df <- data.frame(c = counts, rank = mids)
  ggplot(df) + geom_bar(stat = "identity", aes(x = rank, y = c)) + 
    geom_hline(aes(yintercept = 1/bins), col="red", lty = "dashed") + 
    scale_x_continuous(name = "PIT value", expand = c(0.005, 0.005)) + 
    scale_y_continuous(name = "Rel. Freq.", limits = c(0, ymax), expand = c(0, 0)) +
    ggtitle(title) + 
    theme_bw() +
    theme(legend.title = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank())
}


# function to plot PIT reliability diagrams (from https://github.com/resinj/replication_GR21)
PITreldiag <- function(pit, y, resampling = TRUE, n_resamples = 1000, region_level = 0.9){
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
    
    resamples = sapply(1:n_resamples,function(i) runif(length(y)))
    
    t = seq(0,1,0.001)
    dist_resamples_t = apply(resamples, 2, function(s) ecdf(s)(t))
    dist_resamples_t_sorted = apply(dist_resamples_t,1,sort)
    
    color_step_poly(c(0,1),t,t,dist_resamples_t_sorted[up,],dist_resamples_t_sorted[low,])
    points(t, dist_resamples_t_sorted[low,],type = "l",lty = 1,col = "lightblue2")
    points(t, dist_resamples_t_sorted[up,],type = "l",lty = 1,col = "lightblue2")
  }
  
  points(c(0,1),c(0,1),type = "l",lty = 2,col = "grey")
  plot(dist_PIT,do.points = FALSE,xlim = c(0,1),col.01line = NULL,verticals = TRUE,add=TRUE)
}


# functional to calculate the threshold-weighted CRPS for a normal distribution
twcrps_norm <- function(y, mean, sd, t){
  y_star <- pmax(y, t)
  Ft <- pnorm(t, mean = mean, sd = sd)
  Fy <- pnorm(y, mean = mean, sd = sd)
  Fys <- pnorm(y_star, mean = mean, sd = sd)
  Ft_sq2 <- pnorm(t, mean = mean, sd = sd/sqrt(2))
  dt <- dnorm((t - mean)/sd)
  dys <- dnorm((y_star - mean)/sd)
  
  o1 <- y_star*(2*Fy - 1) - mean*(2*Fys - 1) 
  o2 <- 2*sd*dys - sd/sqrt(pi)
  o3 <- 2*t*(Fys - Fy) - (t - mean)*(Ft^2)
  o4 <- sd*Ft_sq2/sqrt(pi) - 2*sd*Ft*dt
  
  obj <- o1 + o2 + o3 + o4
  return(obj)
}


# functional to calculate the outcome-weighted CRPS for a normal distribution
owcrps_norm <- function(y, mean, sd, t, BS = T){
  y_star <- pmax(y, t)
  
  w_y <- as.numeric(y > t)
  
  obj <- crps_tnorm(y, mean, sd, lower = t)*w_y
  
  if(BS){
    p <- 1 - pnorm(t, mean, sd)
    
    obj <- obj + (p - w_y)^2
  }
  
  return(obj)
}


# functional to calculate the vertically re-scaled CRPS for a normal distribution
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
    obj2 <- mean*(1 - Ft - 2*F0) + (sd^2)*(2*d0 - dt)
  }else{
    obj2 <- mean*(1 - Ft) + (sd^2)*dt
  }
  obj2 <- (obj2 - abs(y)*w_y)*((1 - Ft) - w_y)
  
  obj <- obj1 + obj2
  
  return(obj)
}



################################################################################
# Figure 1

x <- seq(-3.5, 3.5, 0.01)
t <- 0.5
o <- 2

ylims <- c(0, 0.8)

mu <- 0; sigma <- 1

## CRPS
d <- dnorm(x, mu, sigma)

df <- data.frame(x = x, d = d)
ggplot(df) + geom_area(aes(x = x, y = d), fill = "lightblue") + 
  geom_vline(aes(xintercept = o), col = "red") +
  scale_x_continuous(expand = c(0, 0), breaks = c(t, o), labels = c("t", "y")) +
  scale_y_continuous(limits = ylims, expand = c(0, 0), name = "f(x)") +
  theme_bw() +
  theme(legend.title = element_blank(), panel.grid = element_blank()) +
  ggtitle("CRPS")


## twCRPS
d0 <- pnorm(t, mu, sigma)

df <- data.frame(x = x[x >= t], d = d[x >= t])
ggplot() + geom_area(data = df, aes(x = x, y = d), fill = "lightblue") +
  geom_line(aes(x = c(t, t), y = c(0, d0)), col = "lightblue", lwd = 1.5) +
  geom_vline(aes(xintercept = o), col = "red") +
  scale_x_continuous(limits = range(x), expand = c(0, 0), breaks = c(t, o), labels = c("t", "y")) +
  scale_y_continuous(limits = ylims, expand = c(0, 0), name = "f(x)") +
  theme_bw() +
  theme(legend.title = element_blank(), panel.grid = element_blank()) +
  ggtitle("twCRPS")


## owCRPS
c <- 1 - pnorm(t, mu, sigma)
d_trunc <- d/c

df <- data.frame(x = x[x >= t], d = d_trunc[x >= t])
ggplot() + geom_area(data = df, aes(x = x, y = d), fill = "lightblue") +
  geom_vline(aes(xintercept = o), col = "red") + 
  scale_x_continuous(limits = range(x), expand = c(0, 0), breaks = c(t, o), labels = c("t", "y")) +
  scale_y_continuous(limits = ylims/c, expand = c(0, 0), name = "f(x)") +
  theme_bw() +
  theme(legend.title = element_blank(), panel.grid = element_blank()) +
  ggtitle("owCRPS")



################################################################################
# Figure 2

y <- seq(-2, 2, 0.01)
mu <- 0; sigma <- 1

t <- 0.5

s <- crps_norm(y, mu, sigma)
tws <- twcrps_norm(y, mu, sigma, t)
ows <- owcrps_norm(y, mu, sigma, t, BS = F)
ows_bs <- owcrps_norm(y, mu, sigma, t)
vrs <- vrcrps_norm(y, mu, sigma, t)

t_ind <- which(y == t)
lt_t <- 1:t_ind
gt_t <- (t_ind + 1):length(y)


df_lt <- data.frame(s = c(s[lt_t], tws[lt_t], ows[lt_t], ows_bs[lt_t], vrs[lt_t]), y = y[lt_t], 
                    mth = rep(c("CRPS", "twCRPS", "owCRPS", "owCRPS + BS", "vrCRPS"), each = length(y[lt_t])))
df_gt <- data.frame(s = c(s[gt_t], tws[gt_t], ows[gt_t], ows_bs[gt_t], vrs[gt_t]), y = y[gt_t], 
                    mth = rep(c("CRPS", "twCRPS", "owCRPS", "owCRPS + BS", "vrCRPS"), each = length(y[gt_t])))

ggplot() + geom_line(data = df_lt, aes(x = y, y = s, col = mth, linetype = mth)) +
  geom_line(data = df_gt, aes(x = y, y = s, col = mth, linetype = mth)) +
  scale_linetype_manual(values = c("dotted", rep("solid", 4))) +
  scale_color_manual(values = c("black", "red", "blue", "black", "green4")) +
  geom_vline(aes(xintercept = t), col = "grey") +
  scale_x_continuous(breaks = c(-2:2, t), labels = c(-2:2, "t")) +
  scale_y_continuous(name = "S(F, y)") +
  theme_bw() +
  theme(legend.title = element_blank(), panel.grid = element_blank(),
        legend.justification = c(0, 1), legend.position = c(0.3, 0.99)) +
  guides(linetype = guide_legend(nrow = 5))



################################################################################
# Figure 3

n <- 100000 # sample size
bins <- 20 # bins in PIT histogram
t <- 1 # threshold in weight function 

sig2 <- 1/3
sig <- sqrt(sig2)
mu <- rnorm(n, 0, sqrt(1 - sig2))

y <- rnorm(n, mu, sig) # observations
w <- (y > t) # weight function


## compute PIT values
F_y <- pnorm(y, mu, sig) # PIT values
F_y_restricted <- pnorm(y[w], mu[w], sig) # restricted PIT values
F_y_conditional <- (pnorm(y[w], mu[w], sig) - pnorm(t, mu[w], sig))/
  (1 - pnorm(t, mu[w], sig)) # conditional PIT values


## plot PIT values

PIThist(F_y, bins, title = "PIT")
PIThist(F_y_restricted, bins, title = "Restricted PIT")
PIThist(F_y_conditional, bins, title = "Conditional PIT")



################################################################################
# Figure 4

n <- 1000000 # sample size
bins <- 20 # bins in PIT histogram
t <- 2 # threshold in weight function 

log_c <- sqrt(3)/pi # variance normalisation constant

sig2 <- 1/3
sig <- sqrt(sig2)
mu <- rnorm(n, 0, sqrt(1 - sig2))

y <- rlogis(n, mu, sig*sqrt(3)/pi) # observations
w <- (y > t) # weight function


## PIT values

F_x1 <- plogis(y, mu, sig*log_c) 
F_x2 <- pt((y - mu)/sig, df = 5) 
F_x3 <- pnorm(y, mu, sig)


## PIT histograms

PIThist(F_x1, bins, title = "Logistic")
PIThist(F_x2, bins, title = "Student's t")
PIThist(F_x3, bins, title = "Normal")


## PIT reliability diagrams

samp_ind <- sample(1:n, 10000) # take random sample of data (for speed)
PITreldiag(F_x1[samp_ind], y[samp_ind], region_level = 0.99)
PITreldiag(F_x2[samp_ind], y[samp_ind], region_level = 0.99)
PITreldiag(F_x3[samp_ind], y[samp_ind], region_level = 0.99)


## cPIT values

cF_x1 <- (plogis(y[w], mu[w], sig*log_c) - plogis(t, mu[w], sig*log_c))/
  (1 - plogis(t, mu[w], sig*log_c))
cF_x2 <- (pt((y[w] - mu[w])/sig, df = 5) - pt((t - mu[w])/sig, df = 5))/
  (1 - pt((t - mu[w])/sig, df = 5))
cF_x3 <- (pnorm(y[w], mu[w], sig) - pnorm(t, mu[w], sig))/
  (1 - pnorm(t, mu[w], sig))


## cPIT histograms

PIThist(cF_x1, bins, title = "Logistic")
PIThist(cF_x2, bins, title = "Student's t")
PIThist(cF_x3, bins, title = "Normal")


## cPIT reliability diagrams

PITreldiag(cF_x1, y[w], region_level = 0.99)
PITreldiag(cF_x2, y[w], region_level = 0.99)
PITreldiag(cF_x3, y[w], region_level = 0.99)


## standard reliability diagrams

# probability of exceedance
pF_x1 <- 1 - plogis(t, mu, sig*sqrt(3)/pi)
pF_x2 <- 1 - pt((t - mu)/sig, df = 5)
pF_x3 <- 1 - pnorm(t, mu, sig)

reliabilitydiag(X = pF_x1, y = as.numeric(w), xtype = "continuous")
reliabilitydiag(X = pF_x2, y = as.numeric(w), xtype = "continuous")
reliabilitydiag(X = pF_x3, y = as.numeric(w), xtype = "continuous")


