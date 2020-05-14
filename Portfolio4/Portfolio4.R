#----------- Data cleaning and inference test---------

# Set up working directory and load in data
pacman::p_load(R2jags, polspline)

setwd("C:/Users/Simon/Google Drev/Uni/Cognitive modelling/IGT_rawdata")

control = read.table('IGTdata_healthy_control.txt', header = TRUE)

# Get list of unique subjects and number of trials and subjects
subIDs <- unique(control$subjID)
nsubs <- length(subIDs)
ntrials_max <- 100

# Extract choices, outcome, rewards and losses 
choice = control$deck
loss = control$loss
gain = control$gain
outcome = gain +loss

# Create empty arrays to fill 
ntrials_all <- array(0,c(nsubs))
choice_all <- array(0, c(nsubs, ntrials_max))
loss_all <- array(0, c(nsubs, ntrials_max))
gain_all <- array(0, c(nsubs, ntrials_max))
outcome_all <- array(0, c(nsubs, ntrials_max))

# Create a loop that populates these arrays
for (s in 1:nsubs){
  
  ntrials_all[s] <- length(choice[control$subjID==subIDs[s]])
  
  choice_sub <- choice[control$subjID==subIDs[s]]
  length(choice_sub) = ntrials_max
  
  loss_sub <- loss[control$subjID==subIDs[s]]
  length(loss_sub) = ntrials_max
  
  gain_sub <- gain[control$subjID==subIDs[s]]
  length(gain_sub) = ntrials_max
  
  outcome_sub <- outcome[control$subjID==subIDs[s]]
  length(outcome_sub) = ntrials_max
  
  # Populate arrays
  
  choice_all[s,] = choice_sub
  loss_all[s,]= loss_sub
  gain_all[s,] = gain_sub
  outcome_all[s,] = outcome_sub
  
}

# Scale
loss_scale = loss_all/100
gain_scale = gain_all/100
outcome_scale = outcome_all/100

# Apply PVL-model to one subject to check that it works
setwd("C:/Users/Simon/Google Drev/Uni/Cognitive modelling")

# Extract data from one participant
x=choice_all[1,]
r = outcome_all[1,]
ntrials = ntrials_all[1]

# Run JAGS
data <- list('x', 'r', 'ntrials')
params <- c('w', 'A', 'a', 'theta', 'p')
samples <- jags.parallel(data, inits = NULL, params, 
                         model.file = "PVLD.txt",
                         n.chains = 3, n.iter = 5000, n.burnin = 1000, n.thin = 1)

# Check the traceplot & convergence diagnostics
traceplot(samples)
sammples

# Plot choice probability on trial 32 
p_post = samples$BUGSoutput$sims.list$p

par(mfrow = c(2,2))
plot(density(p_post[,32,1]))
plot(density(p_post[,32,2])) # Deck 2 has the highest probability of being chosen
plot(density(p_post[,32,3]))
plot(density(p_post[,32,4]))

# ------ Posterior predictive checks-------------

#PVL-Delta model

# Loop that predicts choices for all subjects on all trials
pred_succes <- array(c(nsubs))

for (s in 1:nsubs){
  
  x=choice_all[s,]
  r=outcome_scale[s,]
  
  ntrials <- ntrials_all[s]
  
  #JAGS
  data <- list('x', 'r', 'ntrials')
  params <- c('w', 'A', 'a', 'theta', 'p')
  samples <- jags.parallel(data, inits = NULL, params, 
                           model.file = "PVLD.txt",
                           n.chains = 3, n.iter = 5000, n.burnin = 1000, n.thin = 1)
  
  p_post <- samples$BUGSoutput$sims.list$p
  
  
  x_predict <- array(c(ntrials))
  for (t in 1:ntrials){
    
    p_predict <- c(
      density(p_post[,t,1])$x[which(density(p_post[,t,1])$y==max(density(p_post[,t,1])$y))],
      density(p_post[,t,1])$x[which(density(p_post[,t,2])$y==max(density(p_post[,t,2])$y))],
      density(p_post[,t,1])$x[which(density(p_post[,t,3])$y==max(density(p_post[,t,3])$y))],
      density(p_post[,t,1])$x[which(density(p_post[,t,4])$y==max(density(p_post[,t,4])$y))])
    
    x_predict[t] <- which.max(p_predict)
  }
  
  pred_succes[s] <- sum(x_predict == x, na.rm = TRUE)
  print (s)
}

PP_PVL = pred_succes
mean(PP_PVL) # Mean 27.5 % (unscaled) and 25.35 % (scaled)

# ORL model

# Loop that predicts choices for all subjects on all trials
pred_succes <- array(c(nsubs))

for (s in 1:nsubs){
  
  x=choice_all[s,]
  r=outcome_scale[s,]
  
  ntrials <- ntrials_all[s]
  
  #JAGS
  data <- list('x', 'r', 'ntrials')
  params <- c("Arew", "Apun", "K", "wp", "wf", "p")
  samples <- jags.parallel(data, inits = NULL, params, 
                           model.file = "ORL.txt",
                           n.chains = 3, n.iter = 5000, n.burnin = 1000, n.thin = 1)
  
  p_post <- samples$BUGSoutput$sims.list$p
  
  
  x_predict <- array(c(ntrials))
  for (t in 1:ntrials){
    
    p_predict <- c(
      density(p_post[,t,1])$x[which(density(p_post[,t,1])$y==max(density(p_post[,t,1])$y))],
      density(p_post[,t,1])$x[which(density(p_post[,t,2])$y==max(density(p_post[,t,2])$y))],
      density(p_post[,t,1])$x[which(density(p_post[,t,3])$y==max(density(p_post[,t,3])$y))],
      density(p_post[,t,1])$x[which(density(p_post[,t,4])$y==max(density(p_post[,t,4])$y))])
    
    x_predict[t] <- which.max(p_predict)
  }
  
  pred_succes[s] <- sum(x_predict == x, na.rm = TRUE)
  print (s)
}

PP_ORL = pred_succes
mean(PP_ORL) # 43.625 % (scaled)

# VSE model model

# Loop that predicts choices for all subjects on all trials
pred_succes <- array(c(nsubs))

for (s in 1:nsubs){
  
  x=choice_all[s,]
  r=gain_scale[s,]
  l=loss_scale[s,]
  l = abs(l)
  
  ntrials <- ntrials_all[s]
  
  #JAGS
  data <- list('x', 'r', 'l', 'ntrials')
  params <- c("theta", "decay" , "alfa", "phi", "C", "p")
  samples <- jags.parallel(data, inits = NULL, params, 
                           model.file = "VSE.txt",
                           n.chains = 3, n.iter = 5000, n.burnin = 1000, n.thin = 1)
  
  p_post <- samples$BUGSoutput$sims.list$p
  
  
  x_predict <- array(c(ntrials))
  for (t in 1:ntrials){
    
    p_predict <- c(
      density(p_post[,t,1])$x[which(density(p_post[,t,1])$y==max(density(p_post[,t,1])$y))],
      density(p_post[,t,1])$x[which(density(p_post[,t,2])$y==max(density(p_post[,t,2])$y))],
      density(p_post[,t,1])$x[which(density(p_post[,t,3])$y==max(density(p_post[,t,3])$y))],
      density(p_post[,t,1])$x[which(density(p_post[,t,4])$y==max(density(p_post[,t,4])$y))])
    
    x_predict[t] <- which.max(p_predict)
  }
  
  pred_succes[s] <- sum(x_predict == x, na.rm = TRUE)
  print (s)
}

PP_VSE = pred_succes
mean(PP_VSE) # Mean accuracy is 34 % (scaled) 

# ------ Hierachical modeling --------------

# PVL-Delta
x <- choice_all
r <- outcome_all

data <- list('x', 'r', 'ntrials_all', 'nsubs')
params <- c("mu_w", "mu_A", "mu_theta", "mu_a", "l_w", "l_A", "l_theta", "l_a")
PVL_samples <- jags.parallel(data, inits = NULL, params, 
                         model.file = "Hier_PVL.txt",
                         n.chains = 3, n.iter = 5000, n.burnin = 1000, n.thin = 1)

# Check convergence
PVL_samples

par(mfrow = c(2,2))
plot(density(PVL_samples$BUGSoutput$sims.list$mu_A), main = "mu_A") 
plot(density(PVL_samples$BUGSoutput$sims.list$mu_w), main = "mu_w")
plot(density(PVL_samples$BUGSoutput$sims.list$mu_a), main = "mu_a")
plot(density(PVL_samples$BUGSoutput$sims.list$mu_theta), main = "mu_theta")

#ORL model

x <- choice_all
r <- outcome_scale

data <- list('x', 'r', 'ntrials_all', 'nsubs')
params <- c("mu_Arew", "mu_Apun", "mu_K", "mu_wp", "mu_wf", "l_Arew", "l_Apun", "l_K", "l_wp", "l_wf")
ORL_samples <- jags.parallel(data, inits = NULL, params, 
                         model.file = "Hier_ORL_NEW.txt",
                         n.chains = 3, n.iter = 5000, n.burnin = 1000, n.thin = 1)

# Check convergence diagnostics
ORL_samples
# Plot group posteriors
par(mfrow = c(3,2))
plot(density(ORL_samples$BUGSoutput$sims.list$mu_Arew), main = "mu Arew")
plot(density(ORL_samples$BUGSoutput$sims.list$mu_Apun), main = "mu Apun")
plot(density(ORL_samples$BUGSoutput$sims.list$mu_K), main = "mu K")
plot(density(ORL_samples$BUGSoutput$sims.list$mu_wp), main = "mu wp")
plot(density(ORL_samples$BUGSoutput$sims.list$mu_wf), main = "mu wf")

# VSE model

x <- choice_all
r <- gain_scale
l <- abs(loss_scale)

data <- list('x', 'r', 'l', 'ntrials_all', 'nsubs')
params <- c("mu_theta", "mu_decay", "mu_alfa", "mu_phi", "mu_C", "l_theta", "l_decay", "l_alfa", "l_phi", "l_C")
VSE_samples <- jags.parallel(data, inits = NULL, params, 
                         model.file = "Hier_VSE.txt",
                         n.chains = 3, n.iter = 5000, n.burnin = 1000, n.thin = 1)
# Check convergence
VSE_samples

par(mfrow = c(3,2))
plot(density(VSE_samples$BUGSoutput$sims.list$mu_theta), main = "Theta")
plot(density(VSE_samples$BUGSoutput$sims.list$mu_decay), main = "Decay")
plot(density(VSE_samples$BUGSoutput$sims.list$mu_alfa), main = "Alfa")
plot(density(VSE_samples$BUGSoutput$sims.list$mu_phi), main = "Phi")
plot(density(VSE_samples$BUGSoutput$sims.list$mu_C), main = "C")

 #--------------Group comparison on trail-by-trial data -----------------

# Load in all data

# Control data
X_ctr = outcome_scale
ntrials_ctr = ntrials_all
nsubs_ctr = nsubs

# Opiate group
setwd("C:/Users/Simon/Google Drev/Uni/Cognitive modelling/IGT_rawdata")
Opi = read.table('IGTdata_heroin.txt', header = TRUE)

subIDs <- unique(Opi$subjID)
nsubs_opi <- length(subIDs)
ntrials_max <- 100

# Extract choices, outcome, rewards and losses 
choice = Opi$deck
loss = Opi$loss
gain = Opi$gain
outcome = gain +loss

# Create empty arrays to fill 
ntrials_opi <- array(0,c(nsubs_opi))
choice_opi <- array(0, c(nsubs_opi, ntrials_max))
loss_opi <- array(0, c(nsubs_opi, ntrials_max))
gain_opi <- array(0, c(nsubs_opi, ntrials_max))
outcome_opi <- array(0, c(nsubs_opi, ntrials_max))

# Create a loop that populates these arrays
for (s in 1:nsubs_opi){
  
  ntrials_opi[s] <- length(choice[Opi$subjID==subIDs[s]])
  
  choice_sub <- choice[Opi$subjID==subIDs[s]]
  length(choice_sub) = ntrials_max
  
  loss_sub <- loss[Opi$subjID==subIDs[s]]
  length(loss_sub) = ntrials_max
  
  gain_sub <- gain[Opi$subjID==subIDs[s]]
  length(gain_sub) = ntrials_max
  
  outcome_sub <- outcome[Opi$subjID==subIDs[s]]
  length(outcome_sub) = ntrials_max
  
  # Populate arrays
  
  choice_opi[s,] = choice_sub
  loss_opi[s,]= loss_sub
  gain_opi[s,] = gain_sub
  outcome_opi[s,] = outcome_sub
  
}

# Scale
loss_scale_opi = loss_opi/100
gain_scale_opi = gain_opi/100
outcome_scale_opi = outcome_opi/100

# Amphetamine group
setwd("C:/Users/Simon/Google Drev/Uni/Cognitive modelling/IGT_rawdata")
Amp = read.table('IGTdata_amphetamine.txt', header = TRUE)

subIDs <- unique(Amp$subjID)
nsubs_amp <- length(subIDs)
ntrials_max <- 100

# Extract choices, outcome, rewards and losses 
choice = Amp$deck
loss = Amp$loss
gain = Amp$gain
outcome = gain +loss

# Create empty arrays to fill 
ntrials_amp <- array(0,c(nsubs_amp))
choice_amp <- array(0, c(nsubs_amp, ntrials_max))
loss_amp <- array(0, c(nsubs_amp, ntrials_max))
gain_amp <- array(0, c(nsubs_amp, ntrials_max))
outcome_amp <- array(0, c(nsubs_amp, ntrials_max))

# Create a loop that populates these arrays
for (s in 1:nsubs_amp){
  
  ntrials_amp[s] <- length(choice[Amp$subjID==subIDs[s]])
  
  choice_sub <- choice[Amp$subjID==subIDs[s]]
  length(choice_sub) = ntrials_max
  
  loss_sub <- loss[Amp$subjID==subIDs[s]]
  length(loss_sub) = ntrials_max
  
  gain_sub <- gain[Amp$subjID==subIDs[s]]
  length(gain_sub) = ntrials_max
  
  outcome_sub <- outcome[Amp$subjID==subIDs[s]]
  length(outcome_sub) = ntrials_max
  
  # Populate arrays
  
  choice_amp[s,] = choice_sub
  loss_amp[s,]= loss_sub
  gain_amp[s,] = gain_sub
  outcome_amp[s,] = outcome_sub
  
}

# Scale
loss_scale_amp = loss_amp/100
gain_scale_amp = gain_amp/100
outcome_scale_amp = outcome_amp/100

# Opiate group vs. controls
setwd("C:/Users/Simon/Google Drev/Uni/Cognitive modelling")
X_add <- outcome_scale_opi
ntrials_add <- ntrials_opi
nsubs_add <- nsubs_opi

data <- list("X_ctr", 'ntrials_ctr', 'nsubs_ctr',
             'X_add', 'ntrials_add', 'nsubs_add')
params <- c('mu', 'alpha', 'Smu_ctr', 'Smu_add')
Opi_ctr <- jags.parallel(data, inits = NULL, params, 
                         model.file = "GroupCompare.txt",
                         n.chains = 3, n.iter = 5000, n.burnin = 1000, n.thin = 1)

#traceplot(Opi_ctr)

# Plot density of mu for comparison
plot(density(Opi_ctr$BUGSoutput$median$Smu_add), col = 'blue', main = "Control (red) vs. Opiate addicts (blue)", xlim = c(-0.2, 0.1))
lines(density(Opi_ctr$BUGSoutput$median$Smu_ctr), col = 'red')

# Alpha plot
plot(density(Opi_ctr$BUGSoutput$sims.list$alpha), main = "Posterior for alpha (Opi vs. ctr)")

# Bayes factor 
fit_pos <- logspline(Opi_ctr$BUGSoutput$sims.list$alpha)
null_Opi_pos <- dlogspline(0, fit_pos)

null_prior <- dnorm(0, 0, (1/sqrt(.1)))

opi_BF <- null_Opi_pos/null_prior
opi_BF # 29.14

# Amphetamine group vs. controls

setwd("C:/Users/Simon/Google Drev/Uni/Cognitive modelling")
X_add <- outcome_scale_amp
ntrials_add <- ntrials_amp
nsubs_add <- nsubs_amp

data <- list("X_ctr", 'ntrials_ctr', 'nsubs_ctr',
             'X_add', 'ntrials_add', 'nsubs_add')
params <- c('mu', 'alpha', 'Smu_ctr', 'Smu_add')
Amp_ctr <- jags.parallel(data, inits = NULL, params, 
                         model.file = "GroupCompare.txt",
                         n.chains = 3, n.iter = 5000, n.burnin = 1000, n.thin = 1)
# Plot density of mu for comparison
plot(density(Amp_ctr$BUGSoutput$median$Smu_add), col = 'blue', main = "Control (red) vs. Amphetamine addicts (blue)", xlim = c(-0.2, 0.1))
lines(density(Amp_ctr$BUGSoutput$median$Smu_ctr), col = 'red')

# Alpha plot
plot(density(Amp_ctr$BUGSoutput$sims.list$alpha), main = "Posterior for alpha (Amp vs. ctr)")

# Bayes factor 
fit_pos <- logspline(Amp_ctr$BUGSoutput$sims.list$alpha)
null_amp_pos <- dlogspline(0, fit_pos)

null_prior <- dnorm(0, 0, (1/sqrt(.1)))

amp_BF <- null_amp_pos/null_prior
amp_BF # 43.51

#-------------Hierarchical Group Comparison Using the ORL model --------------- 
# Opiate group

r_ctr = outcome_scale
r_add = outcome_scale_opi
x_ctr = choice_all
x_add = choice_opi
ntrials_ctr = ntrials_all
nsubs_ctr = nsubs
ntrials_add = ntrials_opi
nsubs_add = nsubs_opi

data <- list("x_ctr", 'r_ctr', 'ntrials_ctr', 'nsubs_ctr',
             'x_add', 'r_add', 'ntrials_add', 'nsubs_add')
params <- c("alpha_Arew", 'alpha_Apun', 'alpha_K', 'alpha_wp', "alpha_wf",
            "Arew_add", "Apun_add", "K_add", "wp_add", "wf_add",
            "Arew_ctr", "Apun_ctr", "K_ctr", "wp_ctr", "wf_ctr")
ORL_opiVSctr<- jags(data, inits = NULL, params, 
                         model.file = "ORLCompare.txt",
                         n.chains = 3, n.iter = 5000, n.burnin = 1000, n.thin = 1)

# Plot alpha
par(mfrow = c(3,2))
plot(density(ORL_opiVSctr$BUGSoutput$sims.list$alpha_Apun), main = "Apun")
lines(density((rnorm(10000, 0, 1))), col = "purple")
plot(density(ORL_opiVSctr$BUGSoutput$sims.list$alpha_Arew), main = "Arew")
lines(density((rnorm(10000, 0, 1))), col = "purple")
plot(density(ORL_opiVSctr$BUGSoutput$sims.list$alpha_K), main = "K")
lines(density((rnorm(10000, 0, 1))), col = "purple")
plot(density(ORL_opiVSctr$BUGSoutput$sims.list$alpha_wp), main = "wp")
lines(density((rnorm(10000, 0, 1))), col = "purple")
plot(density(ORL_opiVSctr$BUGSoutput$sims.list$alpha_wf), main = "wf")
lines(density((rnorm(10000, 0, 1))), col = "purple")

# Plot mu
par(mfrow = c(3,2))
plot(density(ORL_opiVSctr$BUGSoutput$sims.list$Apun_ctr), main = "Apun", ylim = c(0,40))
lines(density(ORL_opiVSctr$BUGSoutput$sims.list$Apun_add), col = 'red')
plot(density(ORL_opiVSctr$BUGSoutput$sims.list$Arew_ctr), main = "Arew")
lines(density(ORL_opiVSctr$BUGSoutput$sims.list$Arew_add), col = 'red')
plot(density(ORL_opiVSctr$BUGSoutput$sims.list$K_ctr), main = "K")
lines(density(ORL_opiVSctr$BUGSoutput$sims.list$K_add), col = 'red')
plot(density(ORL_opiVSctr$BUGSoutput$sims.list$wp_ctr), main = "wp")
lines(density(ORL_opiVSctr$BUGSoutput$sims.list$wp_add), col = 'red')
plot(density(ORL_opiVSctr$BUGSoutput$sims.list$wf_ctr), main = "wf", ylim = c(0,0.5))
lines(density(ORL_opiVSctr$BUGSoutput$sims.list$wf_add), col = 'red')


# Amp group 
r_add = outcome_scale_amp
x_add = choice_amp
ntrials_add = ntrials_amp
nsubs_add = nsubs_amp

data <- list("x_ctr", 'r_ctr', 'ntrials_ctr', 'nsubs_ctr',
             'x_add', 'r_add', 'ntrials_add', 'nsubs_add')
params <- c("alpha_Arew", 'alpha_Apun', 'alpha_K', 'alpha_wp', "alpha_wf",
            "Arew_add", "Apun_add", "K_add", "wp_add", "wf_add",
            "Arew_ctr", "Apun_ctr", "K_ctr", "wp_ctr", "wf_ctr")
ORL_ampVSctr<- jags.parallel(data, inits = NULL, params, 
                             model.file = "ORLCompare.txt",
                             n.chains = 3, n.iter = 5000, n.burnin = 1000, n.thin = 1)

#Plot alpha
par(mfrow = c(3,2))
plot(density(ORL_ampVSctr$BUGSoutput$sims.list$alpha_Apun), main = "Apun")
lines(density((rnorm(10000, 0, 1))), col = "green")
plot(density(ORL_ampVSctr$BUGSoutput$sims.list$alpha_Arew), main = "Arew")
lines(density((rnorm(10000, 0, 1))), col = "green")
plot(density(ORL_ampVSctr$BUGSoutput$sims.list$alpha_K), main = "K")
lines(density((rnorm(10000, 0, 1))), col = "green")
plot(density(ORL_ampVSctr$BUGSoutput$sims.list$alpha_wp), main = "wp")
lines(density((rnorm(10000, 0, 1))), col = "green")
plot(density(ORL_ampVSctr$BUGSoutput$sims.list$alpha_wf), main = "wf")
lines(density((rnorm(10000, 0, 1))), col = "green")


# Plot mu
par(mfrow = c(3,2))
plot(density(ORL_ampVSctr$BUGSoutput$sims.list$Apun_ctr), main = "Apun", ylim = c(0,30))
lines(density(ORL_ampVSctr$BUGSoutput$sims.list$Apun_add), col = 'blue')
plot(density(ORL_ampVSctr$BUGSoutput$sims.list$Arew_ctr), main = "Arew", ylim = c(0,10))
lines(density(ORL_ampVSctr$BUGSoutput$sims.list$Arew_add), col = 'blue')
plot(density(ORL_ampVSctr$BUGSoutput$sims.list$K_ctr), main = "K")
lines(density(ORL_ampVSctr$BUGSoutput$sims.list$K_add), col = 'blue')
plot(density(ORL_ampVSctr$BUGSoutput$sims.list$wp_ctr), main = "wp")
lines(density(ORL_ampVSctr$BUGSoutput$sims.list$wp_add), col = 'blue')
plot(density(ORL_ampVSctr$BUGSoutput$sims.list$wf_ctr), main = "wf", ylim = c(0,0.5))
lines(density(ORL_ampVSctr$BUGSoutput$sims.list$wf_add), col = 'blue')
