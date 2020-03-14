# Set working directory and load libraries
setwd("C:/Users/Simon/Google Drev/Uni/Cognitive modelling")
set.seed(1995)
library(pacman)
p_load(R2jags, extraDistr, magrittr, dplyr, ggplot2)
library(cvms)

#--------------------Task context -------------------------

# Generate a payoff matrix for bandit task
# Choice of bandit a = 30 % chance of payoff 2 otherwise 0
# Choice of bandit b = 70 % chance of payoff 1 otherwise 0

ntrials = 100

# Decide on probabilities and reward
Aprob <- 0.3
Arew <- 2
Bprob <- 0.7
Brew <- 1

# Make a payoff matrix
payoff <- cbind(rbinom(ntrials, 1, Aprob)*Arew, rbinom(ntrials, 1, Bprob)*Brew) 

# Check which bandit is better 
colSums(payoff)

# ------------------ Q-learning --------------------------

# Define function
qlearn <- function(payoff, ntrials, alfa, beta){
    
  x <- array(0, c(ntrials)) #Choice
  r <- array(0, c(ntrials)) #Reward
  Q <- array(0, c(ntrials, 2)) #"Expected Value"
  Qupdate <- array(0,c(ntrials,2))
  exp_p <- array(0, c(ntrials,2)) # Exponential prob
  p <- array(0, c(ntrials,2)) # Probability
  
  # Trial 1
  Q[1,1] <- 1
  Q[1,2] <- 1
  
  x[1] <- rcat(1, c(1/2, 1/2))
  r[1] <- payoff[1, x[1]]
  
  # Trial 2 and onwards
  for (t in 2:ntrials){
    
    for (k in 1:2){
      #Only update chosen option
      Qupdate[t,k] <- Q[t-1,k] + alfa*(r[t-1]-Q[t-1,k])
      Q[t,k] <- ifelse(k==x[t-1], Qupdate[t,k], Q[t-1,k])
      
      exp_p[t,k] <- exp(beta*Q[t,k])
    }
    
    for (k in 1:2){
      p[t,k] <- exp_p[t,k]/sum(exp_p[t,])
    }
    
    x[t] <- rcat(1, p[t,])
    r[t] <- payoff[t,x[t]]
  }
  
  result <- list(x=x, r=r, Q=Q)
  return(result)
}

# Set values
ntrials = 100
alfa = .4
beta = 3
  
# Run function
q.sims=qlearn(payoff, ntrials, alfa, beta)

# Plotting parameter
par(mfrow=c(3,1))

plot(q.sims$Q[,1])
plot(q.sims$Q[,2])
plot(q.sims$x)

#Jags stuff
x <- q.sims$x
r <- q.sims$r

data <- list("x", "r", "ntrials")

params <- c("alfa", "beta")

samples <- jags.parallel(data, inits = NULL, params, 
                model.file = "qlearn.txt",
                n.chains = 3, n.iter = 5000, n.burnin = 1000, n.thin = 1)

traceplot(samples)
par(mfrow=c(1,1))
plot(density(samples$BUGSoutput$sims.list$alfa))
plot(density(samples$BUGSoutput$sims.list$beta))


# Parameter recovery
niter = 100
true_alfa <- array(0, c(niter))
true_beta <- array(0, c(niter))

infer_alfa <- array(0, c(niter))
infer_beta <- array(0, c(niter))

for (i in 1:niter){
  
  # true parameters
  alfa <- runif(1, 0, 1)
  beta <- runif(1, 0, 5)
  
  # Run function and extract responses
  q.sims <- qlearn(payoff, ntrials, alfa, beta)
  x <- q.sims$x
  r <- q.sims$r
  
  data <- list("x", "r", "ntrials")
  
  params <- c("alfa", "beta")
  
  samples <- jags.parallel(data, inits = NULL, params, 
                           model.file = "qlearn.txt",
                           n.chains = 3, n.iter = 5000, n.burnin = 1000, n.thin = 1)
  
  true_alfa[i] <- alfa
  true_beta[i] <- beta
  
  
  alfa_pos = samples$BUGSoutput$sims.list$alfa
  beta_pos = samples$BUGSoutput$sims.list$beta
  
  infer_alfa[i] <- density(alfa_pos)$x[which(density(alfa_pos)$y==max(density(alfa_pos)$y))]
  infer_beta[i] <- density(beta_pos)$x[which(density(beta_pos)$y==max(density(beta_pos)$y))]
  
}

plot(true_alfa, infer_alfa)
abline(0,1)
title(main = "Parameter recovery for alpha")
plot(true_beta, infer_beta) 
abline(0,1)
title(main = "Parameter recovery for beta")
# -------------------Choice kernal ------------------------

# Create function
choice_kernal = function(payoff, ntrials, a, beta){
  x <- array(0, c(ntrials)) #Choice
  r <- array(0, c(ntrials)) #Reward
  CK <- array(0, c(ntrials, 2)) 
  CKchosen <- array(0,c(ntrials,2))
  CKunchosen <- array(0,c(ntrials,2))
  exp_p <- array(0, c(ntrials,2)) # Exponential prob
  p <- array(0, c(ntrials,2)) # Probability
  
  # Trial 1
  CK[1,1] <- 1
  CK[1,2] <- 1
  
  x[1] <- rcat(1, c(1/2, 1/2))
  r[1] <- payoff[1, x[1]]
  
  # Trial 2 and onwards
  for (t in 2:ntrials){
    
    for (k in 1:2){
      #Only update chosen option
      CKchosen[t,k] <- CK[t-1,k]+a*(1-CK[t-1,k])
      CKunchosen[t,k] <- CK[t-1,k]+a*(0-CK[t-1,k])
      
      CK[t,k] <- ifelse(k==x[t-1], CKchosen[t,k], CKunchosen[t,k])
      
      exp_p[t,k] <- exp(beta*CK[t,k])
    }
    
    for (k in 1:2){
      p[t,k] <- exp_p[t,k]/sum(exp_p[t,])
    }
    
    x[t] <- rcat(1, p[t,])
    r[t] <- payoff[t,x[t]]
  
  }
  result <- list(x=x, r=r, CK=CK)
  return(result)
}  

# Set values
ntrials = 100
a = .1
beta = 5

# run forward simulation and do some visuals 
choice.sims=choice_kernal(payoff, ntrials, a, beta)
plot(choice.sims$x)
plot(choice.sims$r)
plot(choice.sims$CK[,2])

# Inference - Jags stuff
x <- choice.sims$x
r <- choice.sims$r

data <- list("x", "r", "ntrials")

params <- c("alfa", "beta")

samples <- jags.parallel(data, inits = NULL, params, 
                         model.file = "kernal.txt",
                         n.chains = 3, n.iter = 5000, n.burnin = 1000, n.thin = 1)

traceplot(samples)
par(mfrow=c(1,1))
plot(density(samples$BUGSoutput$sims.list$alfa))
plot(density(samples$BUGSoutput$sims.list$beta))


# Parameter recovery for the choice kernal model
niter = 100
true_alfa <- array(0, c(niter))
true_beta <- array(0, c(niter))

infer_alfa <- array(0, c(niter))
infer_beta <- array(0, c(niter))

for (i in 1:niter){
  
  # true parameters
  alfa <- runif(1, 0, 1)
  beta <- runif(1, 0, 5)
  
  # Run function and extract responses
  choice.sims=choice_kernal(payoff, ntrials, alfa, beta)
  x <- choice.sims$x
  r <- choice.sims$r
  
  data <- list("x", "r", "ntrials")
  
  params <- c("alfa", "beta")
  
  samples <- jags.parallel(data, inits = NULL, params, 
                           model.file = "kernal.txt",
                           n.chains = 3, n.iter = 5000, n.burnin = 1000, n.thin = 1)
  
  true_alfa[i] <- alfa
  true_beta[i] <- beta
  
  
  alfa_pos = samples$BUGSoutput$sims.list$alfa
  beta_pos = samples$BUGSoutput$sims.list$beta
  
  infer_alfa[i] <- density(alfa_pos)$x[which(density(alfa_pos)$y==max(density(alfa_pos)$y))]
  infer_beta[i] <- density(beta_pos)$x[which(density(beta_pos)$y==max(density(beta_pos)$y))]
  
  print(i)
}

plot(true_alfa, infer_alfa)
abline(0,1)
title(main = "Parameter recovery for alpha")
plot(true_beta, infer_beta) 
abline(0,1)
title(main = "Parameter recovery for beta")

# ------------------Model recovery -----------------------
# Number of iterations
niter = 100
# Set up empty arrays
DICs_RW_dat <- array(0, c(niter, 2))
DICs_kernal_dat <- array(0, c(niter,2))

for (i in 1:niter){
  
  alfa <- runif(1,0,1)
  beta <- rgamma(1, 1, 1)
  
  # Run both simulations
  RW.sims = qlearn(payoff, 100, alfa, beta)
  
  kernal.sims = choice_kernal(payoff, 100, alfa, beta)
  
  # RUN RW inference on RWdata
  x <- RW.sims$x
  r <- RW.sims$r
  
  data <- list("x", "r", "ntrials")
  
  params <- c("alfa", "beta")
  
  RWdat_RWmodel <- jags.parallel(data, inits = NULL, params, 
                           model.file = "qlearn.txt",
                           n.chains = 3, n.iter = 5000, n.burnin = 1000, n.thin = 1)
  
  # RUN RW inference on kernaldata
  x <- kernal.sims$x
  r <- kernal.sims$r
  
  data <- list("x", "r", "ntrials")
  
  params <- c("alfa", "beta")
  
  Kernaldat_RWmodel <- jags.parallel(data, inits = NULL, params, 
                           model.file = "qlearn.txt",
                           n.chains = 3, n.iter = 5000, n.burnin = 1000, n.thin = 1)
  
  # RUN kernal inference on kernaldata
  x <- kernal.sims$x
  r <- kernal.sims$r
  
  data <- list("x", "r", "ntrials")
  
  params <- c("alfa", "beta")
  
  Kernaldat_kernalmodel <- jags.parallel(data, inits = NULL, params, 
                                         model.file = "kernal.txt",
                                         n.chains = 3, n.iter = 5000, n.burnin = 1000, n.thin = 1)
  
  # RUN kernal inference on RW data
  x <- RW.sims$x
  r <- RW.sims$r
  
  data <- list("x", "r", "ntrials")
  
  params <- c("alfa", "beta")
  
  RWdat_kernalmodel <- jags.parallel(data, inits = NULL, params, 
                                         model.file = "kernal.txt",
                                         n.chains = 3, n.iter = 5000, n.burnin = 1000, n.thin = 1)
  
  # Fill up DIC arrays
  DICs_RW_dat[i,1] <- RWdat_RWmodel$BUGSoutput$DIC
  DICs_RW_dat[i,2] <- RWdat_kernalmodel$BUGSoutput$DIC
  DICs_kernal_dat[i,1] <- Kernaldat_kernalmodel$DIC
  DICs_kernal_dat[i, 2] <- Kernaldat_RWmodel$BUGSoutput$DIC
  
  print(i)
}

best_RW <- array(0, c(niter))
best_kernal <- array(0, c(niter))

# Find minimum
for (i in 1:niter){
  best_RW[i] <- which.min(DICs_RW_dat[i,])
  best_kernal[i] <- which.min((DICs_kernal_dat[i,]))
}

best_RW2=recode(best_RW, '1' = 'Qlearn', '2' = 'Kernal')
best_kernal2 = recode(best_kernal, '1'= 'Kernal', '2' = 'Qlearn')

#Confusion matrix

# Make predictions
n_sims = 100
dic_df <- data.frame(predictions = c(best_RW2, best_kernal2), true = c(rep('Qlearn', n_sims), rep('Kernal', n_sims))) 

conf_mat <-  confusion_matrix(dic_df$true, dic_df$predictions)
cf_mat <- conf_mat$`Confusion Matrix`[[1]]

plot_confusion_matrix(cf_mat, add_row_percentages = F, add_col_percentages = F, add_normalized = F, counts_on_top = F) + labs( x= "Model from which the data is produced") + ylab("Best fitting model")
