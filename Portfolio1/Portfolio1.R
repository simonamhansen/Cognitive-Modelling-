#Load JAGS
library(R2jags)
library(magrittr)
library(dplyr)
library(cvms)

# Working directory
setwd("C:/Users/Simon/Google Drev/Uni/Cognitive modelling")

#-------- Model 1: Fixed theta model--------------
SimfixedModel=function(ntrials, theta){

Gfixed <- array(0, c(ntrials)) # Empty variable to fill out

  for (t in 1:ntrials){
    Gfixed[t] <- rbinom(1, 1, theta)  
  }
  return(Gfixed)
}
#--------- Model 2: Learning model -----------------
SimlearningModel=function(ntrials, theta1, alpha){

Glearn <- array(0, c(ntrials))
theta <- array(0, c(ntrials))

# First trial
theta[1] <- theta1
Glearn[1] <- rbinom(1, 1, theta[1])

# The rest of the trials
for (t in 2:ntrials) {
  theta[t] <- theta[t-1]^(1/(1+alpha)) # Learning function
  Glearn[t] <- rbinom(1, 1, theta[t])
}
return(Glearn)
}
#-------- Run inference for the fixed model using the fixed data--------
InferfixedModel = function(ntrials, G){
data <- list("G", "ntrials")
params <- c("theta")

samples <- jags(data, inits=NULL, params, model.file = "chick.txt", 
                n.chains = 3, n.iter = 5000, n.burnin = 1000, n.thin =1)
# To get ouput
return(samples)
}

#Rhat is convergence diagnostics. Rhat should be close to 1.
#DIC: Deviance Information Criteria

# Traceplot
traceplot(samples) 

# Posterior for theta
plot(density(samples$BUGSoutput$sims.list$theta))


#-------- Run inference for the learning model using the learning data--------
InferlearningModel = function(ntrials, G){
data <- list("G", "ntrials")
params <- c("theta", "theta1", "alpha")

samples <- jags(data, inits=NULL, params, model.file = "chick_learn3.txt", 
                n.chains = 3, n.iter = 5000, n.burnin = 1000, n.thin =1)

return(samples)
}

# Only look at trial 2 
x <- samples$BUGSoutput$sims.list$theta[,2]

#--------------- Parameter recovery for fixed model---------------------
paramrecoverfixed = function(ntrials, nsims){
  
  truetheta <- array(0,c(nsims))
  inferredtheta <- array(0, c(nsims))
  for (i in 1:nsims){
    theta <- runif(1,0,1)
    G=SimfixedModel(ntrials, theta)
    truetheta[i] <- theta
    
    samples = InferfixedModel(ntrials, G)
    # Theta MAP (
    theta_pos = samples$BUGSoutput$sims.list$theta
    inferredtheta[i] <- density(theta_pos)$x[which(density(theta_pos)$y==max(density(theta_pos)$y))]
  }
  tplot=plot(truetheta, inferredtheta)
  return(tplot)
}

#--------------- Parameter recovery for learning model---------------------
paramrecoverlearning= function(ntrials, nsims){
truealpha <- array(0,c(nsims))
inferredalpha <- array(0, c(nsims))

truetheta1 <- array(0,c(nsims))
inferredtheta <- array(0, c(nsims))

for (i in 1:nsims){
  
  #--------- Model 2: Learning model -----------------
  Glearn <- array(0, c(ntrials))
  theta <- array(0, c(ntrials))
  
  alpha <- runif(1,0,1)
  theta1 <- runif(1,0,1)
  
  # First trial
  theta[1] <- theta1
  Glearn[1] <- rbinom(1, 1, theta[1])
  
  # The rest of the trials
  for (t in 2:ntrials) {
    theta[t] <- theta[t-1]^(1/(1+alpha)) # Learning function
    Glearn[t] <- rbinom(1, 1, theta[t])
  }
  
  truealpha[i] <-alpha
  truetheta1[i] <- theta1
  
  #Inference
  G <- Glearn
  data <- list("G", "ntrials")
  params <- c("theta", "theta1", "alpha")
  
  samples <- jags(data, inits=NULL, params, model.file = "chick_learn3.txt", 
                  n.chains = 3, n.iter = 5000, n.burnin = 1000, n.thin =1)
  
  # Get MAP of alpha
  alpha = samples$BUGSoutput$sims.list$alpha
  inferredalpha[i] <- density(alpha)$x[which(density(alpha)$y==max(density(alpha)$y))]
  
  
  # Theta MAP (averaged)
  theta_pos = samples$BUGSoutput$sims.list$theta1
  inferredtheta[i] <- density(theta_pos)$x[which(density(theta_pos)$y==max(density(theta_pos)$y))]
  
}

results = cbind(truetheta1, inferredtheta, truealpha, inferredalpha)

return(results)
}

# ---------------------- Model recovery -------------------------------
# We use information criteria to evaluate model recovery as we cannot evaluate against a ground truth
# Simulate data from each model. Fit the model to data from both models

modelrecover = function(niter){
  
  fixtofix_DIC=array(0, c(niter))
  fixtolearn_DIC=array(0, c(niter))
  learntofix_DIC=array(0, c(niter))
  learntolearn_DIC=array(0, c(niter))
  
  for (i in 1:niter){
    theta <- runif(1,0,1)
    fix_dat= SimfixedModel(100, theta)
    alpha <- runif(1,0,1)
    learn_dat = SimlearningModel(100, theta, alpha)
    
    fixtofix=InferfixedModel(100, fix_dat)
    fixtolearn=InferfixedModel(100, learn_dat)
    learntofix=InferlearningModel(100, fix_dat)
    learntolearn=InferlearningModel(100, learn_dat)
    
    fixtofix_DIC[i]=fixtofix$BUGSoutput$DIC
    fixtolearn_DIC[i]=fixtolearn$BUGSoutput$DIC
    learntofix_DIC[i]=learntofix$BUGSoutput$DIC
    learntolearn_DIC[i]=learntolearn$BUGSoutput$DIC
  }

  results=data.frame(fixtofix_DIC, fixtolearn_DIC, learntofix_DIC, learntolearn_DIC)
  return(results)
}

#Run function
test=modelrecover(100)

#Confusion matrix
#Create dataframe for each model
dic_fixed_dat <- data.frame(f_mod = test$fixtofix_DIC, lmod = test$learntofix_DIC)
dic_learning_dat <- data.frame(f_mod = test$fixtolearn_DIC, lmod = test$learntolearn_DIC)
# Find minimum
min_dic_fdat <- apply(dic_fixed_dat, 1, which.min)
min_dic_ldat <- apply(dic_learning_dat, 1, which.min)

# Make predictions
n_sims = 100
dic_df <- data.frame(predictions = c(min_dic_fdat, min_dic_ldat), true = c(rep('Fixed', n_sims), rep('Learning', n_sims))) %>% 
  mutate(predictions = as.character(predictions)) %>%
  mutate(predictions = recode(predictions, '1' = 'Fixed', '2' = 'Learning'))

conf_mat <-  confusion_matrix(dic_df$true, dic_df$predictions)
cf_mat <- conf_mat$`Confusion Matrix`[[1]]

plot_confusion_matrix(cf_mat, add_row_percentages = F, add_col_percentages = F, add_normalized = F, counts_on_top = F)

