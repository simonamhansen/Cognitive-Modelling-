# Set seed, working directory and load relevant libaries 
set.seed(1995)
pacman::p_load(R2jags, polspline, extraDistr)
setwd("C:/Users/Simon/Google Drev/Uni/Cognitive modelling")

# -------------- Environment for EWA ----------------
# Contextual variables // evironment variables
ntrials = 15
nagents = 3
ntokens = 20
niterations = 100
pi = 1.5

# True variables
true = c()
true$delta = array(0, c(nagents, niterations))
true$rho = array(0, c(nagents, niterations))
true$phi = array(0, c(nagents, niterations))
true$lambda = array(0, c(nagents, niterations))

# Inferred variables
MAP = c()
MAP$delta = array(0, c(nagents, niterations))
MAP$rho = array(0, c(nagents, niterations))
MAP$phi = array(0, c(nagents, niterations))
MAP$lambda = array(0, c(nagents, niterations))


#--------------- EWA function and JAGS -----------------------

EWA_fun = function(nagents, ntrials, ntokens, pi, parameters){
  
  # Empty arrays
  c = array(0, c(nagents, ntrials)) # Contributed amount
  N = array(0, c(nagents, ntrials)) # Experience 
  A = array(0, c(nagents, ntrials, ntokens)) # Attractions
  expA = array(0, c(nagents, ntrials, ntokens)) # Softmax input 
  P = array(0, c(nagents, ntrials, ntokens)) # Probability
  
  
  delta = array(0, c(nagents))
  rho = array(0, c(nagents))
  phi = array(0, c(nagents))
  lambda = array(0, c(nagents))
  
  # Parameters
  for (n in 1:nagents){  
    delta[n] = parameters$delta[n]
    rho[n] = parameters$rho[n]
    phi[n] = parameters$phi[n]
    lambda[n] = parameters$lambda[n]
  }
  
  # Trial 1
  
  N[,1] = 1
  A[,1,] = 0
  c[,1] = 20
  
  # Trial 2 and onwards
  
  for (t in 2:ntrials){
    
    for (n in 1:nagents){
      
      N[n,t] = (rho[n]*N[n,t-1]) + 1 # Experience updating
      
      for (j in 1:ntokens){
        
        A[n,t,j] = (
          
          (phi[n]*N[n,t-1]*A[n,t-1,j]) + # Prior attractions
          (delta[n] + ((1-delta[n])*(c[n,t-1]==j))) * # Was Jth token chosen
          ((((j+sum(c[-n,t-1]))*pi)/nagents)-j) # Payoff for each possible contribution
          
          )/
          N[n,t] # experience weighting
        
        expA[n,t,j] = exp(lambda[n]*A[n,t,j])
      }
      
      for (j in 1:ntokens){
        
        P[n,t,j] = expA[n,t,j]/sum(expA[n,t,])
      }
      
      c[n,t] = rcat(1, P[n,t,])
      
    }
  }
  
  results = list(c=c, N=N)
  
  return(results)
}

# Test the EWA function
parameters = data.frame(delta = c(0.2, 0.2, 0.2), rho = c(0.2, 0.2, 0.2), phi = c(0.8, 0.8, 0.8), lambda = c(0.1, 0.1, 0.1))

EWA_sims=EWA_fun(nagents, ntrials, ntokens, pi, parameters)

# Postprocessing
c = EWA_sims$c

Gc = array(0, c(nagents, ntrials))
for (n in 1:nagents){
  Gc[n,] <- colSums(c[-n,])
}

c_choice_index = c

# JAGS
data = list('nagents', 'ntrials', 'ntokens', 'pi', 'c', 'Gc', 'c_choice_index')
params <- c('delta', 'rho', 'phi', 'lambda')
samples <- jags(data, inits = NULL, params,
                model.file = 'EWA.txt',
                n.chains =3, n.iter = 5000, n.burnin =  1000, n.thin = 1)

traceplot(samples)

plot(density(samples$BUGSoutput$sims.list$rho))

# ------------ Parameter recovery for EWA -----------------

delta_i = array(0, c(nagents, 12000))
rho_i = array(0, c(nagents, 12000))
phi_i = array(0, c(nagents, 12000))
lambda_i = array(0, c(nagents, 12000))

for (i in 1:niterations){
  
  # True parameters
  true$delta[,i] = runif(3, 0, 1)
  true$rho[,i] = runif(3, 0, 1)
  true$phi[,i] = runif(3, 0, 1)
  true$lambda[,i] = runif(3, 0, 5)
  
  # Forward simulation
  parameters = data.frame(delta = true$delta[,i], rho = true$rho[,i], phi = true$phi[,i], lambda = true$lambda[,i])
  
  EWA_sims=EWA_fun(nagents, ntrials, ntokens, pi, parameters)
  
  # Postprocessing
  c = EWA_sims$c
  
  Gc = array(0, c(nagents, ntrials))
  for (n in 1:nagents){
    Gc[n,] <- colSums(c[-n,])
  }
  
  c_choice_index = c
  
  # JAGS
  data = list('nagents', 'ntrials', 'ntokens', 'pi', 'c', 'Gc', 'c_choice_index')
  params <- c('delta', 'rho', 'phi', 'lambda')
  samples <- jags.parallel(data, inits = NULL, params,
                  model.file = 'EWA.txt',
                  n.chains =3, n.iter = 5000, n.burnin =  1000, n.thin = 1)
  
  
  # MAP 
  for (n in 1:nagents){
    delta_i[n,] = samples$BUGSoutput$sims.list$delta[,n]
    rho_i[n,] = samples$BUGSoutput$sims.list$rho[,n]
    phi_i[n,] = samples$BUGSoutput$sims.list$phi[,n]
    lambda_i[n,] = samples$BUGSoutput$sims.list$lambda[,n]
  
    MAP$delta[n,i]=density(delta_i[n,])$x[which(density(delta_i[n,])$y==max(density(delta_i[n,])$y))]
    MAP$rho[n,i] = density(rho_i[n,])$x[which(density(rho_i[n,])$y==max(density(rho_i[n,])$y))]
    MAP$phi[n,i] = density(phi_i[n,])$x[which(density(phi_i[n,])$y==max(density(phi_i[n,])$y))]
    MAP$lambda[n,i] = density(lambda_i[n,])$x[which(density(lambda_i[n,])$y==max(density(lambda_i[n,])$y))]
  }
  print(i)
}

par(mfrow = c(2,2))
plot(true$delta, MAP$delta) + abline(0,1) + title("delta")
plot(true$rho, MAP$rho) + abline(0,1) + title("rho")
plot(true$phi, MAP$phi) + abline(0,1) + title("phi")
plot(true$lambda, MAP$lambda) + abline(0,1) + title("lambda")


# ------- Conditional cooperation model. Forward simulation and JAGS------

# Environment for CC model 
ntrials = 15
nagents = 3
niterations = 100
vals = seq(1,20,1) # possible values

# Generated/True parameters
trueCC = c()
trueCC$omega1 = array(0, c(nagents, niterations)) # Weighting for beliefs vs preferences on first trial
trueCC$lambda = array(0, c(nagents, niterations)) # Decay rate for beliefs
trueCC$gamma = array(0, c(nagents, niterations)) # Learning rate for beliefs
trueCC$pbeta = array(0, c(nagents, niterations)) # Slope relating the contribution of others to your own preference
trueCC$p0 = array(0, c(nagents, niterations))

# Inferred parameters
MAPCC = c()
MAPCC$omega1 = array(0, c(nagents, niterations))
MAPCC$lambda = array(0, c(nagents, niterations))
MAPCC$gamma = array(0, c(nagents, niterations))
MAPCC$pbeta = array(0, c(nagents, niterations))
MAPCC$p0 = array(0, c(nagents, niterations))

CC_sim = function(nagents, ntrials, vals, parameters){
  # Free parameters
  Gb1 = parameters$Gb1
  omega1 =parameters$omega1
  lambda = parameters$lambda
  gamma = parameters$gamma
  p0 = parameters$p0
  pbeta = parameters$pbeta
  
  # Simulation arrays
  Ga = array(0, c(ntrials))
  Gb = array(0, c(nagents, ntrials))
  p = array(0, c(nagents,ntrials))
  omega = array(0, c(nagents, ntrials))
  c = array(0, c(nagents, ntrials))
  
  pvals = array(0, c(nagents, length(vals)))
  for (n in 1:nagents){
    pvals[n,] = p0[n] + (pbeta[n]*vals)
  }
  
  omega[,1] = omega1
  
  Gb[,1] = Gb1
  
  c[,1] = Gb1
  Ga[1] = mean(Gb1)
  
  for (t in 2:ntrials){
    
    for (n in 1:nagents){
      
      Gb[n,t] = gamma[n]*Gb[n,t-1]+(1-gamma[n])*Ga[t-1]
      
      p[n,t] = pvals[n, round(Gb[n,t])]
      
      omega[n,t] = omega[n,t-1]*(1-lambda[n])
      
      c[n,t] = ceiling(omega[n,t]*Gb[n,t] + (1-omega[n,t])*p[n,t])
    }
    
    Ga[t] = sum(c[,t])/nagents
    
  }
  
  results = list(c=c, omega = omega, Ga = Ga)
  
  return(results)
}

# Test forward simulation
parameters = data.frame(omega1 = c(0.2, 0.2, 0.2), lambda = c(0.2, 0.2, 0.2), gamma = c(0.8, 0.8, 0.8), p0 = c(10, 10, 10), pbeta = c(0.5, 0.5, 0.5), Gb1 = c(10, 10, 10))

CC = CC_sim(nagents, ntrials, vals, parameters) 
# Print contribution of each agent
CC$c

# JAGS
c = CC$c
Ga = CC$Ga

data = list('nagents', 'ntrials', 'vals', 'c', 'Ga')
params <- c('omega1', 'lambda', 'gamma', 'p0', 'pbeta', 'c', 'omega')
samples <- jags(data, inits = NULL, params,
                model.file = 'CC.txt',
                n.chains =3, n.iter = 5000, n.burnin =  1000, n.thin = 1)

traceplot(samples)

plot(density(samples$BUGSoutput$sims.list$lambda[,1]))

# -------- Parameter recovery for CC model ------
omega1_i = array(0, c(nagents, 12000))
lambda_i = array(0, c(nagents, 12000))
gamma_i = array(0, c(nagents, 12000))
pbeta_i = array(0, c(nagents, 12000))
p0_i= array(0, c(nagents, 12000))

for (i in 1:niterations){
  
  # True parameters
  trueCC$omega1[,i] = runif(3, 0, 1)
  trueCC$lambda[,i] = runif(3, 0, 1)
  trueCC$gamma[,i] = runif(3, 0, 1)
  trueCC$p0[,i] = runif(3, 0, 20)
  trueCC$pbeta[,i] = runif(3, 0, (1-0.05*trueCC$p0[,i]))
  
  # Forward simulation
  parameters = data.frame(omega1 = trueCC$omega1[,i], lambda = trueCC$lambda[,i], gamma = trueCC$gamma[,i], p0 = trueCC$p0[,i], pbeta = trueCC$pbeta[,i], Gb1 = c(10, 10, 10))
  
  sims=CC_sim(nagents, ntrials, vals, parameters)
  
  # Extract choices
  c=sims$c
  Ga = sims$Ga
  
  # JAGS
  data = list('nagents', 'ntrials', 'vals', 'c', 'Ga')
  params <- c('omega1', 'lambda', 'gamma', 'p0', 'pbeta', 'c', 'omega')
  samples <- jags(data, inits = NULL, params,
                  model.file = 'CC.txt',
                  n.chains =3, n.iter = 5000, n.burnin =  1000, n.thin = 1)
  
  for (n in 1:nagents){
    omega1_i[n,] = samples$BUGSoutput$sims.list$omega1[,n]
    lambda_i[n,] = samples$BUGSoutput$sims.list$lambda[,n]
    gamma_i[n,] = samples$BUGSoutput$sims.list$gamma[,n]
    p0_i[n,] = samples$BUGSoutput$sims.list$p0[,n]
    pbeta_i[n,] = samples$BUGSoutput$sims.list$pbeta[,n]
    
    
    MAPCC$omega1[n,i]=density(omega1_i[n,])$x[which(density(omega1_i[n,])$y==max(density(omega1_i[n,])$y))]
    MAPCC$lambda[n,i] = density(lambda_i[n,])$x[which(density(lambda_i[n,])$y==max(density(lambda_i[n,])$y))]
    MAPCC$gamma[n,i] = density(gamma_i[n,])$x[which(density(gamma_i[n,])$y==max(density(gamma_i[n,])$y))]
    MAPCC$p0[n,i] = density(p0_i[n,])$x[which(density(p0_i[n,])$y==max(density(p0_i[n,])$y))]
    MAPCC$pbeta[n,i] = density(pbeta_i[n,])$x[which(density(pbeta_i[n,])$y==max(density(pbeta_i[n,])$y))]
    }
  print(i)
  
}

par(mfrow=c(3,2))
plot(trueCC$omega1, MAPCC$omega1) + abline(0,1) + title("Omega(t=1)")
plot(trueCC$lambda, MAPCC$lambda) + abline(0,1) + title("Lambda")
plot(trueCC$gamma, MAPCC$gamma) + abline(0,1) + title("Gamma")
plot(trueCC$pbeta, MAPCC$pbeta) + abline(0,1) + title("pbeta")
plot(trueCC$p0, MAPCC$p0) + abline(0,1) + title("p0")
