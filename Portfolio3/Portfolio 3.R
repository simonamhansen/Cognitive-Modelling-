pacman::p_load(extraDistr, R2jags, ggplot2)
setwd("C:/Users/Simon/Google Drev/Uni/Cognitive modelling")
set.seed(1995)


#---------------------------- Create a task environment -------------------
# Create rewards and losses for each deck (10 trials)
A_R <- rep(100, 10)
A_L <- c(rep(-250,5), rep(0, 5))

B_R <- rep(100, 10)
B_L <- c(-1250, rep(0, 9))

C_R <- rep(50, 10)
C_L <- c(rep(-50, 5), rep(0, 5))
  
D_R <- rep(50, 10)
D_L <- c(rep(-250, 1), rep(0,9))

# Empty lists for each deck
A <- c()
B <- c()
C <- c()
D <- c()

# Turn the blocks of 10 into 100
for (i in 1:10){
  A <- append(A, A_R+sample(A_L))
  B <- append(B, B_R+sample(B_L))
  C <- append(C, C_R+sample(C_L))
  D <- append(D, D_R+sample(D_L))
}

payoff = cbind(A, B, C, D)

#------------------- PVL-Delta model ------------------

# Set parameters

w <- 0.68
A <- 0.88
theta <- 3.04
a <- 0.25

ntrials <- 100

PVL_D = function(payoff, ntrials, w, A, a, theta){
  x <- array(0, c(ntrials)) #choice
  r <- array(0, c(ntrials)) #outcome
  u <- array(0, c(ntrials, 4)) #utility for each deck
  Ev <- array(0, c(ntrials, 4)) # Expected utility
  Ev_update <- array(0, c(ntrials, 4)) # Trick to update expected utility 
  exp_p <- array(0, c(ntrials, 4)) # Exponentiated Ev
  p <- array(0, c(ntrials, 4)) # probabilities
  
  x[1] <- rcat(1, c(0.25, 0.25, 0.25, 0.25))
  r[1] <- payoff[1, x[1]]
  
  #------------ plot prospect theory ---------------
  #w <- 2 
  #A <- .5
  #objective_val <- seq(-100, 100, 1)
  #subjective_ut <- ifelse(objective_val > 0, objective_val^A, -w*abs(objective_val)^A)
  #plot(objective_val,subjective_ut)
  
  for (t in 2:ntrials){
    
    for (d in 1:4){
      u[t,d] <- ifelse(r[t-1] < 0, -w*abs(r[t-1])^A, (abs(r[t-1]))^A)
      
      Ev_update[t,d] <- Ev[t-1,d]+ a*(u[t,d] - Ev[t-1,d])
      
      Ev[t,d] <- ifelse(x[t-1] ==d, Ev_update[t,d], Ev[t-1,d])
      
      exp_p[t,d] <- exp(theta*Ev[t,d])
    }
    
    for (d in 1:4){
      
      p[t,d] <- exp_p[t,d]/sum(exp_p[t,])
      
    }
    
    x[t] <- rcat(1, p[t,])
    r[t] <- payoff[t,x[t]]
  }
  
  results <- list(x=x, r=r, Ev=Ev)
  
}


payoff_s = payoff/100

PVL.sims=PVL_D(payoff_s, ntrials, w, A, a, theta)

plot(PVL.sims$x)
plot(PVL.sims$r)


par(mfrow = c(1,1))
plot(PVL.sims$Ev[,1])
plot(PVL.sims$Ev[,2])
plot(PVL.sims$Ev[,3])
plot(PVL.sims$Ev[,4])


# JAGS
x <- PVL.sims$x
r <- PVL.sims$r

data <- list("x", "r", "ntrials")
params <- c("w", "A", "theta", "a")

samples <- jags.parallel(data, inits = NULL, params, 
                                       model.file = "PVLD.txt",
                                       n.chains = 3, n.iter = 5000, n.burnin = 1000, n.thin = 1)

#Traceplot
traceplot(samples)


# -------------- ORL model ------------------

ORL = function(payoff, ntrials, Arew, Apun, K, wp, wf, theta){
  x <- array(0, c(ntrials)) #choice
  r <- array(0, c(ntrials)) #outcome
  signX <- array(0,c(ntrials)) #sign of outcome
  Ev <- array(0, c(ntrials, 4)) # Expected utility
  Ev_update <- array(0, c(ntrials, 4)) # Trick to update expected utility 
  Ef <- array(0, c(ntrials, 4)) # Expected frequency
  Ef_chosen <- array(0, c(ntrials, 4)) # Expected frequency for chosen deck
  Ef_notchosen <- array(0, c(ntrials, 4)) # Same as above for unchosen decks
  PS <- array(0, c(ntrials, 4)) # Perseverance value for each deck
  v <- array(0, c(ntrials, 4)) # Valence (combined score of EV, Ef and perseverance)
  exp_p <- array(0, c(ntrials, 4)) # Exponentiated Valence
  p <- array(0, c(ntrials, 4)) # probabilities
  
  x[1] <- rcat(1, c(0.25, 0.25, 0.25, 0.25))
  r[1] <- payoff[1, x[1]]
  
  for (t in 2:ntrials){
    signX[t] <-ifelse(r[t-1]==0, 0, ifelse(r[t-1]<0, -1, 1)) 
    for (d in 1:4){
      
      Ev_update[t,d] <- ifelse(r[t-1] >= 0, Ev[t-1,d]+ Arew*(r[t-1] - Ev[t-1,d]), Ev[t-1,d]+ Apun*(r[t-1] - Ev[t-1,d]))
      
      Ev[t,d] <- ifelse(x[t-1] ==d, Ev_update[t,d], Ev[t-1,d])
      
      Ef_chosen[t,d] <- ifelse(r[t-1] >= 0, Ef[t-1,d]+ Arew*(signX[t] - Ef[t-1,d]), Ef[t-1,d]+ Apun*(signX[t] - Ef[t-1,d]))
      
      Ef_notchosen[t,d] <- ifelse(r[t-1] >= 0, Ef[t-1,d]+ Arew*((-signX[t])/3 - Ef[t-1,d]), Ef[t-1,d]+ Apun*((-signX[t])/3 - Ef[t-1,d])) 
      
      Ef[t,d] <- ifelse(x[t-1] ==d, Ef_chosen[t,d], Ef_notchosen[t,d])
      
      PS[t,d] <- ifelse(x[t-1] ==d, 1/(1+K), PS[t-1,d]/(1+K))
      
      v[t,d] <- Ev[t,d] + wf*Ef[t,d] + wp*PS[t,d]
      
      exp_p[t,d] <- exp(theta*v[t,d])
      
    }
    
    for (d in 1:4){
      
      p[t,d] <- exp_p[t,d]/sum(exp_p[t,])
      
    }
    x[t] <- rcat(1, p[t,])
    r[t] <- payoff[t,x[t]]
  }
  results <- list(x=x, r=r, Ev=Ev, Ef = Ef, PS = PS , p = p, exp_p = exp_p)
}

# Set parameters
Arew = 0.9
Apun = 0.9
K = 0.5
wp = -3
wf = 1
theta = 1

ntrials = 100

payoff_s = payoff/100


ORL.sims=ORL(payoff_s, ntrials, Arew, Apun, K, wp, wf, theta)

par(mfrow= c(2,2))

plot(ORL.sims$PS[,1])
plot(ORL.sims$PS[,2])
plot(ORL.sims$PS[,3])
plot(ORL.sims$PS[,4])

ORL.sims$x
plot(ORL.sims$Ev[,4])
plot(ORL.sims$Ef[,1])
plot(ORL.sims$PS[,1])


# JAGS
x <- ORL.sims$x
r <- ORL.sims$r

data <- list("x", "r", "ntrials")
params <- c("Arew", "Apun", "K", "wp", "wf")

samples <- jags.parallel(data, inits = NULL, params, 
                         model.file = "ORL.txt",
                         n.chains = 3, n.iter = 5000, n.burnin = 1000, n.thin = 1)

#Traceplot
traceplot(samples)


# ------------ VSE model ---------------

# Set up the VSE payoff matrix
# extract rewards from reward trials
RA <- A
RA[A < 0] = 100
RB <- B
RB[B < 0] = 100
RC <- C
RC[C <= 0] = 50
RD <- D
RD[D < 0] = 50

# extract losses from reward trials
LA <- A
LA[A > 0] = 0
LA[LA < 0] = -250
LB <- B
LB[B > 0] = 0
LB[LB < 0] = -1250
LC <- C
LC[C == 0] = -50
LC[LC > 1] = 0
LD <- D
LD[D > 0] = 0
LD[LD < 0] = -250

R <- cbind(RA,RB,RC,RD)
L <- cbind(LA,LB,LC,LD)
L <- abs(L)
R <- R

ntrials = 100

# Set parameter values 
theta = 0.5 # Value sensitivity / Risk parameter
decay = 0.5 # Decay parameter 
alfa = 0.2 # Learning rate
phi =  2 # Exploration bonus 
C = 3 # Consistency parameter

# Set up function
VSE = function(R, L, ntrials, theta, decay, alfa, phi, C) {
  x <- array(0, c(ntrials)) #choice
  r <- array(0, c(ntrials)) #reward
  l <- array(0, c(ntrials)) #loss
  v <- array(0, c(ntrials)) #utility for the last outcome
  Exploit <- array(0, c(ntrials, 4)) # How much you want to exploit each deck
  Explore <- array(0, c(ntrials, 4)) # How much you want to explore each deck
  exp_p <- array(0, c(ntrials, 4)) # Exponentiated Ev
  p <- array(0, c(ntrials, 4)) # probabilities
  
  x[1] <- rcat(1, c(0.25, 0.25, 0.25, 0.25))
  r[1] <- R[1, x[1]]
  l[1] <- L[1, x[1]]
  
  for (t in 2:ntrials){
    
    v[t-1] = r[t-1]^theta - l[t-1]^theta
    
    for (d in 1:4){
      
      Exploit[t,d] = ifelse(x[t-1] == d, Exploit[t-1,d]*decay + v[t-1], Exploit[t-1,d]*decay)
      
      Explore[t,d] = ifelse(x[t-1] == d, 0, Explore[t-1, d]+ alfa*(phi -Explore[t-1, d]))
      
      exp_p[t,d] <- exp((Exploit[t,d] + Explore[t,d])*C)
    }
    
    for (d in 1:4){
      
      p[t,d] <- exp_p[t,d]/sum(exp_p[t,])
      
    }
    x[t] <- rcat(1, p[t,])
    r[t] <- R[t,x[t]]
    l[t] <- L[t,x[t]]
    
  }
  
  results <- list(x=x, r=r, l=l, Exploit=Exploit, Explore = Explore)
}

VSE.sims = VSE(R, L, ntrials, theta, decay, alfa, phi, C)

plot(VSE.sims$Explore[,3])

# JAGS stuff
x <- VSE.sims$x
r <- VSE.sims$r
l <- VSE.sims$l

data <- list("x", "r", 'l', "ntrials")
params <- c("theta", "decay", "alfa", "phi", "C")

samples <- jags.parallel(data, inits = NULL, params, 
                         model.file = "VSE.txt",
                         n.chains = 3, n.iter = 5000, n.burnin = 1000, n.thin = 1)

#Traceplot
traceplot(samples)

# ------------Parameter recovery function--------------

param_recover = function(niter,data, params, paramlist, sims.function, model.filename){
  
  # Empty arrays
  true_param = array(0, c(niter, length(paramlist)))
  infer_param = array(0, c(niter, length(paramlist)))
  
  # Loop through each iteration
  for (n in 1:niter){
  
    # Save true parameter values
    for (i in 1:length(paramlist)){
      run = paramlist[i]
      true_param[n,i] = eval(parse(text = run))
    }
    
    # run function
    sims.list = eval(sims.function)
  
    x = sims.list$x
    r = sims.list$r
    #l = sims.list$l
  
    samples <- jags.parallel(data, inits = NULL, params, 
                           model.file = model.filename,
                           n.chains = 3, n.iter = 5000, n.burnin = 1000, n.thin = 1)
    
    for (i in 1:length(params)){
    
      nam = params[i]
      
      i_pos = samples$BUGSoutput$sims.list[[nam]]
    
      infer_param[n, i] <- density(i_pos)$x[which(density(i_pos)$y==max(density(i_pos)$y))]
    }
    print(n)
  }
  
  results = list(true_param, infer_param)
  return(results)
}


# Perform parameter recovery for each mdodel 

#PVL-D
niter = 100
data = list("x", "r", "ntrials")
params = c("w", "A", "theta", "a")
paramlist = c('w = runif(1, 0, 5)', 'A = runif(1, 0, 1)', 'theta = runif(1, 0, 5)', 'a = runif(1, 0, 1)')
true_param = array(0, c(niter, 4))

# Run function
PVL_recover=param_recover(niter, data, params, paramlist, parse(text = "PVL_D(payoff_s, ntrials, true_param[n,1], true_param[n,2], true_param[n,4], true_param[n,3])"), "PVLD.txt")

# Make plots
true_PVL = PVL_recover[1]
infer_PVL = PVL_recover[2]

true_PVL2 = true_PVL[[1]]
infer_PVL2 = infer_PVL[[1]]

par(mfrow=c(2,2))
plot(true_PVL2[,1], infer_PVL2[,1]) + abline(0,1) + title('w, loss aversion')
plot(true_PVL2[,2], infer_PVL2[,2]) + abline(0,1) + title('A, risk preference')
plot(true_PVL2[,3], infer_PVL2[,3]) + abline(0,1) + title('theta, inverse heat')
plot(true_PVL2[,4], infer_PVL2[,4]) + abline(0,1) + title('a, learning rate')


#ORL 
niter = 100
data = list("x", "r", "ntrials")
params = c("Arew", "Apun", "K", "wp", "wf")
paramlist = c('Arew = runif(1, 0, 1)', 'Apun = runif(1,0,1)', 'K = runif(1, 0, 5)', 'wp = runif(1, -5, 5)', 'wf = runif(1, -5, 5)')
true_param = array(0, c(niter, 5))

ORL_recover = param_recover(niter, data, params, paramlist, parse(text = 'ORL(payoff_s, ntrials, true_param[n,1], true_param[n,2], true_param[n,3], true_param[n,4], true_param[n,5], 1)'), 'ORL.txt')

# Make plots
true_ORL = ORL_recover[1]
infer_ORL = ORL_recover[2]

true_ORL2 = true_ORL[[1]]
infer_ORL2 = infer_ORL[[1]]

par(mfrow=c(3,2))
plot(true_ORL2[,1], infer_ORL2[,1]) + abline(0,1) + title("Arew")
plot(true_ORL2[,2], infer_ORL2[,2]) + abline(0,1) + title("Apun")
plot(true_ORL2[,3], infer_ORL2[,3]) + abline(0,1) + title("K")
plot(true_ORL2[,4], infer_ORL2[,4]) + abline(0,1) + title("wp")
plot(true_ORL2[,5], infer_ORL2[,5]) + abline(0,1) + title("wf")

# VSE
niter = 100
data = list("x", "r", 'l', "ntrials")
params = c("theta", "decay" , "alfa", "phi", "C")
paramlist = c('theta = runif(1, 0, 1)', 'decay = runif(1, 0, 1)', 'alfa = runif(1, 0, 1)', 'phi = runif(1, -2.5, 2.5)', 'C = runif(1, 0, 5)')
true_param = array(0, c(niter, 5))
R=R/100
L=L/100

VSE_recover=param_recover(niter, data, params, paramlist, parse(text = 'VSE(R, L, ntrials, true_param[n,1], true_param[n,2], true_param[n,3], true_param[n,4], true_param[n,5])'), 'VSE.txt')

# Make plots
true_VSE = VSE_recover[1]
infer_VSE = VSE_recover[2]

true_VSE2 = true_VSE[[1]]
infer_VSE2 = infer_VSE[[1]]

par(mfrow=c(3,2))
plot(true_VSE2[,1], infer_VSE2[,1]) + abline(0,1) + title("theta")
plot(true_VSE2[,2], infer_VSE2[,2]) + abline(0,1) + title("decay")
plot(true_VSE2[,3], infer_VSE2[,3]) + abline(0,1) + title("alfa")
plot(true_VSE2[,4], infer_VSE2[,4]) + abline(0,1) + title("phi")
plot(true_VSE2[,5], infer_VSE2[,5]) + abline(0,1) + title("C")
