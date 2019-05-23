########### ALTIAR Mussels ###############
setwd("C:/Users/vasot/Dropbox/Research/Texas_State/TXST_Research/Dissertation_Data/Mussels/Ricke_EnvModel/ALT_QUPE/CatModel")

# Front-end -----
# . Library load ----
library(R2jags)

# . Function definition ----
# Function for inverting the logit
inv.logit <- function(p)exp(p)/(1+exp(p))

#load("ALT_QUHO_DisPhi.RData")

# Data manipulation ----
# . Data read and definitions ----
# Capture histories
d<-read.csv("ALT_CH_QUPE.csv",header=F)

# Discharge data, NEEDS TO BE SAME LENGTH AS NUMBER OF MUSSELS
dis<-read.csv("MinDischargeDummy_ALT.csv",header=F)

str(dis)

# Make a matrix out of desired columns  

dis<-as.matrix(dis)


discharge <- array(data=dis,dim=c(120,4,3))
discharge[1:10,1:4,2]



#log
#discharge<-log(discharge)

#scale
#discharge<-scale(c(discharge))
#discharge <- structure(discharge, dim=c(308,4))

#standardize

dimnames(discharge) <- NULL

#################### NO CHANGES BELOW ###########################


# Make ch into a matrix
ch <- as.matrix(d)
dimnames(ch) <- NULL

# Get first captures for primaries
firsts<-function(x=x){min(which(x>0))}
fc <- apply(ch, 1, firsts)
fc <- ceiling(fc/3)

# Mark-recapture sampling regime
n.sec.occ <- 3                  # Secondary occasions
n.years <- 5                    # Primary occasions
n <- nrow(ch)                   # Individuals

# . Ragged array formatting ----  
### DSS: This is backward right now,
### but the rest of the code is, too.
### See comment below about yes/no matrices 

# Create 3-d array to hold observations
obs <- array(NA, c(n, n.years, n.sec.occ))

# Fill the observation array
# from the full ch
for (i in 1:(n)){
  for (t in fc[i]:n.years){
    for (j in 1:n.sec.occ){
      obs[i,t,j] <- ch[i, (t+n.years*(j-1))]
    }
  }
}
obs[is.na(obs)] <- 0

# . Primary period capture histories ----  
# Create aggregate capture histories
# for primary periods

ch2 <- cbind(apply(ch[,1:3], 1, sum),
             apply(ch[,4:6], 1, sum),
             apply(ch[,7:9], 1, sum),
             apply(ch[,10:12], 1, sum),
             apply(ch[,13:15], 1, sum))
ch2[ch2>0] <- 1


# . Secondary period capture histories ----
# Summarize detection for each
# secondary across all individuals.
test <- matrix(NA, n, n.years)
for (i in 1:nrow(test)){
  for (j in 1:ncol(test)){
    test[i,j] <- sum(obs[i,j,1:n.sec.occ])
  }
}

# Create seen/missed matrices used to build
# the yes/no/total matrices for binomial 
# detection model.
seen <- array(0, c(n, n.years, n.sec.occ))
missed <- array(0, c(n, n.years, n.sec.occ))

# Fill in the arrays  
for (i in 1:nrow(test)){
  for (t in 1:ncol(test)){
    for (j in 1:n.sec.occ){
      if(test[i,t] >= 1 & obs[i,t,j] == 1){seen[i,t,j] <- 1}
      if(test[i,t] >= 1 & obs[i,t,j] == 0){missed[i,t,j] <- 1}
    }
  }
}

# Create 2-d array to hold yes/no data
#yes <- matrix(apply(ch, 2, sum), ncol=3, byrow = TRUE)
yes <- matrix(NA, n.years, n.sec.occ)
no <- matrix(NA, n.years, n.sec.occ)


### THIS SEEMS TO BE FILLING BACKWARD
### LIKELY BECAUSE I BUILD RAGGED CH
### ARRAY BACKWARD AT START
for (i in 1:n.years){
  for (j in 1:n.sec.occ){
    yes[i,j] <- sum(seen[,i,j])
    no[i,j] <- sum(missed[,i,j])
  }
}
### SO, I AM UNLISTING AND RE-ORDERING FOR NOW
yes <- matrix(c(yes), ncol=3, byrow = TRUE)
no <- matrix(c(no), ncol=3, byrow = TRUE)
total <- yes + no


# Model specification ----
modelString = "
model {
  # Secondary occasions p
  for (t in 1:n.years){
    for (j in 1:(n.sec[t])){
      p[t,j] ~ dunif(0,1)
      yes[t,j] ~ dbin(p[t,j], total[t,j])
    }
  }   
  # Primary occasions p
  for (t in 1:n.years){
    pstar[t] <- 1 - prod(1 - p[t,])
  }
  # Time-invariant prior on gamma
    alpha.gamma ~ dnorm(0,0.001)
  
  ### environmental priors
  # Shared priors on the parameters of
  # linearized model for phi
    alpha.phi ~ dnorm(0, 0.001)
    beta.1 ~ dnorm(0, 0.001)
    beta.2 ~ dnorm(0, 0.001)
    beta.3 ~ dnorm(0, 0.001)
  
  # Linear predictors of phi and gamma
    for (i in 1:n.ind){
    for (t in 1:(n.years-1)){
       # Linearized models on phi and gamma
       # using logit link function. All parameters
       # returned will be on the link scale
       logit(phi[i,t]) <- alpha.phi + beta.1 * discharge[i,t,1] + beta.2 * discharge[i,t,2] + beta.3 * discharge[i,t,3]
       logit(gamma[i,t]) <- alpha.gamma
    }}
  
  # Likelihood
  for (i in 1:n.ind){
    
    z[i,first[i]] <- ch[i,first[i]]
    
    for (t in (first[i]+1):n.years){
      mu1[i,t] <- z[i,t-1] * phi[i,t-1]
      mu2[i,t] <- z[i,t] * gamma[i,t-1] * pstar[t]
      z[i,t] ~ dbern(mu1[i,t])
      ch[i,t] ~ dbern(mu2[i,t])
    } 
    
  } # end likelihood
  
} ## end model
"

# Model calibration ----
# . Bundle data ----
# Get first capture for primaries
first <- apply(ch2, 1, firsts)  



#################### NO CHANGES ABOVE ###########################

# Bundle the data for JAGS




##### NEED TO CHANGE DISCHARGE TO MATCH LENGTH  
dat <- list(first=first,
            ch = ch2,
            n.sec = rep(n.sec.occ, n.years),
            n.years = ncol(ch2),
            n.ind = nrow(ch2),
            yes = yes,
            total = total,
            discharge=discharge
)

# . Initial Values ----
z.init <- matrix(NA, nrow(ch2), ncol(ch2))
for (i in 1:nrow(ch2)){
  if(first[i] < ncol(z.init)){
    z.init[i, (first[i] + 1):ncol(z.init)] <- 1
  }
}

inits <- function(){list(z = z.init)}  

# . Parameters to monitor ----
pars <- c('pstar', 'p', 'alpha.phi','beta.1','beta.2','beta.3',
          'alpha.gamma')

# . MCMC settings ----
nc <- 3
nt <- 25
ni <- 30000
nb <- 5000


# . Call JAGS and run model ----
rdt <- jags(dat, inits=inits, pars,
            model.file=textConnection(modelString), 
            n.chains = nc, n.iter = ni,
            n.burnin = nb, n.thin = nt,
            refresh = 0.001
)


# . Print model summary ----

print(rdt, digits=3)

#rdt.min<-rdt
rdt.max<-rdt
#rdt.med<-rdt
