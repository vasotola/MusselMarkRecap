# Package load ----
library(tidyr)
library(dplyr)
library(reshape)
library(R2jags)
library(lubridate)

# Data manipulation ----
# . Tagging data ----
# Read in the tagging data
tags <- read.csv('data_alt_tags.csv', stringsAsFactors = FALSE)
check <- vector(mode = 'character', length=nrow(tags))

# Replace slashes in floy_id with letters
for(i in 1:nrow(tags)){
  tags$floy_id[i] <- gsub(
    pattern = "/", 
    x=tags$floy_id[i],
    replacement=substr(tags$floy_id[i], start=1, stop=1)
    )  
}

# Assign number of floy tags (1 or 2)
# Based on the length of the floy_id
# variable. All mussels have at least one.
n_floy_tags <- 1
n_floy_tags[nchar(tags$floy_id) > 4] <- 2
n_floy_tags[is.na(n_floy_tags)] <- 1

# Do the same for PIT tags
# Not all mussels were PIT tagged
n_pit_tags <- 0
n_pit_tags[nchar(tags$pit_id) > 0] <- 1
n_pit_tags[is.na(n_pit_tags)] <- 0

# Add an indicator variable for tagging
# configuration. Here, we paste together
# the number of each tag type with first
# letter of tag types to get the 4 unique
# combinations.
tags$tag_config <- paste0(
  n_floy_tags, "F",
  n_pit_tags, "P"
)

# Assign probabilities of tag loss for each
# tagging configuration.
pit_loss <- 0.053
floy_loss <- 0.024                        # "1F0P"
floy_2_loss <- floy_loss^2                # "2F0P"
pit_floy_loss <- pit_loss * floy_loss     # "1F1P"
pit_floy_2_loss <- pit_loss * floy_loss^2 # "2F1P"

# Add individual probabilities of tag loss
# to the tagging data
tags$p_tag_loss <- 0
tags$p_tag_loss[tags$tag_config=="1F0P"] <- floy_loss
tags$p_tag_loss[tags$tag_config=="2F0P"] <- floy_2_loss
tags$p_tag_loss[tags$tag_config=="1F1P"] <- pit_floy_loss
tags$p_tag_loss[tags$tag_config=="2F1P"] <- pit_floy_2_loss

# . Detection data ----
# Read in the floy recapture file
# from grubbing surveys
recap_floy <- read.csv('data_alt_floy.csv', stringsAsFactors = F)

# Read in the PIT detection file
# tag ID is indexed by floy_id
recap_pit <- read.csv('data_alt_pit.csv', stringsAsFactors = F)

# Convert both sets of detection data
# to long format

# Grubbing data
  # Convert
  rf_long <- tidyr::gather(
    data = recap_floy,
    key = date,
    value = ch,
    4:ncol(recap_floy)
    )
  # Replace slashes in floy_id with letters
  for(i in 1:nrow(rf_long)){
    rf_long$floy_id[i] <- gsub(
      pattern = "/", 
      x=rf_long$floy_id[i],
      replacement=substr(rf_long$floy_id[i], start=1, stop=1)
      )  
  }

# PIT data
  # Convert
  rp_long<- tidyr::gather(
    data = recap_pit,
    key = date,
    value = ch,
    4:ncol(recap_pit)
    )
  # Replace slashes in floy_id with letters
  for(i in 1:nrow(rp_long)){
    rp_long$floy_id[i] <- gsub(
      pattern = "/", 
      x=rp_long$floy_id[i],
      replacement=substr(rp_long$floy_id[i], start=1, stop=1)
      )  
  }


# Combine floy and PIT detections
dets <- rbind(rf_long, rp_long)

# There are three floy tags that 
# have IDs switched. 
flopped <- unique(dets$floy_id[!dets$floy_id%in%tags$floy_id])

# This is what they look like corrected
unflopped <- paste0(
  substr(flopped, start=5, stop=nchar(flopped)),
  substr(flopped, start=1, stop=4)
)

# Can replace them now
for(i in 1:length(flopped)){
  dets$floy_id[dets$floy_id==flopped[i]] <- unflopped[i]
}

# . Data merge ----
# Merge the tag data with detection data
ch_test <- merge(
  dets, tags, by = c('floy_id', 'species')
  )
nrow(ch_test) - nrow(dets)

ch_test <- ch_test[ch_test$species %in% c('QUHO', 'QUPE', 'AMPL'), ]

# Replace dates with unambiguous string format
# through serial replacement and string manipulation
# with pipes (%>%)
ch_test$date %>% 
  gsub(pattern='X', replacement = '') %>%
    gsub(pattern='\\.', replacement = '/')%>%
     as.character %>%
      as.Date(format = "%m/%d/%Y") -> ch_test$date

# Summarize covariates by floy_id
# in detection data
covs <- plyr::ddply(
  ch_test,
  c('floy_id', 'species', 'tag_config', 'p_tag_loss'),
  plyr::summarize, check <- 1
  )

# . Capture histories ----
# Make a capture history from combined
# detection data
ch2 <- data.frame(reshape::cast(
  data = ch_test,
  formula = floy_id ~ date,
  value = 'ch',
  fun.aggregate = sum
  ))


# Merge with the covs data for individual
# covariates
ch1 <- merge(ch2, covs, by = 'floy_id')

# Re-organize and convert all values in
# capture history > 1 to 1.
ch0 <- data.frame(
  cbind(
    
    ch1[ , c(1, (ncol(ch1) - 3):(ncol(ch1) - 1))],
    
    apply(
      ch1[ , 2:(ncol(ch1) - 4)],
      2,
      function (x) ifelse(x > 0, 1, 0)
      
    )
  )
)

# Have a look at raw counts for a sanity check
plot(apply(ch0[ ,5:ncol(ch0)], 2, sum), type='l')

# Make ch into a matrix
ch <- as.matrix(ch0[ , 5:ncol(ch0)])
dimnames(ch) <- NULL

# Get first captures for primaries
firsts<-function(x=x){min(which(x>0))}
fc <- apply(ch, 1, firsts)
fc <- ceiling(fc/3)

# Mark-recapture sampling regime
n.sec.occ <- 3                  # Secondary occasions
n.years <- 5                    # Primary occasions
n <- nrow(ch)                   # Individuals


# . Array formatting ----  
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
# Split into primaries
ch_mod <- cbind(apply(ch[,1:3], 1, sum),
             apply(ch[,4:6], 1, sum),
             apply(ch[,7:9], 1, sum),
             apply(ch[,10:12], 1, sum),
             apply(ch[,13:15], 1, sum))

# Designate as one if observed
ch_mod[ch_mod>0] <- 1


# . Secondary period capture histories ----
test <- matrix(NA, n, n.years)
for (i in 1:nrow(test)){
  for (j in 1:ncol(test)){
    test[i,j] <- sum(obs[i,j,1:n.sec.occ])
  }
}

# Create seen/missed matrices used to build
# the yes/no/total matrices for binomial 
# detection model (secondary occasions).
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
yes <- matrix(NA, n.years, n.sec.occ)
no <- matrix(NA, n.years, n.sec.occ)

### THIS SEEMS TO BE FILLING BACKWARD
### LIKELY BECAUSE I BUILD RAGGED CH
### ARRAY BACKWARD AT START. OKAY FOR
### ANALYSIS, JUST TOUGH TO PARALLEL
### REICKE (2016) CODE.
for (i in 1:n.years){
  for (j in 1:n.sec.occ){
    yes[i,j] <- sum(seen[,i,j])
    no[i,j] <- sum(missed[,i,j])
  }
}

# Unlist and re-order
yes <- matrix(c(yes), ncol=3, byrow = TRUE)
yesarray <- array(
  rep(yes, nrow(ch0)),
  dim = c(nrow(ch0), nrow(yes), ncol(yes))
  )
no <- matrix(c(no), ncol=3, byrow = TRUE)
total <- yes + no
totalarray <- array(
  rep(total, nrow(ch0)),
  dim = c(nrow(ch0), nrow(total), ncol(total))
  )

# Get counts for each species in each year for abundance
ch_test$primary <- paste0(
  lubridate::year(ch_test$date), "-",
  lubridate::month(ch_test$date)
)

# Species-specific counts in primary periods
# Assign a variable to hold primary period (month/season)
for(i in 1:nrow(ch_test)){
  if(nchar(ch_test$primary[i])<7){
    ch_test$primary[i] <- paste0(
      substr(ch_test$primary[i], start = 1, stop = 5),
      "0",
      substr(ch_test$primary[i], start = 6, stop = 6)
      )
  }
}

# Collapse individual detections on individuals,
# convert to presence/absence, and then collapse
# on species x primary period.
counts <- plyr::ddply(
  ch_test,
  c('species', 'primary', 'floy_id'),
  summarize,
  counter = sum(ch)
  )
counts$counter[counts$counter>1] <- 1
counts1 <- plyr::ddply(
  counts,
  c('species', 'primary'),
  summarize,
  counter = sum(counter)
  )
sp_totes <- data.frame(
  reshape::cast(counts1, formula = species~primary, value='counter')[ , -1]
)

# Get total number of individuals in each primary 
# period by species
sp_totes <- sp_totes[ , c(order(names(sp_totes)))]

# Get total number of individuals tagged for each
# species as a starting guess
sp_unique <- plyr::ddply(
  ch_test, 'species', plyr::summarise, totes=length(unique(floy_id))
  )

# Model specification ----
modelString = "
model {
  # Priors
    # Survival and emigration
      for(i in 1:(n.ind)){
        for(t in 1:(n.years-1)){
          # Logit-scale predictor of apparent survival (species)
            logit(r_phi[i,t]) <- l_phi_i[i,t]
            l_phi_i[i,t] <-  l_phi[species[i],t]
          
          # Logit-scale predictor of site fidelity (species)
            logit(r_gamma[i,t]) <- l_gamma_i[i,t]
            l_gamma_i[i,t] <- l_gamma[species[i],t]
        }
      }
    
    # Species- and time-specific priors on
    # phi and gamma
      for(s in 1:n.species){
        for(t in 1:(n.years-1)){
          l_phi[s,t] ~ dnorm(0, 0.001)
          l_gamma[s,t] ~ dnorm(0, 0.001)

          logit(phi[s,t]) <- l_phi[s,t]
          logit(gamma[s,t]) <- l_gamma[s,t]
        }
      }
    
    # Secondary occasions p
      # Logit-scale linear predictor on
      # p by year and secondary occasion
        for (i in 1:n.ind){
          for (t in 1:n.years){
            for (j in 1:(n.sec[t])){
                logit(p[i,t,j]) <- lp[i,t,j]  
                lp[i,t,j] <- inprod(beta[t,j,], X[i,])*(1-p_loss[i])
                yes[i,t,j] ~ dbin(p[i,t,j], total[i,t,j])
            }
          }
        }
      # Betas (intercepts & 'slopes') for linear predictor of logit-scale
      # p, varying by time-period
      for(t in 1:n.years){
        for(j in 1:(n.sec[t])){
          for(n in 1:nX){
            beta[t,j,n] ~ dnorm(0, 0.001)
          }
        }
      }

    # Primary occasions p
      for(i in 1:n.ind){
        for(t in 1:n.years){        
          pstar[i,t] <- 1 - prod(1 - p[i,t,])
        }
      }
  
  # State (presence/absence) and observation (detection)
  # process likelihoods
    for (i in 1:n.ind){

      # Initialize first capture
        z[i,first[i]] <- ch[i,first[i]]
  
      # Likelihood for state process
        for (t in (first[i]+1):n.years){
          mu1[i,t] <- z[i,t-1] * r_phi[i, t-1]
          mu2[i,t] <- z[i,t] * r_gamma[i, t-1] * pstar[i,t]
          z[i,t] ~ dbern(mu1[i,t])
          ch[i,t] ~ dbern(mu2[i,t])
        }
    }
    
  # Count process likelihood (abundance)
  # Population abundance estimated using n-mixture
  # approach. We assume that observed abundance
  # is the outcome (C) of a binomial process with some
  # unknown probability of success (pstar) and some
  # starting number of trials (N)
  for(s in 1:n.species){
    # Initialize counts for each species 
    N[s, 1] ~ dpois(maxes[s]/mean(pstar[,])) 

    log(lambda[s]) <- mu[s]
    mu[s] ~ dunif(-10, 10)

    # Total number alive at time t as function of
    # number alive at t-1, survival, and emigration
    for(t in 2:n.years){
      N[s,t] ~ dpois(N[s, t-1]*phi[s, t-1]*gamma[s,t-1])
    }

    for(t in 1:n.years){
      C[s, t] ~ dbin(mean(pstar[,t]), N[s,t])
    }

  }
  
}
"

# Model calibration ----
#. Bundle data ----
# Get first capture for primaries
first <- apply(ch_mod, 1, firsts)  

# Make a model matrix for tag configuration
# or other covariates of interest
dummy <- ch0[ , c(3:ncol(ch0))]
betas <- model.matrix(glm(p_tag_loss~tag_config, data = ch0))
attributes(betas)[2:3] <- NULL

# Package into a list for JAGS
dat <- list(
  first=first,
  ch = ch_mod,
  n.sec = rep(n.sec.occ, n.years),
  n.years = ncol(ch_mod),
  n.ind = nrow(ch_mod),
  n.species = length(unique(ch0$species)),
  species = as.numeric(as.factor(ch0$species)),
  X = betas,
  nX = ncol(betas),
  yes = yesarray,
  total = totalarray,
  C = sp_totes,
  maxes = sp_unique[,2],
  p_loss = ch0$p_tag_loss
)

# . Initial Values ----
z.init <- matrix(NA, nrow(ch_mod), ncol(ch_mod))
for (i in 1:nrow(ch_mod)){
  if(first[i] < ncol(z.init)){
    z.init[i, (first[i] + 1):ncol(z.init)] <- 1
  }
}

yes.init <- yes

inits <- function(){list(
  z = z.init,
  N = matrix(1000, nrow = nrow(sp_totes), ncol = ncol(sp_totes))
  )}  

# . Parameters to monitor ----
pars <- c('beta', 'gamma', 'phi', 'N')

# . MCMC settings ----
nc <- 3
nt <- 10
ni <- 1500
nb <- 500

# . Call JAGS and run model ----
rdt <- jags(dat, inits=inits, pars,
            model.file = textConnection(modelString), 
            n.chains = nc, n.iter = ni,
            n.burnin = nb, n.thin = nt,
            refresh = 0.001, DIC=TRUE
            )

# . Print model summary ----
print(rdt, digits=3)

# Result ----
# . Extract posteriors ----
sims <- rdt$BUGSoutput$sims.list

# . Detection ----
# .. Detection by primary period ----
# Detection for each MCMC iteration in 
# each primary period from posts
p_t <- boot::inv.logit(apply(sims$beta, c(1,2), mean))

# Melt data across remainig dimensions
# of 1: iteration, 2: primary period
p_t <- reshape2::melt(p_t)
names(p_t) <- c('iteration', 'primary_period', 'p')

# Get mean and 95% HDI of detection by primary period
p_t_summary <- plyr::ddply(.data = p_t,
                           .variables = 'primary_period',
                           plyr::summarize,
                           est = mean(p), 
                           lwr = quantile(p, probs=0.025),
                           upr = quantile(p, probs=0.975)
                           )

boxplot(p~primary_period, data=p_t,
        boxwex=.25, outline=FALSE, medlwd=1, col='gray87',
        names = c(1,2,3,4,5), ylab = 'Detection probability (p)',
        whisklty=1, yaxt='n', xlab='Primary period')
axis(side=2, las=TRUE)

# .. Detection by secondary period within primaries ----
# Mean detection in each secondary period within primary
# periods.
p_j <- boot::inv.logit(apply(sims$beta, c(1,2,3), mean))

# Melt into long-form data
p_j <- reshape2::melt(p_j)
names(p_j) <- c('iteration', 'primary_period', 'secondary_period', 'p')

# Get mean and 95% HDI of detection by primary period
p_j_summary <- plyr::ddply(.data = p_j,
                           .variables = c('primary_period', 'secondary_period'),
                           plyr::summarize,
                           est = mean(p), 
                           lwr = quantile(p, probs=0.025),
                           upr = quantile(p, probs=0.975)
                           )
# Make a boxplot
par(mar = c(5,5,1,1))
boxplot(p~primary_period+secondary_period, data=p_j,
        boxwex=.25, outline=FALSE, medlwd=1, col='gray87',
        names = rep(c(1,2,3), 5), ylab = 'Detection probability (p)',
        whisklty=1, yaxt='n', xlab='Secondary period')
abline(v=seq(3.5,3.5*4, 3))
axis(side = 2, las = TRUE)
axis(side = 3, at=seq(2, 2*7, 3), labels=c(1, 2, 3, 4, 5), tick=FALSE)
mtext(side = 3, line= 4, text = 'Primary period')

# .. Detection by tag configuration ----
# Summarize posteriors by tag configuration (dimension 4 of beta)
p_config <- boot::inv.logit(apply(sims$beta, c(1,4), mean))

# Make a plot
boxplot(p_config, ylim=c(0.35,.55), ylab = 'Detection probability (p)',
        boxwex=.25, outline=FALSE, medlwd=1, col='gray87',
        names = c("1F0P", "1F1P", "2F0P", "2F1P"),
        whisklty=1, yaxt='n', xlab='Tagging configuration')
axis(side=2, las=TRUE)

# . Abundance, phi, gamma checks
# .. Abundance by primary period ----
# Melt data across remainig dimensions
# of 1: iteration, 2: species, 3: primary period
n_t <- reshape2::melt(sims$N)
names(n_t) <- c('iteration', 'species', 'primary_period', 'N')

phi_t <- reshape2::melt(sims$phi)
names(phi_t) <- c('iteration', 'species', 'primary_period', 'Phi')

gamma_t <- reshape2::melt(sims$gamma)
names(gamma_t) <- c('iteration', 'species', 'primary_period', 'gamma')

species = 1

par(mfrow=c(3,1))
boxplot(N~primary_period, data=n_t[n_t$species==species,],
        boxwex=.25, outline=FALSE, medlwd=1, col='gray87',
        names = c(1,2,3,4,5), ylab = 'Abundance (N)',
        whisklty=1, yaxt='n', xlab='Primary period')
axis(side=2, las=TRUE)

boxplot(Phi~primary_period, data=phi_t[phi_t$species==species,],
        boxwex=.25, outline=FALSE, medlwd=1, col='gray87',
        names = c(1,2,3,4), ylab = 'Apparent survival',
        whisklty=1, yaxt='n', xlab='Primary period')
axis(side=2, las=TRUE)

boxplot(gamma~primary_period, data=gamma_t[gamma_t$species==species,],
        boxwex=.25, outline=FALSE, medlwd=1, col='gray87',
        names = c(1,2,3,4), ylab = 'Site fidelity',
        whisklty=1, yaxt='n', xlab='Primary period')
axis(side=2, las=TRUE)

