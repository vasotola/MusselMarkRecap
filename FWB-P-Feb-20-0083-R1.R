# Package load ----
library(tidyr)
library(dplyr)
library(reshape)
library(R2jags)
library(lubridate)

# Data manipulation ----
# . Tagging data ----
# Read in the tagging data
tags <- read.csv('data_tags.csv', stringsAsFactors = FALSE)
check <- vector(mode = 'character', length=nrow(tags))

# Replace slashes in floy_id with letters
for(i in 1:nrow(tags)){
  tags$floy_id[i] <- gsub(
    pattern = "/", 
    x=tags$floy_id[i],
    replacement=substr(tags$floy_id[i], start=1, stop=1),
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
      replacement=substr(rf_long$floy_id[i], start=1, stop=1),
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
      replacement=substr(rp_long$floy_id[i], start=1, stop=1),
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
    ch_test$primary[i] <- paste0("0", ch_test$primary[i])
  }
}

# Cross-tab into a 3-d array with
sp_totes <- xtabs(
  formula = ch ~ species + primary, data = ch_test
)


# Model specification ----
modelString = "
model {
  # Priors
    # Survival and emigration
      for(i in 1:(n.ind)){
        for(t in 1:(n.years-1)){
          # Logit-scale predictor of apparent survival (species)
            logit(r_phi[i, t]) <- l_phi[i, t]
            l_phi[i, t] <-  phi[species[i], t]
          
          # Logit-scale predictor of site fidelity (species)
            logit(r_gamma[i, t]) <- l_gamma[i, t]
            l_gamma[i, t] <- gamma[species[i], t]
        }
      }
    
    # Species- and time-specific priors on
    # phi and gamma
      for(s in 1:n.species){
        for(t in 1:(n.years-1)){
          phi[s, t] ~ dunif(0, 1)
          gamma[s, t] ~ dunif(0, 1)
        }
      }
    
    # Secondary occasions p
      # Logit-scale linear predictor on
      # p by year and secondary occasion
        for (i in 1:n.ind){
          for (t in 1:n.years){
            for (j in 1:(n.sec[t])){
                logit(p[i,t,j]) <- lp[i,t,j]  
                lp[i,t,j] <- alpha[t] + beta[1]*config[i,1] + beta[2]*config[i,2] +
                             beta[3]*config[i,3] + beta[4]*config[i,4]
                yes[i,t,j] ~ dbin(p[i,t,j], total[i,t,j])
            }
          }
        }
      # Alphas (intercepts) for linear predictor of logit-scale
      # p, varying by time-period
      for(t in 1:n.years){
        alpha[t] ~ dunif(0, 1)
      }
      
      # Betas for linear predictor of
      # logit-scale p
        for(c in 1:n.config){
          beta[c] ~ dnorm(0, 1)
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
    for(t in 1:n.years){
      C[s,t] ~ dbinom(mean(pstar[,t]), N[s,t])
      N[s,t] ~ dpois(lambda[s,t])
      log(lambda[s,t]) <- mu[s,t]
      mu[s,t] ~ dnorm(0, 0.001)
    }
  }
  
  # Derived quantities
  # Detection probability from year-specific intercepts (alpha) and
  # tag-specific intercept adjustments (beta)
  for(t in 1:n.years){
    pl.out[1 ,t] <- alpha[t] + beta[1]
  }
  for(c in 2:n.config){
    for(t in 1:n.years){
      pl.out[c, t] <- alpha[t] + beta[1] + beta[c]
    }
  }
  for(c in 1:n.config){
    for(t in 1:n.years){
      logit(p_est[c,t]) <- pl.out[c,t]
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
config <- model.matrix(glm(X2017.04.28~tag_config, data = ch0))
attributes(config)[2:3] <- NULL

# Package into a list for JAGS
dat <- list(
  first=first,
  ch = ch_mod,
  n.sec = rep(n.sec.occ, n.years),
  n.years = ncol(ch_mod),
  n.ind = nrow(ch_mod),
  n.config = length(unique(ch0$tag_config)),
  n.species = length(unique(ch0$species)),
  species = as.numeric(as.factor(ch0$species)),
  config = config,
  yes = yesarray,
  total = totalarray,
  C = sp_totes
)

# . Initial Values ----
z.init <- matrix(NA, nrow(ch_mod), ncol(ch_mod))
for (i in 1:nrow(ch_mod)){
  if(first[i] < ncol(z.init)){
    z.init[i, (first[i] + 1):ncol(z.init)] <- 1
  }
}

yes.init <- array(1, dim=unlist(dim(yesarray)), dimnames = NULL)
inits <- function(){list(
  z = z.init,
  N = matrix(1000, nrow = nrow(sp_totes), ncol = ncol(sp_totes))
  )}  

# . Parameters to monitor ----
pars <- c('alpha', 'beta', 'gamma', 'phi', 'N', 'p_est')

# . MCMC settings ----
nc <- 3
nt <- 10
ni <- 15000
nb <- 5000

# . Call JAGS and run model ----
rdt <- jags(dat, inits=inits, pars,
            model.file = textConnection(modelString), 
            n.chains = nc, n.iter = ni,
            n.burnin = nb, n.thin = nt,
            refresh = 0.001, DIC=TRUE
)

# . Print model summary ----
print(rdt, digits=3)









