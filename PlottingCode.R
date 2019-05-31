setwd("C:/Users/vasot/Dropbox/Research/Texas_State/TXST_Research/Dissertation_Data/Mussels/Ricke_EnvModel/PitData")

load("FinalCapRecap_PIT_ADJ.RData")

#nc <- 3
#nt <- 10
#ni <- 15000
#nb <- 5000

############### PIT Tags #####################

####  SS
phi.qupe <- rdt.ss.qupe.adj$BUGSoutput$sims.list$phi
pstar.qupe <- rdt.ss.qupe.adj$BUGSoutput$sims.list$pstar
abund.qupe <- rdt.ss.qupe.adj$BUGSoutput$sims.list$Nin

phi.quho <- rdt.ss.quho.adj$BUGSoutput$sims.list$phi
pstar.quho <- rdt.ss.quho.adj$BUGSoutput$sims.list$pstar
abund.quho <- rdt.ss.quho.adj$BUGSoutput$sims.list$Nin


### ALT
phi.qupe <- rdt.alt.qupe.adj$BUGSoutput$sims.list$phi
pstar.qupe <- rdt.alt.qupe.adj$BUGSoutput$sims.list$pstar
abund.qupe <- rdt.alt.qupe.adj$BUGSoutput$sims.list$Nin

phi.quho <- rdt.alt.quho.adj$BUGSoutput$sims.list$phi
pstar.quho <- rdt.alt.quho.adj$BUGSoutput$sims.list$pstar
abund.quho <- rdt.alt.quho.adj$BUGSoutput$sims.list$Nin


## Figure 2 (Lower petrina and pustulosa estimates)

jpeg("Lower_PIT_Figure2.jpg",width=8,height=8,units="in",res=1000,family="serif")
par(mfrow=c(3,2))


## detection prob
# QUPE
par(mar=c(4,5,2,1))

# False boxplot to get access to
# statistics matrix
h.qupe <- boxplot(pstar.qupe, plot=FALSE)

# Overwrite whisker stats with
# 95% CRIs for posteriors
h.qupe$stats[c(1,5),] <- 
  apply(pstar.qupe, 2, quantile,
        probs=c(0.025, 0.975)
  )[c(1,2),]

# Make the plot with 95% CRIs
bxp(h.qupe, xaxt='n', yaxt='n',
    ylim=c(0,1),
    outline=FALSE,
    whiskcol='gray60', whisklwd=2, whisklty=1, 
    staplelty=0, staplecol='gray60',
    boxcol='gray60', boxlwd=2, boxfill='gray87',
    boxwex=.33,
    medcol='gray60'
)

# X (1) and Y (2) axes
axis(1)
axis(2, las=2)

# X and Y axis labels
mtext(side = 1, 'Primary period (t)', line=2.5)
mtext(side = 2, line=3,
      expression(paste(
        'Detection probability (p','*'['t'],')'))
)    
mtext(side = 3,at = c(.5), 'A')
mtext(side=3,expression(italic("C. petrina")))

# QUHO
par(mar=c(4,5,2,1))

# False boxplot to get access to
# statistics matrix
h.quho <- boxplot(pstar.quho, plot=FALSE)

# Overwrite whisker stats with
# 95% CRIs for posteriors
h.quho$stats[c(1,5),] <- 
  apply(pstar.quho, 2, quantile,
        probs=c(0.025, 0.975)
  )[c(1,2),]

# Make the plot with 95% CRIs
bxp(h.quho, xaxt='n', yaxt='n',
    ylim=c(0,1),
    outline=FALSE,
    whiskcol='gray60', whisklwd=2, whisklty=1, 
    staplelty=0, staplecol='gray60',
    boxcol='gray60', boxlwd=2, boxfill='gray87',
    boxwex=.33,
    medcol='gray60'
)

# X (1) and Y (2) axes
axis(1)
axis(2, las=2)

# X and Y axis labels
mtext(side = 1, 'Primary period (t)', line=2.5)
#mtext(side = 2, line=3,
 #     expression(paste('p','*'['t']))
#)  
mtext(side = 3,at = c(.5), 'B')
mtext(side=3,expression(italic("C. pustulosa")))

## estimated abundance
# QUPE
par(mar=c(4,5,1,1))

# False boxplot to get access to
# statistics matrix
h.qupe <- boxplot(abund.qupe, plot=FALSE)

# Overwrite whisker stats with
# 95% CRIs for posteriors
h.qupe$stats[c(1,5),] <- 
  apply(abund.qupe, 2, quantile,
        probs=c(0.025, 0.975)
  )[c(1,2),]

# Make the plot with 95% CRIs
bxp(h.qupe, xaxt='n', yaxt='n',
    #alt 
    ylim=c(0,125),
    #ss
    #ylim=c(0,400),
    outline=FALSE,
    whiskcol='gray60', whisklwd=2, whisklty=1, 
    staplelty=0, staplecol='gray60',
    boxcol='gray60', boxlwd=2, boxfill='gray87',
    boxwex=.33,
    medcol='gray60'
)

# X (1) and Y (2) axes
axis(1)
axis(2, las=2)

# X and Y axis labels
mtext(side = 1, 'Primary period (t)', line=2.5)
mtext(side = 2, line=3,'Estimated abundance (N)')    
mtext(side = 3,at = c(.5), 'C')

# QUHO
par(mar=c(4,5,1,1))

# False boxplot to get access to
# statistics matrix
h.quho <- boxplot(abund.quho, plot=FALSE)

# Overwrite whisker stats with
# 95% CRIs for posteriors
h.quho$stats[c(1,5),] <- 
  apply(abund.quho, 2, quantile,
        probs=c(0.025, 0.975)
  )[c(1,2),]

# Make the plot with 95% CRIs
bxp(h.quho, xaxt='n', yaxt='n',
    #alt
    ylim=c(0,300),
    #ss
    #ylim=c(0,20),
    outline=FALSE,
    whiskcol='gray60', whisklwd=2, whisklty=1, 
    staplelty=0, staplecol='gray60',
    boxcol='gray60', boxlwd=2, boxfill='gray87',
    boxwex=.33,
    medcol='gray60'
)

# X (1) and Y (2) axes
axis(1)
axis(2, las=2)

# X and Y axis labels
mtext(side = 1, 'Primary period (t)', line=2.5)
#mtext(side = 2, line=3,'Estimated Abundance (N)')     
mtext(side = 3,at = c(.5), 'D')


## apparent survival (phi)
# QUPE
par(mar=c(4,5,1,1))

# False boxplot to get access to
# statistics matrix
h.qupe <- boxplot(phi.qupe, plot=FALSE)

# Overwrite whisker stats with
# 95% CRIs for posteriors
h.qupe$stats[c(1,5),] <- 
  apply(phi.qupe, 2, quantile,
        probs=c(0.025, 0.975)
  )[c(1,2),]

# Make the plot with 95% CRIs
bxp(h.qupe, xaxt='n', yaxt='n',
    ylim=c(0,1),
    outline=FALSE,
    whiskcol='gray60', whisklwd=2, whisklty=1, 
    staplelty=0, staplecol='gray60',
    boxcol='gray60', boxlwd=2, boxfill='gray87',
    boxwex=.33,
    medcol='gray60'
)

# X (1) and Y (2) axes
axis(1)
axis(2, las=2)

# X and Y axis labels
mtext(side = 1, 'Interval (t to t+1)', line=2.5)
mtext(side = 2, line=3,
      expression(paste(
        'Apparent survival ( ', phi['t'], ' )'
      ))
)
mtext(side = 3,at = c(.5), 'E')

# QUHO
par(mar=c(4,5,1,1))

# False boxplot to get access to
# statistics matrix
h.quho <- boxplot(phi.quho, plot=FALSE)

# Overwrite whisker stats with
# 95% CRIs for posteriors
h.quho$stats[c(1,5),] <- 
  apply(phi.quho, 2, quantile,
        probs=c(0.025, 0.975)
  )[c(1,2),]

# Make the plot with 95% CRIs
bxp(h.quho, xaxt='n', yaxt='n',
    ylim=c(0,1),
    outline=FALSE,
    whiskcol='gray60', whisklwd=2, whisklty=1, 
    staplelty=0, staplecol='gray60',
    boxcol='gray60', boxlwd=2, boxfill='gray87',
    boxwex=.33,
    medcol='gray60'
)

# X (1) and Y (2) axes
axis(1)
axis(2, las=2)

# X and Y axis labels
mtext(side = 1, 'Interval (t to t+1)', line=2.5)
mtext(side = 3,at = c(.5), 'F')

#mtext(side = 2, line=3,
#      expression(paste(
#        'Apparent survival ( ', phi['t'], ' )'
#      ))
#)

dev.off()


################ GRUBBING TARGET #####################

setwd("C:/Users/vasot/Dropbox/Research/Texas_State/TXST_Research/Dissertation_Data/Mussels/Ricke_EnvModel/GrubbingTarg")

load("FinalCapRecap_Grub_ADJ.RData")

####  SS
#phi.qupe <- rdt.ss.qupe$BUGSoutput$sims.list$phi
#pstar.qupe <- rdt.ss.qupe$BUGSoutput$sims.list$pstar
#abund.qupe <- rdt.ss.qupe$BUGSoutput$sims.list$Nin

#phi.quho <- rdt.ss.quho$BUGSoutput$sims.list$phi
#pstar.quho <- rdt.ss.quho$BUGSoutput$sims.list$pstar
#abund.quho <- rdt.ss.quho$BUGSoutput$sims.list$Nin


### ALT
phi.qupe <- rdt.alt.qupe.adj$BUGSoutput$sims.list$phi
pstar.qupe <- rdt.alt.qupe.adj$BUGSoutput$sims.list$pstar
abund.qupe <- rdt.alt.qupe.adj$BUGSoutput$sims.list$Nin

phi.quho <- rdt.alt.quho.adj$BUGSoutput$sims.list$phi
pstar.quho <- rdt.alt.quho.adj$BUGSoutput$sims.list$pstar
abund.quho <- rdt.alt.quho.adj$BUGSoutput$sims.list$Nin


## Figure 2 (Lower petrina and pustulosa estimates)

jpeg("Lower_Grub_Figure4.jpg",width=8,height=8,units="in",res=1000,family="serif")
par(mfrow=c(3,2))


## detection prob
# QUPE
par(mar=c(4,5,2,1))

# False boxplot to get access to
# statistics matrix
h.qupe <- boxplot(pstar.qupe, plot=FALSE)

# Overwrite whisker stats with
# 95% CRIs for posteriors
h.qupe$stats[c(1,5),] <- 
  apply(pstar.qupe, 2, quantile,
        probs=c(0.025, 0.975)
  )[c(1,2),]

# Make the plot with 95% CRIs
bxp(h.qupe, xaxt='n', yaxt='n',
    ylim=c(0,1),
    outline=FALSE,
    whiskcol='gray60', whisklwd=2, whisklty=1, 
    staplelty=0, staplecol='gray60',
    boxcol='gray60', boxlwd=2, boxfill='gray87',
    boxwex=.33,
    medcol='gray60'
)

# X (1) and Y (2) axes
axis(1)
axis(2, las=2)

# X and Y axis labels
mtext(side = 1, 'Primary period (t)', line=2.5)
mtext(side = 2, line=3,
      expression(paste(
        'Detection probability (p','*'['t'],')'))
)    
mtext(side = 3,at = c(.5), 'A')
mtext(side=3,expression(italic("C. petrina")))

# QUHO
par(mar=c(4,5,2,1))

# False boxplot to get access to
# statistics matrix
h.quho <- boxplot(pstar.quho, plot=FALSE)

# Overwrite whisker stats with
# 95% CRIs for posteriors
h.quho$stats[c(1,5),] <- 
  apply(pstar.quho, 2, quantile,
        probs=c(0.025, 0.975)
  )[c(1,2),]

# Make the plot with 95% CRIs
bxp(h.quho, xaxt='n', yaxt='n',
    ylim=c(0,1),
    outline=FALSE,
    whiskcol='gray60', whisklwd=2, whisklty=1, 
    staplelty=0, staplecol='gray60',
    boxcol='gray60', boxlwd=2, boxfill='gray87',
    boxwex=.33,
    medcol='gray60'
)

# X (1) and Y (2) axes
axis(1)
axis(2, las=2)

# X and Y axis labels
mtext(side = 1, 'Primary period (t)', line=2.5)
#mtext(side = 2, line=3,
#     expression(paste('p','*'['t']))
#)  
mtext(side = 3,at = c(.5), 'B')
mtext(side=3,expression(italic("C. pustulosa")))

## estimated abundance
# QUPE
par(mar=c(4,5,1,1))

# False boxplot to get access to
# statistics matrix
h.qupe <- boxplot(abund.qupe, plot=FALSE)

# Overwrite whisker stats with
# 95% CRIs for posteriors
h.qupe$stats[c(1,5),] <- 
  apply(abund.qupe, 2, quantile,
        probs=c(0.025, 0.975)
  )[c(1,2),]

# Make the plot with 95% CRIs
bxp(h.qupe, xaxt='n', yaxt='n',
    #alt 
    ylim=c(0,100),
    #ss
    #ylim=c(0,400),
    outline=FALSE,
    whiskcol='gray60', whisklwd=2, whisklty=1, 
    staplelty=0, staplecol='gray60',
    boxcol='gray60', boxlwd=2, boxfill='gray87',
    boxwex=.33,
    medcol='gray60'
)

# X (1) and Y (2) axes
axis(1)
axis(2, las=2)

# X and Y axis labels
mtext(side = 1, 'Primary period (t)', line=2.5)
mtext(side = 2, line=3,'Estimated abundance (N)')    
mtext(side = 3,at = c(.5), 'C')

# QUHO
par(mar=c(4,5,1,1))

# False boxplot to get access to
# statistics matrix
h.quho <- boxplot(abund.quho, plot=FALSE)

# Overwrite whisker stats with
# 95% CRIs for posteriors
h.quho$stats[c(1,5),] <- 
  apply(abund.quho, 2, quantile,
        probs=c(0.025, 0.975)
  )[c(1,2),]

# Make the plot with 95% CRIs
bxp(h.quho, xaxt='n', yaxt='n',
    #alt
    ylim=c(0,325),
    #ss
    #ylim=c(0,20),
    outline=FALSE,
    whiskcol='gray60', whisklwd=2, whisklty=1, 
    staplelty=0, staplecol='gray60',
    boxcol='gray60', boxlwd=2, boxfill='gray87',
    boxwex=.33,
    medcol='gray60'
)

# X (1) and Y (2) axes
axis(1)
axis(2, las=2)

# X and Y axis labels
mtext(side = 1, 'Primary period (t)', line=2.5)
#mtext(side = 2, line=3,'Estimated Abundance (N)')     
mtext(side = 3,at = c(.5), 'D')


## apparent survival (phi)
# QUPE
par(mar=c(4,5,1,1))

# False boxplot to get access to
# statistics matrix
h.qupe <- boxplot(phi.qupe, plot=FALSE)

# Overwrite whisker stats with
# 95% CRIs for posteriors
h.qupe$stats[c(1,5),] <- 
  apply(phi.qupe, 2, quantile,
        probs=c(0.025, 0.975)
  )[c(1,2),]

# Make the plot with 95% CRIs
bxp(h.qupe, xaxt='n', yaxt='n',
    ylim=c(0,1),
    outline=FALSE,
    whiskcol='gray60', whisklwd=2, whisklty=1, 
    staplelty=0, staplecol='gray60',
    boxcol='gray60', boxlwd=2, boxfill='gray87',
    boxwex=.33,
    medcol='gray60'
)

# X (1) and Y (2) axes
axis(1)
axis(2, las=2)

# X and Y axis labels
mtext(side = 1, 'Interval (t to t+1)', line=2.5)
mtext(side = 2, line=3,
      expression(paste(
        'Apparent survival ( ', phi['t'], ' )'
      ))
)
mtext(side = 3,at = c(.5), 'E')

# QUHO
par(mar=c(4,5,1,1))

# False boxplot to get access to
# statistics matrix
h.quho <- boxplot(phi.quho, plot=FALSE)

# Overwrite whisker stats with
# 95% CRIs for posteriors
h.quho$stats[c(1,5),] <- 
  apply(phi.quho, 2, quantile,
        probs=c(0.025, 0.975)
  )[c(1,2),]

# Make the plot with 95% CRIs
bxp(h.quho, xaxt='n', yaxt='n',
    ylim=c(0,1),
    outline=FALSE,
    whiskcol='gray60', whisklwd=2, whisklty=1, 
    staplelty=0, staplecol='gray60',
    boxcol='gray60', boxlwd=2, boxfill='gray87',
    boxwex=.33,
    medcol='gray60'
)

# X (1) and Y (2) axes
axis(1)
axis(2, las=2)

# X and Y axis labels
mtext(side = 1, 'Interval (t to t+1)', line=2.5)
mtext(side = 3,at = c(.5), 'F')

#mtext(side = 2, line=3,
#      expression(paste(
#        'Apparent survival ( ', phi['t'], ' )'
#      ))
#)

dev.off()



#################### NON TARGET ####################

setwd("C:/Users/vasot/Dropbox/Research/Texas_State/TXST_Research/Dissertation_Data/Mussels/Ricke_EnvModel/NonTarget")


load("FinalCapRecap_NT_ADJ.RData")


## ALT AMPL
phi.ampl <- rdt.alt.ampl.adj$BUGSoutput$sims.list$phi
pstar.ampl <- rdt.alt.ampl.adj$BUGSoutput$sims.list$pstar
abund.ampl <- rdt.alt.ampl.adj$BUGSoutput$sims.list$Nin

## SS QUVE
phi.quve <- rdt.ss.quve.adj$BUGSoutput$sims.list$phi
pstar.quve <- rdt.ss.quve.adj$BUGSoutput$sims.list$pstar
abund.quve <- rdt.ss.quve.adj$BUGSoutput$sims.list$Nin


## Figure 2 (Lower petrina and pustulosa estimates)

jpeg("NonTarget_Grub_Figure6.jpg",width=8,height=8,units="in",res=1000,family="serif")
par(mfrow=c(3,2))


## detection prob
# AMPL
par(mar=c(4,5,2,1))

# False boxplot to get access to
# statistics matrix
h.ampl <- boxplot(pstar.ampl, plot=FALSE)

# Overwrite whisker stats with
# 95% CRIs for posteriors
h.ampl$stats[c(1,5),] <- 
  apply(pstar.ampl, 2, quantile,
        probs=c(0.025, 0.975)
  )[c(1,2),]

# Make the plot with 95% CRIs
bxp(h.ampl, xaxt='n', yaxt='n',
    ylim=c(0,1),
    outline=FALSE,
    whiskcol='gray60', whisklwd=2, whisklty=1, 
    staplelty=0, staplecol='gray60',
    boxcol='gray60', boxlwd=2, boxfill='gray87',
    boxwex=.33,
    medcol='gray60'
)

# X (1) and Y (2) axes
axis(1)
axis(2, las=2)

# X and Y axis labels
mtext(side = 1, 'Primary period (t)', line=2.5)
mtext(side = 2, line=3,
      expression(paste(
        'Detection probability (p','*'['t'],')'))
)    
mtext(side = 3,at = c(.5), 'A')
mtext(side=3,expression(italic("A. plicata")))

# QUVE
par(mar=c(4,5,2,1))

# False boxplot to get access to
# statistics matrix
h.quve <- boxplot(pstar.quve, plot=FALSE)

# Overwrite whisker stats with
# 95% CRIs for posteriors
h.quve$stats[c(1,5),] <- 
  apply(pstar.quve, 2, quantile,
        probs=c(0.025, 0.975)
  )[c(1,2),]

# Make the plot with 95% CRIs
bxp(h.quve, xaxt='n', yaxt='n',
    ylim=c(0,1),
    outline=FALSE,
    whiskcol='gray60', whisklwd=2, whisklty=1, 
    staplelty=0, staplecol='gray60',
    boxcol='gray60', boxlwd=2, boxfill='gray87',
    boxwex=.33,
    medcol='gray60'
)

# X (1) and Y (2) axes
axis(1)
axis(2, las=2)

# X and Y axis labels
mtext(side = 1, 'Primary period (t)', line=2.5)
#mtext(side = 2, line=3,
#     expression(paste('p','*'['t']))
#)  
mtext(side = 3,at = c(.5), 'B')
mtext(side=3,expression(italic("T. verrucosa")))

## estimated abundance
# QUPE
par(mar=c(4,5,1,1))

# False boxplot to get access to
# statistics matrix
h.ampl <- boxplot(abund.ampl, plot=FALSE)

# Overwrite whisker stats with
# 95% CRIs for posteriors
h.ampl$stats[c(1,5),] <- 
  apply(abund.ampl, 2, quantile,
        probs=c(0.025, 0.975)
  )[c(1,2),]

# Make the plot with 95% CRIs
bxp(h.ampl, xaxt='n', yaxt='n',
    #alt 
    ylim=c(0,300),
    outline=FALSE,
    whiskcol='gray60', whisklwd=2, whisklty=1, 
    staplelty=0, staplecol='gray60',
    boxcol='gray60', boxlwd=2, boxfill='gray87',
    boxwex=.33,
    medcol='gray60'
)

# X (1) and Y (2) axes
axis(1)
axis(2, las=2)

# X and Y axis labels
mtext(side = 1, 'Primary period (t)', line=2.5)
mtext(side = 2, line=3,'Estimated abundance (N)')    
mtext(side = 3,at = c(.5), 'C')

# quve
par(mar=c(4,5,1,1))

# False boxplot to get access to
# statistics matrix
h.quve <- boxplot(abund.quve, plot=FALSE)

# Overwrite whisker stats with
# 95% CRIs for posteriors
h.quve$stats[c(1,5),] <- 
  apply(abund.quve, 2, quantile,
        probs=c(0.025, 0.975)
  )[c(1,2),]

# Make the plot with 95% CRIs
bxp(h.quve, xaxt='n', yaxt='n',
    ylim=c(0,125),
    outline=FALSE,
    whiskcol='gray60', whisklwd=2, whisklty=1, 
    staplelty=0, staplecol='gray60',
    boxcol='gray60', boxlwd=2, boxfill='gray87',
    boxwex=.33,
    medcol='gray60'
)

# X (1) and Y (2) axes
axis(1)
axis(2, las=2)

# X and Y axis labels
mtext(side = 1, 'Primary period (t)', line=2.5)
#mtext(side = 2, line=3,'Estimated Abundance (N)')     
mtext(side = 3,at = c(.5), 'D')


## apparent survival (phi)
# QUPE
par(mar=c(4,5,1,1))

# False boxplot to get access to
# statistics matrix
h.ampl <- boxplot(phi.ampl, plot=FALSE)

# Overwrite whisker stats with
# 95% CRIs for posteriors
h.ampl$stats[c(1,5),] <- 
  apply(phi.ampl, 2, quantile,
        probs=c(0.025, 0.975)
  )[c(1,2),]

# Make the plot with 95% CRIs
bxp(h.ampl, xaxt='n', yaxt='n',
    ylim=c(0,1),
    outline=FALSE,
    whiskcol='gray60', whisklwd=2, whisklty=1, 
    staplelty=0, staplecol='gray60',
    boxcol='gray60', boxlwd=2, boxfill='gray87',
    boxwex=.33,
    medcol='gray60'
)

# X (1) and Y (2) axes
axis(1)
axis(2, las=2)

# X and Y axis labels
mtext(side = 1, 'Interval (t to t+1)', line=2.5)
mtext(side = 2, line=3,
      expression(paste(
        'Apparent survival ( ', phi['t'], ' )'
      ))
)
mtext(side = 3,at = c(.5), 'E')

# quve
par(mar=c(4,5,1,1))

# False boxplot to get access to
# statistics matrix
h.quve <- boxplot(phi.quve, plot=FALSE)

# Overwrite whisker stats with
# 95% CRIs for posteriors
h.quve$stats[c(1,5),] <- 
  apply(phi.quve, 2, quantile,
        probs=c(0.025, 0.975)
  )[c(1,2),]

# Make the plot with 95% CRIs
bxp(h.quve, xaxt='n', yaxt='n',
    ylim=c(0,1),
    outline=FALSE,
    whiskcol='gray60', whisklwd=2, whisklty=1, 
    staplelty=0, staplecol='gray60',
    boxcol='gray60', boxlwd=2, boxfill='gray87',
    boxwex=.33,
    medcol='gray60'
)

# X (1) and Y (2) axes
axis(1)
axis(2, las=2)

# X and Y axis labels
mtext(side = 1, 'Interval (t to t+1)', line=2.5)
mtext(side = 3,at = c(.5), 'F')

#mtext(side = 2, line=3,
#      expression(paste(
#        'Apparent survival ( ', phi['t'], ' )'
#      ))
#)

dev.off()
