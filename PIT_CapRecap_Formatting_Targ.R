setwd("C:/Users/vasot/Dropbox/Research/Texas_State/TXST_Research/Dissertation_Data/Mussels/Ricke_EnvModel/UpdatedModel_Review1/DataFormatting/Targ_PIT")

####### caphist Lower (Altair) PIT Tag Data ############
library(reshape)

pitreader_apr2017<-read.csv("Altair_April2017_ReaderFloy.csv",header=T)
pitreader_Aug2017<-read.csv("Altair_Aug2017_ReaderFloy.csv",header=T)
pitreader_Nov2017<-read.csv("Altair_Nov2017_ReaderFloy.csv",header=T)
pitreader_apr2018<-read.csv("Altair_April2018_ReaderFloy.csv",header=T)
pitreader_Aug2018<-read.csv("Altair_Aug2018_ReaderFloy.csv",header=T)


all_pitdata <- rbind(pitreader_apr2017,pitreader_Aug2017,pitreader_Nov2017,pitreader_apr2018,pitreader_Aug2018)
str(all_pitdata)

ALT_floy_caphist <- cast(all_pitdata,floy_id~date)
spp <- as.data.frame(all_pitdata[,c(9:10)])
spp[1:10,1:2]
str(spp)

spp <- subset(spp,!duplicated(spp$floy_id))


ALT_floy_caphist_spp <- merge(spp,ALT_floy_caphist[,c("floy_id",setdiff(colnames(ALT_floy_caphist),colnames(spp)))],by="floy_id")
dim(ALT_floy_caphist_spp)
ALT_floy_caphist_spp[1:10,1:17]
str(ALT_floy_caphist_spp)

caphist <- ALT_floy_caphist_spp[,c(3:17)]
ALT_floy_caphist_spp <- ALT_floy_caphist_spp[,-c(3:17)]

dim(caphist)
dim(ALT_floy_caphist_spp)

caphist[1:10,1:10]
caphist[caphist > "1"] <- 1
caphist[1:10,1:10]


ALT_floy_caphist_spp <- cbind(ALT_floy_caphist_spp,caphist)


write.csv(ALT_floy_caphist_spp,file="LWR_PIT_CH_Targ.csv")

caphist <- read.csv("LWR_PIT_CH_Targ.csv",header=T)
levels(caphist$species)
cypu <- subset(caphist,caphist$species=="QUHO")
cype <- subset(caphist,caphist$species=="QUPE")
#trma <- subset(caphist,caphist$species=="TRMA")

write.csv(cypu,file="LWR_CYPU_PIT_CH.csv")
write.csv(cype,file="LWR_CYPE_PIT_CH.csv")
#write.csv(trma,file="TRMA_flow_caphist.csv")


######## caphist Upper (San Saba) PIT ############
library(reshape)

pitreader_Apr2018<-read.csv("SanSaba_April2018_ReaderFloy.csv",header=T)
pitreader_April2019<-read.csv("SanSaba_April2019_ReaderFloy.csv",header=T)
pitreader_Aug2017<-read.csv("SanSaba_Aug2017_ReaderFloy.csv",header=T)
pitreader_Aug2018<-read.csv("SanSaba_Aug2018_ReaderFloy.csv",header=T)
pitreader_Nov2017<-read.csv("SanSaba_Nov2017_ReaderFloy.csv",header=T)

all_pitdata <- rbind(pitreader_Aug2017,pitreader_Nov2017,pitreader_Apr2018,pitreader_Aug2018,pitreader_April2019)
str(all_pitdata)
SS_floy_caphist <- cast(all_pitdata,floy_id~date)
spp <- as.data.frame(all_pitdata[,c(9:10)])
spp[1:10,1:2]
str(spp)

spp <- subset(spp,!duplicated(spp$floy_id))

SS_floy_caphist_spp <- merge(spp,SS_floy_caphist[,c("floy_id",setdiff(colnames(SS_floy_caphist),colnames(spp)))],by="floy_id")
dim(SS_floy_caphist_spp)
SS_floy_caphist_spp[1:10,1:15]
str(SS_floy_caphist_spp)

caphist <- SS_floy_caphist_spp[,c(3:17)]
SS_floy_caphist_spp <- SS_floy_caphist_spp[,-c(3:17)]

head(caphist)

caphist[caphist > "1"] <- 1
SS_floy_caphist_spp <- cbind(SS_floy_caphist_spp,caphist)

write.csv(SS_floy_caphist_spp,file="UPR_PIT_CH_Targ.csv")

caphist <- read.csv("UPR_PIT_CH_Targ.csv",header=T)
levels(caphist$species)

cypu <- subset(caphist,caphist$species=="QUHO")
cype <- subset(caphist,caphist$species=="QUPE")
#trma <- subset(caphist,caphist$species=="TRMA")

write.csv(cypu,file="UPR_CYPU_PIT_CH.csv")
write.csv(cype,file="UPR_CYPE_PIT_CH.csv")

