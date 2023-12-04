#=========================================#
#          Phenotype for GWAS
#=========================================#

library(lme4)
library(dplyr)
library(plyr)
library(nephro)
library(ggplot2)

#-----------------------------------------------------#
#------------------- Serum Creatinine ----------------
#-----------------------------------------------------#

# Serumcreatinine Instrument controlled
chris$SerumCreatinine.Ins <- as.factor(chris$SerumCreatinine.Ins)
mSCr                      <- lm(SerumCreatinine ~ SerumCreatinine.Ins, data = chris)
myresidmSCr               <- data.frame(resm = mSCr$residuals)
chris$SerumCreatinine.Res <- NA #appending residuals
chris$SerumCreatinine.Res[as.numeric(rownames(myresidmSCr))] <- myresidmSCr$resm
chris$SerumCreatinine.Std <- 
  mean(chris$SerumCreatinine, na.rm = TRUE) + chris$SerumCreatinine.Res
#---------#

chris$Sex          <- as.factor(chris$Sex)
chris$week         <- as.factor(chris$week)
mSCr.Mixed         <- lmer(SerumCreatinine.Std ~ Age + Sex + (1|week), data = chris, na.action = "na.omit")#length(resid(mSCr.Mixed))
myresidmSCr.Mixed  <- data.frame(resm = resid(mSCr.Mixed))
chris$SerumCreatinine.Std.Mixed.Res <- NA
#---------#

# >>>>>>> 1st phenotype: This should go for GWAS <<<<<<<

chris$SerumCreatinine.Std.Mixed.Res[as.numeric(rownames(myresidmSCr.Mixed))] <- myresidmSCr.Mixed$resm
#Should not have used SCr.Std.Cleaned for GWAS
#chris$SerumCreatinine.Std.Cleaned <- mean(chris$SerumCreatinine.Std, na.rm = TRUE) + chris$SerumCreatinine.Std.Mixed.Res

#log(serum creatinine) ~ age + sex + (1|week)
chris$SerumCreatinine.Std.log           <- log(chris$SerumCreatinine.Std)
mSCrlog                                 <- lmer(SerumCreatinine.Std.log ~ Age + Sex + (1|week), data = chris, na.action = "na.omit")
myresidmSCrlog                          <- data.frame(resm = resid(mSCrlog))
chris$SerumCreatinine.Std.log.mixed.Res <- NA
chris$SerumCreatinine.Std.log.mixed.Res[as.numeric(rownames(myresidmSCrlog))] <- myresidmSCrlog$resm
chris$SerumCreatinine.Std.log.Cleaned   <- 
  mean(chris$SerumCreatinine.Std.log, na.rm = TRUE) + chris$SerumCreatinine.Std.log.mixed.Res


#-----------------------------------------------------#
#----------------------- eGFR ------------------------
#-----------------------------------------------------#

# Controlling Instruments effect on eGFR

table(chris$Sex) #SEX needs to be recoded
chris$Sex       <- dplyr::recode(chris$Sex, `1` = 1L, `2` = 0L) #1:Males 0:Females
#---------#

# Using CKDEpi.creat function in nephro package

chris$ethnicity <- 0 # 0 for Europeans, and 1 for Africans
chris$Sex       <- as.numeric(chris$Sex) #class(chris$Sex) is integer
#---------#

#eGFR adjusted by SerumCreat instrument
chris$eGFR      <- nephro::CKDEpi.creat(chris$SerumCreatinine.Std, chris$Sex, chris$Age, chris$ethnicity)
#---------#

#Mixed-model Regression
chris$Sex       <- as.factor(chris$Sex)
meGFR           <- lmer(eGFR ~ Age + Sex + (1|week), data = chris, na.action = "na.omit")
myresidmeGFR    <- data.frame(resm = resid(meGFR))
chris$eGFR.Res  <- NA
#---------#

#>>>>>>> 2nd phenotype: This should go for GWAS <<<<<<<
chris$eGFR.Res[as.numeric(rownames(myresidmeGFR))] <- myresidmeGFR$resm
#chris$eGFR.Std <- mean(chris$eGFR, na.rm = TRUE) + chris$eGFR.Res
#---------#

#log(eGFR)
chris$eGFR.log  <- log(chris$eGFR)   #I should not have used log(chris$eGFR.Std) as log(eGFR)
meGFRlog        <- lmer(chris$eGFR.log ~ Age + Sex + (1|week), data = chris, na.action = "na.omit")#length(resid(m4))
myresidmeGFRlog <- data.frame(resm = resid(meGFRlog))
chris$eGFR.log.Res <- NA
#>>>>>>> 3rd phenotype: This should go for GWAS <<<<<<<
chris$eGFR.log.Res[as.numeric(rownames(myresidmeGFRlog))] <- myresidmeGFRlog$resm
#chris$eGFR.log.Std <- mean(chris$eGFR.log, na.rm = TRUE) + chris$eGFR.log.Res



#-----------------------------------------------------#
#------------------------ UACR -----------------------
#-----------------------------------------------------#

# saving residuals
# appending residuals to chris dataset
#---------#

# Second way of controlling Instruments
#m2 <- lm(chris$UACR ~ chris$UACR.Ins)
#chris$UACR.Ins.Std <- m2$resid + mean(chris$UACR)
#---------#

# Third UACR computation
# This approach is no longer being used in our pipeline
# Look at UACR calculation section below of the scripts


#UACR
chris$Sex            <- as.factor(chris$Sex)
mUACR                <- lmer(UACR ~ Age + Sex + (1|week), data = chris, na.action = "na.omit")
myresidmUACR         <- data.frame(resm = resid(mUACR))
chris$UACR.Res       <- NA
chris$UACR.Res[as.numeric(rownames(myresidmUACR))] <- myresidmUACR$resm
chris$UACR.Std       <- mean(chris$UACR, na.rm = TRUE) + chris$UACR.Res
#---------#

#log(UACR)
mUACRlog             <- lmer(UACR.log ~ Age + Sex + (1|week), data = chris, na.action = "na.omit")
myresidmUACRlog      <- data.frame(resm = resid(mUACRlog))
chris$UACR.log.Res   <- NA
chris$UACR.log.Res[as.numeric(rownames(myresidmUACRlog))] <- myresidmUACRlog$resm
chris$UACR.log.Std   <- mean(chris$UACR.log, na.rm = TRUE) + chris$UACR.log.Res
#---------#

#log(UACR.log)
mUACR2xlog           <- lmer(UACR.2xlog ~ Age + Sex + (1|week), data = chris, na.action = "na.omit")
myresidmUACR2xlog    <- data.frame(resm = resid(mUACR2xlog))
chris$UACR.2xlog.Res <- NA
chris$UACR.2xlog.Res[as.numeric(rownames(myresidmUACR2xlog))] <- myresidmUACR2xlog$resm
chris$UACR.2xlog.Std <- mean(chris$UACR.2xlog, na.rm = TRUE) + chris$UACR.2xlog.Res
#---------#

summary(chris[, c("UACR", "UACR.Std", "UACR.log.Std", "UACR.2xlog.Std")])



#-----------------------------------------------------#
#--------------------- UACR-GFR ----------------------
#-----------------------------------------------------#


#UACR controlled Age, Sex, eGFR & (already Instruments)
chris$Sex          <- as.factor(chris$Sex)
m8                 <- lm(UACR ~ Age + Sex + eGFR, data = chris)
myresidm8          <- data.frame(resm = m8$residuals)
chris$UACR.GFR.Res <- NA
chris$UACR.GFR.Res[as.numeric(rownames(myresidm8))] <- myresidm8$resm
chris$UACR.GFR.Std <- mean(chris$UACR, na.rm = TRUE) + chris$UACR.GFR.Res
#---------#

# log(UACR)> ln(UACR-GFR)
# Natural log or ln equals to log() or logb(chris$UACR, base = exp(1))
chris$UACR.log         <- log(chris$UACR) 
m9                     <- lm(UACR.log ~ Age + Sex + eGFR, data = chris)
myresidm9              <- data.frame(resm = m9$residuals)
chris$UACR.GFR.log.Res <- NA
chris$UACR.GFR.log.Res[as.numeric(rownames(myresidm9))] <- myresidm9$resm
chris$UACR.GFR.log.Std <- mean(chris$UACR.log, na.rm = TRUE) + chris$UACR.GFR.log.Res
#---------#

#ln(ln(UACR)) > 2xln(UACR-GFR)
chris$UACR.2xlog         <- log(chris$UACR.log)
m10                      <- lm(UACR.2xlog ~ Age + Sex + eGFR, data = chris)
myresidm10 <- data.frame(resm = m10$residuals)
chris$UACR.GFR.2xlog.Res <- NA
chris$UACR.GFR.2xlog.Res[as.numeric(rownames(myresidm10))] <- myresidm10$resm
chris$UACR.GFR.2xlog.Std <- mean(chris$UACR.2xlog, na.rm = TRUE) + chris$UACR.GFR.2xlog.Res



#-----------------------------------------------------#
#-------------------- Proportions --------------------
#-----------------------------------------------------#

#Serum to Urine Creatinine
chris$Serum2Urine.Creat     <- chris$SCr.Ins.Std/chris$UrinaryCreatinine.Std
mS2UCreat                   <- lmer(Serum2Urine.Creat ~ Age + Sex+ (1|week), data = chris)
myresidmS2UCreat            <- data.frame(resm = mS2UCreat$residuals)
chris$S2UCreat.Res          <- NA
chris$S2UCreat.Res[as.numeric(rownames(myresidmS2UCreat))] <- myresidmS2UCreat$resm
chris$S2UCreat.Std          <- mean(chris$Serum2Urine.Creat, na.rm = TRUE) + chris$S2UCreat.Res
#---------#

#log(S2UCreat)
chris$Serum2Urine.Creat.log <- log(chris$Serum2Urine.Creat)
mS2UCreatlog                <- lmer(Serum2Urine.Creat.log ~ Age + Sex+ (1|week), data = chris)
myresidmS2UCreatlog         <- data.frame(resm = mS2UCreatlog$residuals)
chris$S2UCreat.log.Res      <- NA
chris$S2UCreat.log.Res[as.numeric(rownames(myresidmS2UCreatlog))] <- myresidmS2UCreatlog$resm
chris$S2UCreat.log.Std      <- mean(chris$Serum2Urine.Creat.log, na.rm = TRUE) + chris$S2UCreat.log.Res
#---------#

#Serum to Urine Albumin
chris$Serum2Urine.Album     <- chris$SerumAlbumin.Std/chris$UrinaryAlbumin.Std
mS2UAlbum                   <- lmer(Serum2Urine.Album ~ Age + Sex+ (1|week), data = chris)
myresidmS2UAlbum            <- data.frame(resm = mS2UAlbum$residuals)
chris$S2UAlbum.Res          <- NA
chris$S2UAlbum.Res[as.numeric(rownames(myresidmS2UAlbum))] <- myresidmS2UAlbum$resm
chris$S2UAlbum.Std          <- mean(chris$Serum2Urine.Album, na.rm = TRUE) + chris$S2UAlbum.Res
#---------#

#log(S2UAlbum)
chris$Serum2Urine.Album.log <- log(chris$Serum2Urine.Album)
mS2UAlbumlog                <- lmer(Serum2Urine.Album.log ~ Age + Sex+ (1|week), data = chris)
myresidmS2UAlbumlog         <- data.frame(resm = mS2UAlbumlog$residuals)
chris$S2UAlbum.log.Res      <- NA
chris$S2UAlbum.log.Res[as.numeric(rownames(myresidmS2UAlbumlog))] <- myresidmS2UAlbumlog$resm
chris$S2UAlbum.log.Std      <- mean(chris$Serum2Urine.Album.log, na.rm = TRUE) + chris$S2UAlbum.log.Res



#------------------------------------------------------#
#------------------- eGFR corrected --------------------
#------------------------------------------------------#


#chris$Sex     <- as.factor(chris$Sex)  # -> Later for the lm model
chris$week    <- as.factor(chris$week)
chris$SCr.Ins <- as.factor(chris$SCr.Ins)
#------------------------------------------------------#

#serum creatinine new scheme
mSCr.Std                      <- lmer(SCr ~ SCr.Ins + (1|week), data = chris, na.action = "na.exclude")
myresidmSCr.Std               <- data.frame(resm = resid(mSCr.Std))
chris$SerumCreatinine.Res.Std <- NA 
chris$SerumCreatinine.Res.Std[as.numeric(rownames(myresidmSCr.Std))] <- myresidmSCr.Std$resm
chris$SerumCreatinine.Std     <- mean(chris$SCr, na.rm = TRUE) + chris$SerumCreatinine.Res.Std


#------------------------------------------------------#
# Final scheme for computing eGFR
#------------------------------------------------------#

# Check if SEX needs to be re-coded ->
# In Nephro, males are coded in "1"

table(chris$Sex)

# Start computing eGFR with the 
# above standardized Serum Creatinine

chris$Sex        <- as.numeric(as.character(chris$Sex))
chris$eGFR       <- CKDEpi.creat(chris$SerumCreatinine.Std, chris$Sex, chris$Age, chris$ethnicity)
chris$eGFR.log   <- log(chris$eGFR)
chris$eGFRw      <- chris$eGFR
chris[chris$eGFR < 15 & is.na(chris$eGFR) != TRUE, "eGFRw"] <- 15
chris$eGFRw.log  <- log(chris$eGFRw)
#------------------------------------------------------#

#plots for comaprison

library(ggplot2)

ggplot(chris, aes(eGFRw.log.Res)) + 
  geom_density(color = "steelblue", fill = "gold2", alpha=0.9)+
  theme_minimal() #+ xlim(-10,100)
#geom_vline(aes(xintercept=mean(1.5, na.rm=TRUE)), color="#FC4E07", linetype="dashed", size=1)

ggplot(chris, aes(SerumCreatinine.Stdw.Res)) + 
  geom_density(color = "steelblue", fill = "turquoise3", alpha=0.9) +
  theme_classic()
ggsave("SerumCreatinine.Stdw.Res_cutpoint1.5_density.png",last_plot(), width = 8, height = 5.5, pointsize = 5, dpi = 300, units = "in")

grid.arrange(ggplot(chris, aes(UrinaryAlbumin.Std)) + 
               geom_density(color = "steelblue", fill = "gold2", alpha=0.9) + 
               theme_minimal() +
               geom_vline(aes(xintercept=mean(1.5, na.rm=TRUE)), color="#FC4E07", linetype="dashed", size=1),
             ggplot(chris, aes(UrinaryAlbumin)) + 
               geom_density(color = "steelblue", fill = "gold2", alpha=0.9) +
               theme_minimal() +
               geom_vline(aes(xintercept=mean(1.5, na.rm=TRUE)), color="#FC4E07", linetype="dashed", size=1),
             ncol=1)

ggplot(chris, aes(week, UrinaryAlbumin)) + geom_point()#+ ylim(0, 100)


#------------------------------------------------------#
#>>>>>>> These three variables should go for GWAS <<<<<<<
#------------------------------------------------------#


#01 - winsorized SerumCreatinine
chris$SerumCreatinine.Stdw     <- chris$SerumCreatinine.Std
chris[chris$SerumCreatinine.Std > 1.5 &
      !is.na(chris$SerumCreatinine.Std), "SerumCreatinine.Stdw"] <- 1.5
chris$Sex                      <- as.factor(chris$Sex)
mSCr.Stdw                      <- lm(SerumCreatinine.Stdw ~ Age + Sex, data = chris, na.action = "na.omit")
myresidmSCr.Stdw               <- data.frame(resm = mSCr.Stdw$residuals)
chris$SerumCreatinine.Stdw.Res <- NA
chris$SerumCreatinine.Stdw.Res[as.numeric(rownames(myresidmSCr.Stdw))] <- myresidmSCr.Stdw$resm
#------------------------------------------------------#

#02 - eGFR residuals by new scheme
chris$Sex           <- as.factor(chris$Sex)
meGFRw              <- lm(eGFRw ~ Age + Sex, data = chris)
myresidmeGFRw       <- data.frame(resm = meGFRw$residuals)
chris$eGFRw.Res     <- NA
chris$eGFRw.Res[as.numeric(rownames(myresidmeGFRw))] <- myresidmeGFRw$resm
#------------------------------------------------------#

#03 - eGFR.log residuals by new scheme
chris$Sex           <- as.factor(chris$Sex)
meGFRw.log          <- lm(eGFRw.log ~ Age + Sex, data = chris)
myresidmeGFRw.log   <- data.frame(resm = meGFRw.log$residuals)
chris$eGFRw.log.Res <- NA
chris$eGFRw.log.Res[as.numeric(rownames(myresidmeGFRw.log))] <- myresidmeGFRw.log$resm

#------------------------------------------------------#
#------------------- UACR calculation -----------------
#------------------------------------------------------#

# UACR without and with controlling instruments effect
# chris$UACR       <- (chris$UrinaryAlbumin.Std/10)/(chris$UrinaryCreatinine.Std/1000)
#---------#

# Note: UARC after in-place quantile normalization 
# of the traits (look at 00_data_transformation.R) 
# can be calculated as follows:

chris$UACR       <- (chris$UAlb/10)/(chris$UCr/1000)
chris$UACR.log   <- log(chris$UACR)
chris$UACR.2xlog <- log(chris$UACR.log)

chris$UACR.Ins <- apply(chris[,c("UrinaryAlbumin.Ins","UrinaryCreatinine.Ins")], 1, sum)
chris$eGFR.cat <- ifelse(chris$eGFR < 60, "CKD_eGFR", "NonCKD_eGFR")
chris$UACR.cat <- ifelse(chris$UACR > 30, "CKD_UACR", "NonCKD_UACR")

#table(chris$eGFR.cat, chris$UACR.cat)
#prop.table(table(chris$eGFR.cat, chris$UACR.cat))*100



#-----------------------------------------------------#
#---------------- Export phenotype file --------------
#-----------------------------------------------------#

# Have a look at the example .ped files here.
# PED file is nothing more than just a text file with the suffix ".ped".
# Its only difference with a .txt file is that a "#" is preceding the 
# column's names in order to be readable by EPACTS.
# 

# Example .ped files
testf <- read.delim("C://Users//Dariush_G3//Desktop//Coursera Courses//GWAS course//Eva_Example//CHRIS10K_example.ped", header=TRUE)
testf <- read.delim("C://Users//Dariush_G3//Desktop//Coursera Courses//GWAS course//ChF_Example//hemostasis.inv.ped")

table(testf$X.FAM_ID != testf$IND_ID)
table(CHRIS$AID != chris$AID)
#---------#

chris$Sex <- as.numeric(chris$Sex)

UACRdf <- 
  data.frame(FAM_ID = CHRIS$AID,
                   IND_ID = CHRIS$AID,
                   FAT_ID = rep(0,length(chris$AID)),
                   MOT_ID = rep(0,length(chris$AID)),
                   #SEX    = chris$Sex,
                   #AGE    = chris$Age,
                   #eGFR   = chris$eGFR,
                   #UACR.Ins=chris$UACR.Ins,
                   #UACR   = chris$UACR
                   llUACR.res = chris$llUACR.Res
)

write.table(UACRdf, file = "UACR.ped", sep = "\t", row.names = FALSE, quote=FALSE)
#---------#

eGFRdf <- 
  data.frame(FAM_ID = CHRIS$AID,
                   IND_ID = CHRIS$AID,
                   FAT_ID = rep(0,length(chris$AID)),
                   MOT_ID = rep(0,length(chris$AID)),
                   eGFR.Res = chris$eGFR.Res
)

write.table(eGFRdf, file = "eGFR.res.test.ped", sep = "\t", row.names = FALSE, quote=FALSE)
#---------#

SCrClndf <- 
  data.frame(FAM_ID = CHRIS$AID,
                     IND_ID = CHRIS$AID,
                     FAT_ID = rep(0,length(chris$AID)),
                     MOT_ID = rep(0,length(chris$AID)),
                     SCr.Cleaned = chris$SCr.cleaned, 
                     SCr.Ins.Week.Res = chris$SCr.Ins.Week.Res,
                     log.eGFR.Res = chris$leGFR.Res,
                     log.UACR.Res = chris$lUACR.Res,
                     UACR.Res = chris$UACR.Res
)

write.table(SCrClndf, file = "SCr.Cleaned4.ped", sep = "\t", row.names = FALSE, quote=FALSE)
#---------#

phenodf3 <- 
  data.frame(FAM_ID = chris$AID,
                     IND_ID = chris$AID,
                     FAT_ID = rep(0,nrow(chris)),
                     MOT_ID = rep(0,nrow(chris)),
                     SerumCreatinine.Cleaned   = chris$SerumCreatinine.Std.Cleaned,
                     SerumCreatinine.log.Cleaned  = chris$SerumCreatinine.Std.log.Cleaned
                     #SerumAlbumin.cleaned      = chris$SerumAlbumin.cleaned,
                     #UrinaryCreatinine.cleaned = chris$UrinaryCreatinine.cleaned,
                     #UrinaryAlbumin.cleaned    = chris$UrinaryAlbumin.cleaned,
                     #UACR.Std                  = chris$UACR.Std,
                     #UACR.log.Std              = chris$UACR.log.Std,
                     #UACR.2xlog.Std            = chris$UACR.2xlog.Std,
                     #eGFR.Std                  = chris$eGFR.Std,
                     #eGFR.log.Std              = chris$eGFR.log.Std,
                     #UACR.GFR.Std              = chris$UACR.GFR.Std,
                     #UACR.GFR.log.Std          = chris$UACR.GFR.log.Std,
                     #UACR.GFR.2xlog.Std        = chris$UACR.GFR.2xlog.Std,
                     #S2UCreat.Std              = chris$S2UCreat.Std,
                     #S2UCreat.log.Std          = chris$S2UCreat.log.Std,
                     #S2UAlbum.Std              = chris$S2UAlbum.Std,
                     #S2UAlbum.log.Std          = chris$S2UAlbum.log.Std
)

write.table(phenodf3, file = "phenodf3.ped", sep = "\t", row.names = FALSE, quote=FALSE)
#---------#
phenodf4 <- 
  data.frame(FAM_ID = chris$AID,
             IND_ID = chris$AID,
             FAT_ID = rep(0,nrow(chris)),
             MOT_ID = rep(0,nrow(chris)),
             SerumCreatinine.Stdw.Res = chris$SerumCreatinine.Stdw.Res,
             eGFRw.Res                = chris$eGFRw.Res,
             eGFRw.log.Res            = chris$eGFRw.log.Res)

write.table(phenodf4, file = "phenodf_NEW_scheme_W.ped", sep = "\t", row.names = FALSE, quote=FALSE)

#-----------------------------------------------------#
#-----------------------------------------------------#
#-----------------------------------------------------#




