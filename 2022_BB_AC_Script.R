#load libraries
library(lme4)
library(survival)
library(Hmisc)
library(ggplot2)
library(AICcmodavg)
library(rafalib)
library(MuMIn)
library(MASS)
library(pROC)
library(car)
library(caret)
library(psych)
library(dplyr)
library(tidyr)
library(broom)
library(tidyverse)
library(magrittr)
library(irr)
library(splitstackshape)
library(jtools)
library(ggstance)
library(ggh4x)
library(grid)
library(gridExtra)
library(jtools)
library(interactions)
library(BAMMtools)
library(DescTools)
library(AMR)
library(patchwork)
library(reshape2)
library(AER)
library(glmmTMB)
library(performance)
library(rsq)
library(ggpubr)
library(corrplot)
library(magick)
library(sjPlot)
library(ggbiplot)
library(glmm)
library(ggiraphExtra)
#########################################################################

#Analysis 1 : Does AC affect bear ORD and DD (i.e., does it affect wariness)

#Confirm SD is not a potential confound
SD <- read.csv("SD.csv")
SDCont <- subset(SD, ExptCont == "Control")
SDTreat <- subset(SD, ExptCont == "Treatment")

t.test(SDCont$SD ~ SDCont$PrePost) #t = -1.4308, df = 13.084, p-value = 0.1759
t.test(SDTreat$SD ~ SDTreat$PrePost) #t = 0.38182, df = 21.563, p-value = 0.7063


#Read in data
Dists <- read.csv("ORD_DD.csv")

#Create new columns that show percent difference in ORD and DD
Dists <- Dists %>%
  dplyr::mutate(PercDD = (PostDD - PreDD) / PreDD) %>%
  dplyr::mutate(PercORD = (PostORD - PreORD) / PreORD)

#Convert to percentages
Dists$PercDD <- Dists$PercDD*100
Dists$PercORD <- Dists$PercORD*100

#Determine count of individuals from each year
sum(Dists$Year == "2007") #11
sum(Dists$Year == "2008") #11

#Determine count of individuals from each sex
sum(Dists$Sex == "M") #12
sum(Dists$Sex == "F") #10

#Determine count of individuals from each sound group
sum(Dists$Sound == "Sound") #7
sum(Dists$Sound == "NoSound") #7
sum(Dists$Sound == "Control") #8

#Start with ORD

#Develop a Null model
ORD_null <- lm(PercORD~1, data = Dists)
summary(ORD_null)
AICc(ORD_null) #AICc = 248.0

#Add the treatment predictor to the model
ORD_T <- lm(PercORD~ContTreat, data = Dists)
summary(ORD_T)
confint(ORD_T)
AICc(ORD_T) #AICc = 240.0
AICc(ORD_null) - AICc(ORD_T)
anova(ORD_null, ORD_T, test = "F") #P = 0.001964

#Include sex as only predictor
ORD_S <- lm(PercORD~Sex, data = Dists)
summary(ORD_S) #Worsens model
confint(ORD_S)
AICc(ORD_S) #AICc = 250.7
AICc(ORD_null) - AICc(ORD_S)
anova(ORD_null, ORD_S, test = "F") #P = 0.7906

#Include sex and treatment
ORD_TS <- lm(PercORD~ ContTreat + Sex, data = Dists)
summary(ORD_TS) #Sex is not important
confint(ORD_TS)
AICc(ORD_TS) # AICc = 242.6
AICc(ORD_T) - AICc(ORD_S) #-10.71946
anova(ORD_T, ORD_TS, test = "F") # P = 0.553

#Include sex and treatment as an interaction
ORD_TxS <- lm(PercORD ~ ContTreat*Sex, data = Dists)
summary(ORD_TxS) #Only treatment is a important effect
AICc(ORD_TxS) #AICc = 245.7
anova(ORD_TxS, ORD_T, test = "F") #P > 0.05

#We conclude that treatment significantly increases ORD
#But sex doesn't have an impact


#Repeat for DD

#Develop NUll Model
DD_null <- lm(PercDD~1, data = Dists)
summary(DD_null)
AICc(DD_null) #AICc = 243.7

#Add the treatment predictor to the model
DD_T <- lm(PercDD~ContTreat, data = Dists)
summary(DD_T) #Treatment increases wariness
confint(DD_T) #49.97566 118.85928
AICc(DD_T) #AICc = 229.60
AICc(DD_null) - AICc(DD_T) #14.17816
anova(DD_null, DD_T, test = "F") #P = 0.000108

#Include sex but not treatment
DD_S <- lm(PercDD~Sex, data = Dists)
summary(DD_S) #Sex is not important; worsens model
AICc(DD_S) #AICc = 246.5
anova(DD_null, DD_S, test = "F") #P = 0.8887

#Include sex and treatment
DD_TS <- lm(PercDD~ ContTreat + Sex, data = Dists)
summary(DD_TS) #Sex is not important
confint(DD_TS)
AICc(DD_TS) # AICc = 231.3
AICc(DD_T) - AICc(DD_TS)
anova(DD_T, DD_TS, test = "F") # P = 0.2807

#Include sex and treatment as an interaction
DD_TxS <- lm(PercDD ~ ContTreat*Sex, data = Dists)
summary(DD_TxS) #Only treatment is an important
AICc(DD_TxS) #AICc = 234.7
anova(DD_TxS, DD_T, test = "F") #P = 0.575

#Get Rsquared for ORD_T and DD_T
summary(ORD_T) #Adjusted R-squared:  0.3572
summary(DD_T) #Adjusted R-squared:  0.5125

#Also run t-tests
t.test(Dists$PercORD~Dists$ContTreat) #t = -4.0313, df = 19.708, p-value = 0.0006702
t.test(Dists$PercDD~Dists$ContTreat) #t = -4.9197, df = 15.772, p-value = 0.0001604

###Now consider sound.
#Pull out the treatment bears (to compare sound treated and no sound treated bears)
SoundSS <- subset(Dists, Sound != "Control")

#Create Null (ORD)
SoundSSNull <- lm(PercORD ~ 1, data = SoundSS)
summary(SoundSSNull)

#Create model with sound as predictor (for ORD)
SoundSS1 <- lm(PercORD ~ Sound, data = SoundSS)
summary(SoundSS1) #Sound does not affect ORD

#Create Null (DD)
SoundSSNull2 <- lm(PercDD ~ 1, data = SoundSS)
summary(SoundSSNull2)

#Crate model with sound as predictor (For DD)
SoundSS1_2 <- lm(PercDD ~ Sound, data = SoundSS)
summary(SoundSS1_2) #Sound does not affect DD


#Owing to concerns about pseudoreplication, repeat analysis while randomly 
#selecting either 2007 or 2008 for bears F05, F10, M32, M39

#Generate random year for each:
floor(runif(4, min=2007, max=2009)) #2008 2008 2008 2007

#Remove rows 1, 3, 10, 22
Dists_Reduced <- Dists %>%
  dplyr::slice( -c(1,3,10,22))

#Develop a Null model for ORD
ORD_null2 <- lm(PercORD~1, data = Dists_Reduced)
summary(ORD_null2)
AICc(ORD_null2) #AICc = 202.667

#Add the treatment predictor to the model
ORD_T2 <- lm(PercORD~ContTreat, data = Dists_Reduced)
summary(ORD_T2) #Treatment matters
AICc(ORD_T2) #AICc = 199.7303
anova(ORD_null2, ORD_T2, test = "F") #P = 0.02469 

#Include sex but not treatment
ORD_S2 <- lm(PercORD~Sex, data = Dists_Reduced)
summary(ORD_S2) #Sex doesn't matter
AICc(ORD_S2) #AICc = 205.0148
anova(ORD_null2, ORD_S2, test = "F") #P = 0.4848

#Include sex and treatment
ORD_TS2 <- lm(PercORD~ ContTreat + Sex, data = Dists_Reduced)
summary(ORD_TS2) #Sex doesn't matter
AICc(ORD_TS2) # AICc = 201.8641
anova(ORD_T2, ORD_TS2, test = "F") # P = 0.3196

#Include sex and treatment as an interaction
ORD_TxS2 <- lm(PercORD ~ ContTreat*Sex, data = Dists_Reduced)
summary(ORD_TxS2) #no important effects
AICc(ORD_TxS2) #AICc = 205.0245
anova(ORD_TxS2, ORD_T2, test = "F") #0.4609

#Repeat for DD
#Create Null
DD_null2 <- lm(PercDD~1, data = Dists_Reduced)
summary(DD_null2)
AICc(DD_null2) #AICc = 199.6053

#Add treatment to model
DD_T2 <- lm(PercDD~ContTreat, data = Dists_Reduced)
summary(DD_T2) #Treatment is important
AICc(DD_T2) #AICc = 188.1817
anova(DD_null2, DD_T2, test = "F") #P = 0.0004342

#Include sex but not treatment
DD_S2 <- lm(PercDD~Sex, data = Dists_Reduced)
summary(DD_S2) #Sex is not  important
AICc(DD_S2) #AICc = 202.1792
anova(DD_null2, DD_S2, test = "F") #P = 0.5881

#Include sex and treatment
DD_TS2 <- lm(PercDD~ ContTreat + Sex, data = Dists_Reduced)
summary(DD_TS2) #Sex is not important
AICc(DD_TS2) # AICc = 189.9401
anova(DD_T2, DD_TS2, test = "F") # P = 0.2554

#Include sex and treatment as an interaction
DD_TxS2 <- lm(PercDD ~ ContTreat*Sex, data = Dists_Reduced)
summary(DD_TxS2) #Only treatment is an important effect
AICc(DD_TxS2) #AICc = 193.8585
anova(DD_TxS2, DD_T2, test = "F") #P = 0.5349


#Create charts showing change in ORD and DD for control vs. treatment bears

#First, develop relevant dataframe
Dists %>%
  group_by(ContTreat) %>%
  summarise_at(vars(PercORD), list(name = mean))

Dists %>%
  group_by(ContTreat) %>%
  summarise_at(vars(PercORD), list(name = sd))

Dists %>%
  group_by(ContTreat) %>%
  summarise_at(vars(PercDD), list(name = mean))

Dists %>%
  group_by(ContTreat) %>%
  summarise_at(vars(PercDD), list(name = sd))

Treatment <- c("Control", "Treatment")
meanORD <- c(-32.7, 46.5)
sdORD <- c(35.5, 56.5)
meanDD <- c(-15.4, 69.0)
sdDD <- c(37.5, 40.8)

ORDChartData <- cbind(Treatment, meanORD, sdORD)
ORDChartData <- as.data.frame (ORDChartData)
ORDChartData$meanORD <- as.numeric(ORDChartData$meanORD)
ORDChartData$sdORD <- as.numeric(ORDChartData$sdORD)

DDChartData <- cbind(Treatment, meanDD, sdDD)
DDChartData <- as.data.frame (DDChartData)
DDChartData$meanDD <- as.numeric(DDChartData$meanDD)
DDChartData$sdDD <- as.numeric(DDChartData$sdDD)

#Create bar chart using dataframe (ORD)
ORDChart1 <- ggplot(data = ORDChartData, mapping = aes(x = Treatment, y = meanORD)) +
  geom_bar(stat = "identity", colour = "black", fill = "darkgrey") +
  geom_errorbar(aes(ymin=meanORD-sdORD, ymax=meanORD+sdORD), width=.2,
                position=position_dodge(.9)) +
  geom_hline(yintercept=0) +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black", face = "plain"),
        axis.text.y = element_text(colour = "black", face = "plain")) +
  labs(x="Treatment group", y="Average percent change between\npre and post treatment periods")

ORDChart2 <- ORDChart1 + theme_classic()

ORDChart3 <- ORDChart2 + theme(axis.text.x = element_text(colour = "black", face = "plain", size = 12),
                               axis.text.y = element_text(colour = "black", face = "plain", size = 12),
                               axis.title.y = element_text(colour = "black", face = "plain", size = 12),
                               axis.title.x = element_text(colour = "black", face = "plain", size = 12)) +
  ggtitle("(A) Overt Reaction Distance")

#Create bar chart using dataframe (DD)
DDChart1 <- ggplot(data = DDChartData, mapping = aes(x = Treatment, y = meanDD)) +
  geom_bar(stat = "identity", colour = "black", fill = "darkgrey") +
  geom_errorbar(aes(ymin=meanDD-sdDD, ymax=meanDD+sdDD), width=.2,
                position=position_dodge(.9)) +
  geom_hline(yintercept=0) +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black", face = "plain"),
        axis.text.y = element_text(colour = "black", face = "plain")) +
  labs(x="Treatment group", y="Average percent change between\npre and post treatment periods")

DDChart2 <- DDChart1 + theme_classic()

DDChart3 <- DDChart2 + theme(axis.text.x = element_text(colour = "black", face = "plain", size = 12),
                               axis.text.y = element_text(colour = "black", face = "plain", size = 12),
                               axis.title.y = element_text(colour = "black", face = "plain", size = 12),
                               axis.title.x = element_text(colour = "black", face = "plain", size = 12)) +
  ggtitle("(B) Displacement Distance")


###Further consider wariness by looking at reactions to researchers prior to beginning AC
AC <- read.csv("AC_All.csv") #405 observations

#What factors predicted that bears would leave before AC could be conducted?
AC_ResReact <- AC
AC_ResReact$ResReact <- AC_ResReact$RXN.to.researchers
AC_ResReact <- AC_ResReact %>%
  dplyr::filter(ResReact != "DNR",
                ResReact != "Return") #reaction was recorded for 319 encounters

sum(AC_ResReact$ResReact == "Approach") #7 times 
sum(AC_ResReact$ResReact == "Bluff") #4 times 
sum(AC_ResReact$ResReact == "Indifferent") #61 times 
sum(AC_ResReact$ResReact == "Leave_run") #105 times 
sum(AC_ResReact$ResReact == "Leave_walk") #52 times 
sum(AC_ResReact$ResReact == "Stress") #4 times 
sum(AC_ResReact$ResReact == "Tree") #3 times 
sum(AC_ResReact$ResReact == "Unaware") #5 times 
sum(AC_ResReact$ResReact == "Wary") #78 times 

#Group outcome variable into desired outcome (1) and negative (0)
#1's will be leave-run, leave-walk, and tree
#0s will be approach, bluff, indifferent, stress, wary
#Do not include bears that were unaware
AC_ResReact$ResReact <- ifelse(AC_ResReact$ResReact == "Approach", 0, AC_ResReact$ResReact)
AC_ResReact$ResReact <- ifelse(AC_ResReact$ResReact == "Bluff", 0, AC_ResReact$ResReact)
AC_ResReact$ResReact <- ifelse(AC_ResReact$ResReact == "Indifferent", 0, AC_ResReact$ResReact)
AC_ResReact$ResReact <- ifelse(AC_ResReact$ResReact == "Stress", 0, AC_ResReact$ResReact)
AC_ResReact$ResReact <- ifelse(AC_ResReact$ResReact == "Wary", 0, AC_ResReact$ResReact)
AC_ResReact$ResReact <- ifelse(AC_ResReact$ResReact == "Leave_run", 1, AC_ResReact$ResReact)
AC_ResReact$ResReact <- ifelse(AC_ResReact$ResReact == "Leave_walk", 1, AC_ResReact$ResReact)
AC_ResReact$ResReact <- ifelse(AC_ResReact$ResReact == "Tree", 1, AC_ResReact$ResReact)

AC_ResReact <- AC_ResReact %>%
  dplyr::filter(ResReact != "Unaware") #314 observations
AC_ResReact$ResReact <- factor(AC_ResReact$ResReact, levels = c(0,1))
str(AC_ResReact)

#How many of each outcome
sum(AC_ResReact$ResReact == "0") #154 negative outcomes
sum(AC_ResReact$ResReact == "1") #160 positive outcomes

#Confirm that bear ID is important
ResReactTbl <- table(AC_ResReact$ResReact, AC_ResReact$BearID)
ResReactChi <- chisq.test(ResReactTbl)
ResReactChi #X-squared = 34.083, df = 23, p-value = 0.06397
ResReactChi$expected #many expected values are very small, so fisher test is appropriate
ResReactFish <- fisher.test(ResReactTbl, simulate.p.value = TRUE)
ResReactFish #p-value = 0.02949


#Do univariate tests to determine what affected this response
ResReactNull <- glmer(ResReact ~ (1|BearID) + 1, data = AC_ResReact, family = binomial)
summary(ResReactNull) 
AICc(ResReactNull) # 435.4901

#Does no. prevAC (month) affect it?
ResReactPAC <- glmer(ResReact ~ (1|BearID) + PrevAC_Month, data = AC_ResReact, family = binomial)
summary(ResReactPAC)
confint(ResReactPAC) #0.02435679  0.06574572
AICc(ResReactPAC) #418.5588
AICc(ResReactNull) - AICc(ResReactPAC) #16.9
anova(ResReactPAC, ResReactNull, test = "Chisq") #P = 1.328e-05

#Get effect size
summ(ResReactPAC, digits = 4, exp = TRUE, confint = TRUE) #1.0452 
#Does Sex affect response (univariate)?
ResReactSex <- glmer(ResReact ~ (1|BearID) + Sex, data = AC_ResReact, family = binomial)
summary(ResReactSex) #Sex doesn't matter
confint(ResReactSex) #Overlaps 0
AICc(ResReactNull) - AICc(ResReactSex) #-1.549675
summ(ResReactSex, exp = TRUE, confint = TRUE)
anova(ResReactSex, ResReactNull, test = "Chisq") #P = 1.328e-05

#Does sex improve upon previous model?
ResReactS <- glmer(ResReact ~ (1|BearID) + PrevAC_Month + Sex, data = AC_ResReact, 
                   family = binomial)
summary(ResReactS) #Doesn't matter
AICc(ResReactS) #419.348
AICc(ResReactPAC) - AICc(ResReactS) #-0.7892143
anova(ResReactPAC, ResReactS, test = "Chisq") #P = 0.2611

#Does sound treatment affect it (univariate)?
#First remove cases where sound wasn't recorded
AC_ResReact_Sound <- AC_ResReact %>%
  dplyr::filter(Whistle != "Unk")

#Now make a univariate
ResReactSound <- glmer(ResReact ~ (1|BearID) + Whistle, data = AC_ResReact_Sound,
                       family = binomial) 
summary(ResReactSound) ##beta = -0624, p =  0.0934 .
confint(ResReactSound)
AICc(ResReactSound) #433.2721
summ(ResReactSound, exp = TRUE, confint = TRUE)

#Make a model using sound and prev AC
ResReactSoundPrev <- glmer(ResReact ~ (1|BearID) + Whistle + PrevAC_Month,
                           data = AC_ResReact_Sound, family = binomial) 
confint(ResReactSoundPrev)
confint(ResReactSoundPrev, level = 0.85)
summary(ResReactSoundPrev) ##beta = -0.60284, p =  0.1059 .
AICc(ResReactSoundPrev) #416.9235

#Make a null model for comparison
ResReactSoundNull_L <- glmer(ResReact ~ (1|BearID) + 1, data = AC_ResReact_Sound,
                             family = binomial) 

AICc(ResReactSoundNull_L) - AICc(ResReactSound) #0.853
anova(ResReactSound, ResReactSoundNull_L, test = "Chisq") #P = 0.08905 
summ(ResReactSound, confint = TRUE, exp = TRUE)

#Create model using same dataset but including only prev AC (to act as null)
ResReactSoundNull <- glmer(ResReact ~ (1|BearID) + PrevAC_Month, 
                           data = AC_ResReact_Sound, family = binomial) 
summary(ResReactSoundNull) ##beta = -0.60284, p =  0.1059 .
AICc(ResReactSoundNull) #417.6291

anova(ResReactSoundNull, ResReactSoundPrev, test = "Chisq") #0.09679 
AICc(ResReactSoundNull) - AICc(ResReactSoundPrev) #0.7055302


#Does distace to cover affect outcome (univariate)?
#First remove cases where distcover wasn't recorded
AC_ResReact_DC <- AC_ResReact %>%
  dplyr::filter(DistCover != "DNR")
AC_ResReact_DC$DistCover <- as.numeric(AC_ResReact_DC$DistCover)

#Now make a univariate
ResReactDistC <- glmer(ResReact ~ (1|BearID) + DistCover, data = AC_ResReact_DC,
                       family = binomial)
summary(ResReactDistC) #beta = -0.01125 , P = 0.363
confint(ResReactDistC) #-0.03588762, 0.01297107
AICc(ResReactDistC) #379.9212
summ(ResReactDistC, exp = TRUE, confint = TRUE)
ResReactDistC_L <- glmer(ResReact ~ (1|BearID) + 1, data = AC_ResReact_DC,
                       family = binomial)

anova(ResReactDistC_L, ResReactDistC, test = "Chisq") #0.3617 
AICc(ResReactDistC_L) - AICc(ResReactDistC) #-1.21253


#What about Dist cover is added to the previous model
#Start by making a model where you have prevAC and Bear ID as predictors but
#Use reduced dataset where Dist Cover was recorded
ResReactPAC_2 <- glmer(ResReact ~ (1|BearID) + PrevAC_Month, data = AC_ResReact_DC, 
                       family = binomial)
ResReactDC <- glmer(ResReact ~ (1|BearID) + PrevAC_Month + DistCover, 
                    data = AC_ResReact_DC, family = binomial)
summary(ResReactDC) #dist cover is not a significant term
AICc(ResReactDC) #361.0543
AICc(ResReactPAC_2) - AICc(ResReactDC) #-1.487179
anova(ResReactPAC_2, ResReactDC, test = "Chisq") #P = 0.4492

#What about sex and distance to cover?
ResReactDCS <- glmer(ResReact ~ (1|BearID) + PrevAC_Month + DistCover + Sex, 
                     data = AC_ResReact_DC, family = binomial)
summary(ResReactDCS) #prev ac count is the only significant effect
AICc(ResReactDCS) #363.0621
AICc(ResReactPAC_2) - AICc(ResReactDCS) #-3.494966
anova(ResReactPAC_2, ResReactDCS, test = "Chisq") #P = 0.7261

#Does SD affect outcome (univariate)?
#First remove cases where SD wasn't recorded
AC_ResReact_SD <- AC_ResReact %>%
  dplyr::filter(SD != "NA")
AC_ResReact_SD$SD <- as.numeric(AC_ResReact_SD$SD)

#Create univariate model
ResReactSD <- glmer(ResReact ~ (1|BearID) + SD, data = AC_ResReact_SD,
                       family = binomial)
summary(ResReactSD) #Not a significant effect beta = 0.002397, P = 0.705
confint(ResReactSD) #-0.01027612 0.01559472
summ(ResReactSD, exp = TRUE, confint = TRUE)
ResReactSD_L <- glmer(ResReact ~ (1|BearID) + 1, data = AC_ResReact_SD,
                         family = binomial)

anova(ResReactSD_L, ResReactSD, test = "Chisq") #0.7042
AICc(ResReactSD_L) - AICc(ResReactSD) #--2.066246

#Including sex, distance to cover, sound, start distance all worsen model, so
#proceed with ResReactPAC
r2_nakagawa(ResReactPAC)
r2_nakagawa(ResReactSoundPrev)
#Conditional R2: 0.102 = 10.2% accounted for by bear ID AND PrevAC Count
#Marginal R2: 0.088 = 1.4% attributable to JUST bear ID

#Get summary metrics
summ(ResReactPAC, confint = TRUE, digits = 4)

#Get effect size
summ(ResReactPAC, confint = TRUE, digits = 4, exp = TRUE)
#effect size is 1.04, so for every additional AC event, the likelihood
#of positive outcome increases by 4%

#Create marginal effects plot
PrevAC1 <- plot_model(ResReactPAC, type = "pred", terms = "PrevAC_Month[all]",
                      title = "(C) Previous AC events")
PrevAC2 <- PrevAC1 + theme_classic()
PrevAC3 <- PrevAC2 + 
  theme(axis.text.x = element_text(colour = "black", face = "plain", size = 12),
        axis.text.y = element_text(colour = "black", face = "plain", size = 12),
        axis.title.y = element_text(colour = "black", face = "plain", size = 12),
        axis.title.x = element_text(colour = "black", face = "plain", size = 12)) +
  labs(x="No. AC events (Previous 30 days)", y="Likelihood that bear flees\nbefore AC begins")


#Create the three charts related to wariness into a single panel with three rows
CombinedPlot1 <- ggarrange(ORDChart3, DDChart3, PrevAC3, nrow = 3)
ggsave("Wariness_plot.png", CombinedPlot1, width = 3.5, height = 9, dpi = 700)


##Analysis 2 = predictors of successful AC events (unit of rep = AC event)
#Response variable = BearReactionSage
sum(AC$BearReactionSage != "DNR") #115 events where reaction was recorded
AC_React <- AC %>%
  dplyr::filter(BearReactionSage != "DNR") %>%
  dplyr::rename("Reaction" = BearReactionSage)

AC_React$Reaction <- ifelse(AC_React$Reaction == "Leave_walk", "Leave_Walk", AC_React$Reaction)

sum(AC_React$Reaction == "Leave_Run") #68 events where bear left at run
sum(AC_React$Reaction == "Leave_Walk") #16 events where bear left at walk
sum(AC_React$Reaction == "NoDisplacement") #31 events where bear didn't displace

#code outcome variables
AC_React$Runs <- AC_React$Reaction
AC_React$Runs <- ifelse(AC_React$Reaction == "Leave_Run", 1, 0)
AC_React$Runs <- factor(AC_React$Runs, levels = c(0,1))

#Confirm that individual ID affects outcome
ReactTbl <- table(AC_React$BearID, AC_React$Runs)
ReactChi <- chisq.test(ReactTbl)
ReactChi #X-squared = 34.511, df = 20, p-value = 0.02287
ReactChi$expected #Many small values, which means fisher test > chi square
ReactFisher <- fisher.test(ReactTbl, simulate.p.value=TRUE)
ReactFisher #p-value = 0.007496

#Does sex influence likelihood of '1' outcome?
#Create Null
AC_React_Null <- glmer(Runs ~ (1|BearID) + 1, data = AC_React, family = binomial)
summary(AC_React_Null)
AICc(AC_React_Null) #154.6649

#Create model with sex as only fixed effect
AC_React_Sex <- glmer(Runs~ (1|BearID) + Sex, data = AC_React, family = binomial)
summary(AC_React_Sex) #p = 0.241 for sex = male
confint(AC_React_Sex)
AICc(AC_React_Sex) #155.4376 
AICc(AC_React_Null) - AICc(AC_React_Sex) #-0.7727569
anova(AC_React_Null, AC_React_Sex, test = "Chisq") # P = 0.2477
summ(AC_React_Sex, exp = TRUE, digits = 4, confint = TRUE)

#Create model with previous AC as only fixed effect
AC_React_PAC <- glmer(Runs~ (1|BearID) + PrevAC_Month, data = AC_React, family = binomial)
summary(AC_React_PAC) #significant effect
confint(AC_React_PAC) #significant effect
AICc(AC_React_PAC) #141.4915 
AICc(AC_React_Null) - AICc(AC_React_PAC) #12.49211
anova(AC_React_Null, AC_React_PAC, test = "Chisq") #P = 0.0001328
summ(AC_React_PAC, exp = TRUE, digits = 4, confint = TRUE)

#Does distance to cover influence likelihood of flight?
#Remove instances where dist to cover wasn't recorded
AC_React_Cover <- AC_React %>%
  dplyr::filter(DistCover != "DNR")
AC_React_Cover$DistCover <- as.numeric(AC_React_Cover$DistCover)

AC_React_CoverMod <- glmer(Runs~ (1|BearID) + DistCover, data = AC_React_Cover, 
                           family = binomial)
summary(AC_React_CoverMod) #significant negative effect
confint(AC_React_CoverMod)
AICc(AC_React_CoverMod) #130.4203 
#Get effect size
summ(AC_React_CoverMod, exp = TRUE, digits = 4, confint = TRUE)

#Create appropriate null for comparison
AC_React_Null2 <- glmer(Runs ~ (1|BearID) + 1, data = AC_React_Cover, 
                        family = binomial)

#Compare using LR test and delta aicc
AICc(AC_React_Null2) - AICc(AC_React_CoverMod) #7.291216
anova(AC_React_Null2, AC_React_CoverMod, test = "Chisq") #P = 0.002154 

#Does projectile type influence likelihood of flight (only 37 observations)
AC_React_SG <- AC_React %>%
  dplyr::filter(ProjectileType != "None")
AC_React_SG$ProjectileType <- factor(AC_React_SG$ProjectileType, 
                                     levels = c("Slingshot", "Gun"))

#Create model
React_SG <- glmer(Runs~ (1|BearID) + ProjectileType, data = AC_React_SG, 
                  family = binomial)
summary(React_SG) #no effect of projectile type
confint(React_SG) 
summ(React_SG, exp = TRUE, digits = 4, confint = TRUE) #large effect size, but
#wide CIs; interpret with caution

#Create appropriate null for comparison
AC_React_Null3 <- glmer(Runs ~ (1|BearID) + 1, data = AC_React_SG, 
                        family = binomial)

#Compare using LR test and delta aicc
AICc(AC_React_Null3) - AICc(React_SG) #-0.01
anova(AC_React_Null3, React_SG, test = "Chisq") #P = 0.124 

#Confirm that projectile type isn't important using chi square test
Typetbl <- table(AC_React_SG$ProjectileType, AC_React_SG$Runs)
TypeChi <- chisq.test(Typetbl)
TypeChi #X-squared = 2.2564, df = 1, p-value = 0.1331 
TypeChi$expected


#Create final model including both prev_AC and distance to cover
AC_Final <- glmer(Runs~ (1|BearID) + DistCover + PrevAC_Month, 
                  data = AC_React_Cover, family = binomial)
summary(AC_Final)
summ(AC_Final, exp = TRUE, confint = TRUE, digits = 4)
summ(AC_Final, confint = TRUE, digits = 4)

#Let's compare to the more powerful of the two univariate GLMMs, but we need to
#creae a new one using this reduced (n = 103) data set
FinalNull <- glmer(Runs~ (1|BearID) + PrevAC_Month, 
                   data = AC_React_Cover, family = binomial)
AICc(FinalNull) - AICc(AC_Final) #3.078528
anova(FinalNull, AC_Final, test = "Chisq") #0.02202 

#Gte R2 for final model
r2_nakagawa(AC_Final)
#Conditional R2: 0.382 = 38.2% accounted for by bear ID AND fixed effects
#Marginal R2: 0.245 = 24.5% attributable to JUST fixed effects, 13.7 to bear ID

#Make a chart showing the important effects
PlotAC1 <- plot_model(AC_Final, type = "pred", terms = "PrevAC_Month[all]",
                      title = "(A) Previous AC events")
PlotAC2 <- PlotAC1 + theme_classic()
PlotAC3 <- PlotAC2 + 
  theme(axis.text.x = element_text(colour = "black", face = "plain", size = 12),
        axis.text.y = element_text(colour = "black", face = "plain", size = 12),
        axis.title.y = element_text(colour = "black", face = "plain", size = 12),
        axis.title.x = element_text(colour = "black", face = "plain", size = 12)) +
  labs(x="No. AC events (Previous 30 days)", y="Likelihood that bear leaves at a run")

PlotACD1 <- plot_model(AC_Final, type = "pred", terms = "DistCover[all]",
                      title = "(B) Distance to cover")
PlotACD2 <- PlotACD1 + theme_classic()
PlotACD3 <- PlotACD2 + 
  theme(axis.text.x = element_text(colour = "black", face = "plain", size = 12),
        axis.text.y = element_text(colour = "black", face = "plain", size = 12),
        axis.title.y = element_text(colour = "black", face = "plain", size = 12),
        axis.title.x = element_text(colour = "black", face = "plain", size = 12)) +
  labs(x="Distance to cover (m)", y="Likelihood that bear leaves at a run")


#Let's combine the two charts you've created that exhibit changes in response to AC
CombinedACPlot1 <- ggarrange(PlotAC3, PlotACD3, nrow = 2)
ggsave("ACSuccess_plot.png", CombinedACPlot1, width = 3.5, height = 6, dpi = 700)


#Third Analysis: Spatial and Temporal differences
#read in proportion/ distance data
PropDist <- read.csv("PropsDists2023.csv")

#Convert columns to appropriate data types
PropDist$Time <- factor(PropDist$Time, levels = c("Pre", "Post"))
PropDist$Period <- ifelse(PropDist$Period == "day  ", "Day", PropDist$Period)
PropDist$Period <- ifelse(PropDist$Period == "night", "Night", PropDist$Period)
PropDist$Period <- factor(PropDist$Period, levels = c("Day", "Night"))
PropDist$Treatment <- ifelse(PropDist$Treatment == "AC     ", "Treatment", PropDist$Treatment)
PropDist$Treatment <- ifelse(PropDist$Treatment == "control", "Control", PropDist$Treatment)
PropDist$Treatment <- factor(PropDist$Treatment, levels = c("Control", "Treatment"))

#Make dataframe for each of night and day
PropDistN <- subset(PropDist, PropDist$Period == "Night")
PropDistD <- subset(PropDist, PropDist$Period == "Day")

#During the DAY, is there a difference between pre/ post control/ treatment?
#Create a GLM
DayDistNull <- lmer(DistDev ~ (1|Bear.ID) + 1, data = PropDistD)
summ(DayDistNull)

DayDist1 <- lmer(DistDev ~ (1|Bear.ID) + Treatment*Time, data = PropDistD)
summ(DayDist1, digits = 4) 
confint(DayDist1)

AICc(DayDist1) - AICc(DayDistNull) #-32.88893
anova(DayDistNull, DayDist1, test = "F")

CatPlotDay <- 
  cat_plot(DayDist1, pred = Treatment, modx = Time, geom = "point", point.shape = TRUE,
         colors = c("black", "#696969"), vary.lty = TRUE)

CatPlotDay2 <- CatPlotDay + 
  theme_classic() +
  theme(axis.text.x = element_text(colour = "black", face = "plain", size = 12),
        axis.text.y = element_text(colour = "black", face = "plain", size = 12),
        axis.title.y = element_text(colour = "black", face = "plain", size = 12),
        axis.title.x = element_text(colour = "black", face = "plain", size = 12),
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", face = "plain", size = 12),
        legend.position = "none") +
  ylim(-200, 800) +
  labs(x="Treatment Group", y="Distance to human-\ndominated habitat (m)") +
  ggtitle("(A) Day")

#I will use means to compare among groups (effect sizes)
DistDayMean <- PropDistD %>%
  group_by(Treatment, Time) %>%
  summarise_at(vars(DistDev), list(name = mean))

#Percent diff for treatment
(DistDayMean[4,3] - DistDayMean[3,3])/ DistDayMean[3,3] * 100 #-50.47884
(DistDayMean[2,3] - DistDayMean[1,3])/ DistDayMean[1,3] * 100 #12.02977


#During the NIGHT, is there a difference between pre/ post control/ treatment?
NightDistNull <- lmer(DistDev ~ (1|Bear.ID) + 1, data = PropDistN)
summ(NightDistNull)

NightDist1 <- lmer(DistDev ~ (1|Bear.ID) + Treatment*Time, data = PropDistN)
summ(NightDist1, digits = 4) 
confint(NightDist1)
AICc(NightDist1) - AICc(NightDistNull)
anova(NightDistNull, NightDist1, test = "Chisq")

CatPlotNight <- 
  cat_plot(NightDist1, pred = Treatment, modx = Time, geom = "point", point.shape = TRUE,
           colors = c("black", "#696969"), vary.lty = TRUE)

CatPlotNight2 <- CatPlotNight + 
  theme_classic() +
  theme(axis.text.x = element_text(colour = "black", face = "plain", size = 12),
        axis.text.y = element_text(colour = "black", face = "plain", size = 12),
        axis.title.y = element_text(colour = "black", face = "plain", size = 12),
        axis.title.x = element_text(colour = "black", face = "plain", size = 12),
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", face = "plain", size = 12),
        legend.position = c(0.8,0.9)) +
  ylim(-200, 800) +
  labs(x="Treatment Group", y="Distance to human-\ndominated habitat (m)") +
  ggtitle("(B) Night")

DistNightMean <- PropDistN %>%
  group_by(Treatment, Time) %>%
  summarise_at(vars(DistDev), list(name = mean))

#Percent diff for treatment
(DistNightMean[4,3] - DistNightMean[3,3])/ DistNightMean[3,3] * 100 #-38.9
(DistNightMean[2,3] - DistNightMean[1,3])/ DistNightMean[1,3] * 100 #54.1


#REPEAT FOR PROPORTIONS
#During the DAY, is there a difference between pre/ post control/ treatment?
DayPropNull <- lmer(PropDev ~ (1|Bear.ID) + 1, data = PropDistD)
summ(DayPropNull)

DayProp1 <- lmer(PropDev ~ (1|Bear.ID) + Treatment*Time, data = PropDistD)
summ(DayProp1, digits = 4)
confint(DayProp1)

AICc(DayProp1) - AICc(DayPropNull)
anova(DayPropNull, DayProp1, test = "Chisq")


CatPlotDayProp <- 
  cat_plot(DayProp1, pred = Treatment, modx = Time, geom = "point", point.shape = TRUE,
           colors = c("black", "#696969"), vary.lty = TRUE)

CatPlotDayProp2 <- CatPlotDayProp + 
  theme_classic() +
  theme(axis.text.x = element_text(colour = "black", face = "plain", size = 12),
        axis.text.y = element_text(colour = "black", face = "plain", size = 12),
        axis.title.y = element_text(colour = "black", face = "plain", size = 12),
        axis.title.x = element_text(colour = "black", face = "plain", size = 12),
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", face = "plain", size = 12),
        legend.position = "none") +
  ylim(-0.1, 0.8) +
  labs(x="Treatment Group", y="Proportion human-dominated\nhabitat") +
  ggtitle("(C) Day")

PropDayMean <- PropDistD %>%
  group_by(Treatment, Time) %>%
  summarise_at(vars(PropDev), list(name = mean))

#Percent diff for treatment
(PropDayMean[4,3] - PropDayMean[3,3])/ PropDayMean[3,3] * 100 #-9.3
(PropDayMean[2,3] - PropDayMean[1,3])/ PropDayMean[1,3] * 100 #-20.3



#During the NIGHT, is there a difference between pre/ post control/ treatment?
NightPropNull <- lmer(PropDev ~ (1|Bear.ID) + 1, data = PropDistN)
summ(NightPropNull)

NightProp1 <- lmer(PropDev ~ (1|Bear.ID) + Treatment*Time, data = PropDistN)
summ(NightProp1, digits = 4)
confint(NightProp1)

AICc(NightProp1) - AICc(NightPropNull)
anova(NightPropNull, NightProp1, test = "Chisq")

CatPlotNightProp <- 
  cat_plot(NightProp1, pred = Treatment, modx = Time, geom = "point", point.shape = TRUE,
           colors = c("black", "#696969"), vary.lty = TRUE)

CatPlotNightProp2 <- CatPlotNightProp + 
  theme_classic() +
  theme(axis.text.x = element_text(colour = "black", face = "plain", size = 12),
        axis.text.y = element_text(colour = "black", face = "plain", size = 12),
        axis.title.y = element_text(colour = "black", face = "plain", size = 12),
        axis.title.x = element_text(colour = "black", face = "plain", size = 12),
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", face = "plain", size = 12),
        legend.position = "none") +
  ylim(-0.1, 0.8) +
  labs(x="Treatment Group", y="Proportion human-dominated\nhabitat") +
  ggtitle("(D) Night")

PropNightMean <- PropDistN %>%
  group_by(Treatment, Time) %>%
  summarise_at(vars(PropDev), list(name = mean))

#Percent diff for treatment
(PropNightMean[4,3] - PropNightMean[3,3])/ PropNightMean[3,3] * 100 #52.7
(PropNightMean[2,3] - PropNightMean[1,3])/ PropNightMean[1,3] * 100 #-8.5

#Make composite charts
CombinedCatPlot <- ggarrange(CatPlotDay2, CatPlotNight2, CatPlotDayProp2, CatPlotNightProp2, nrow = 2, ncol = 2)
ggsave("CatPlots.png", CombinedCatPlot, width = 7, height = 6, dpi = 700)
