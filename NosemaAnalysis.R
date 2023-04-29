##

#Analysis for lithium chloride manuscript

#Contact: cprouty@ufl.edu    https://scholar.google.com/citations?user=PpeDx78AAAAJ&hl=en

##

#load packages
library(lme4)
library(lsmeans)
library(multcomp)
library(brglm2)
library(afex)
library(survival)
library(coxme)
library(DHARMa)
library(car)
library(mgcv)
library(fitdistrplus)
library(pscl)
library(AER)
##
setwd("C:/Users/codyp/Desktop/WorkBackup/Research/HBREL/NosemaMS")
#Load datasheets
Consumption <- read.csv("Consumption.csv")
Intensity <- read.csv("Intensity.csv")
Prevalence <- read.csv("Prevalence.csv")
Survival <- read.csv("Survival1.csv")
HiveAlive <- read.csv("HiveAlive.csv")
HiveAliveV <- read.csv("HiveAliveV.csv")
Study1 <- read.csv("Study1Master.csv")
Study1V <- read.csv("Study1Vert.csv")
Study2 <- read.csv("Study2.csv")
###

#Data organization
Prevalence$Event <- as.factor(Prevalence$Event)
Prevalence$Cage <- as.factor(Prevalence$Cage)
Intensity$Event <- as.factor(Intensity$Event)
HiveAlive$MitesBees1 <- (HiveAlive$Mites1/HiveAlive$Bees1)*100
HiveAlive$MitesBees2 <- (HiveAlive$Mites2/HiveAlive$Bees2)*100
HiveAlive$MitesBees3 <- (HiveAlive$Mites3/HiveAlive$Bees3)*100
HiveAlive$MitesBees4 <- (HiveAlive$Mites4/HiveAlive$Bees4)*100
HiveAlive$MitesBees5 <- (HiveAlive$Mites5/HiveAlive$Bees5)*100
HiveAliveV$MitesBees <- (HiveAliveV$Mites/HiveAliveV$Bees) *100
HiveAliveV$Week <- as.factor(HiveAliveV$Week)
HiveAliveV$Spores1 <- HiveAliveV$Spores * 50000
names(Study1V)[1]<- "Treatment"
Study1V$Time <- as.factor(Study1V$Time)
names(Study2)[1] <- "Hive"
Study2$Time <- as.factor(Study2$Time)
names(Consumption)[1] <- "Treatment"
Intensity$Raw <- Intensity$Spores / 50000
Intensity1 <- subset(Intensity, Event == 1)
Intensity2 <- subset(Intensity, Event == 2)
names(Intensity)[1] <- "Treatment"
names(Prevalence)[1] <- "Treatment"
Prevalence1 <- subset(Prevalence, Event == 1)
Prevalence2 <- subset(Prevalence, Event == 2)
names(Study1)[1] <- "Treatment"
###

#Experiment 1 - Laboratory exposure to Nosema and treatment
#Pollen consumption
Pollen <- lmer(Pollen ~ Treatment + (1|Cage), data = Consumption)
anova(Pollen)
###

#Water consumption
Water <- lmer(Water ~ Treatment + (1|Cage), data = Consumption)
anova(Water)
###

#Sucrose syrup consumption
Syrup <- lmer(Syrup ~ Treatment + (1|Cage), data = Consumption)
anova(Syrup)
###

#Effect of treatment on intensity of Nosema
Int <- lm(Spores ~ Treatment*Event, data = Intensity)
anova(Int, test="F")

lsm<-lsmeans (Int, list( ~ Event))
cld(lsm)
###

#Effect of treatment on prevalence of Nosema
Prev <- glmer(Infected ~ Treatment*Event + (1|Cage), family=binomial(link="logit"), data = Prevalence)
Anova(Prev, type=3)



sim_bmond <- simulateResiduals(fittedModel = Prev, n = 250)
hist(sim_bmond)
plot(sim_bmond)
###

#Survival
m1 <- coxph(Surv(Longevity, Censor) ~ Treatment, data = Survival)
anova(m1)
summary(m1)
###

#Survival graph and analysis
tgcPrev <- summarySE(Prevalence, measurevar="Infected", groupvars=c("Treatment", "Event"))
Survival$SurvObj <- with(Survival, Surv(Longevity, Censor == 0))

survdiff(Surv(Longevity, Censor) ~ Treatment, Survival)

Survival2 <- survfit(SurvObj ~ Survival$Treatment, data = Survival)
plot(Survival2, col=c("#66CCEE", "#4477AA", "#CCBB44", "#EE6677"), lty=1,lwd=3,ylim = c(0.35,1), xlim=c(0,30), xlab="Time (days)", ylab="Proportion Survived")
legend("bottomleft", c("Control", "Fumagilin", "HoneyBHealthy", "Nozevit"), col=c("#66CCEE", "#4477AA", "#CCBB44", "#EE6677"), lty=1, lwd=3, title="Treatment")
###

#Experiment 2 - Field treatment: Hive Alive 
#Mites on treatment and time
Bees <- lmer(Bees ~ Treatment+(1|Hive), data = HiveAliveV)
anova(Bees)
#Bees do not differ significantly between treatments, so we can use the poisson distribution for raw mite counts
HiveAliveV$MitesBees.t <- HiveAliveV$MitesBees + 1
HiveAliveV$MitesBees.L <- log10(HiveAliveV$MitesBees.t)
Mites <- glmer(MBBin ~ Treatment*Week +(1|Hive), family=binomial(link="logit"), data = HiveAliveV)
Anova(Mites)

lsm<-lsmeans (Mites, list( ~ Treatment))
cld(lsm)
###

#Nosema on treatment and time
HiveAliveV$Spores1.t <- HiveAliveV$Spores1 + 1
HiveAliveV$Spores1.L <- log10(HiveAliveV$Spores1.t)
Nos <- lmer(Spores1.L ~ Treatment*Week + (1|Hive), data = HiveAliveV)
Anova(Nos)

lsm<-lsmeans (Nos, list( ~ Treatment*Week))
cld(lsm)
###

#Experiment 3 - Field treatment Fall vs Spring 1

#Brood estimation
Brood <- lm(cm2.brood ~ Treatment, data = Study1)
anova(Brood, test="F")
###

#Bees estimation
Bees <- lm(number.of.bees ~ Treatment, data = Study1)
anova(Bees, test="F")
###

#Number of Spores
Study1V$Spores.t <- Study1V$Spores + 1
Study1V$Spores.L <- log10(Study1V$Spores.t)
Spores <- lmer(Spores.L ~ Ingred*Time+(1|Colony), data = Study1V)
anova(Spores)

lsm<-lsmeans (Spores, list( ~ Time + Ingred))
cld(lsm)

sim_bmond <- simulateResiduals(fittedModel = Spores, n = 250)
hist(sim_bmond)
plot(sim_bmond)

#MCs are a mess, break datasheet by time sampled
Time0 <- subset(Study1V, Time == 0)
Time1 <- subset(Study1V, Time == 1)
Time2 <- subset(Study1V, Time == 2)
Time3 <- subset(Study1V, Time == 3)

Sp <- lm(Spores.L ~ Ingred, data = Time1)
anova(Sp)
summary(Sp)

lsm<-lsmeans (Sp, list( ~ Ingred))
cld(lsm)

#There were no significant differences within sampling times.

###

#Prevalence of Spores

#There is complete separation in a few groups (Prevalence is 1), so a bias-reduced glm is used for the analysis of prevalence
Prev <- glm(Prev ~ Treatment*Time, 
          family = binomial, Study1V, method=brglmFit)
anova(Prev, test = "Chisq")
###

#Experiment 4 - Field treatment Fall vs Spring 2

#Brood estimation
Brood <- lmer(Brood ~ Treatment*Time + (1|Hive), data = Study2)
anova(Brood, test="F")
###

#Bees estimation
Bees <- lmer(Bees ~ Treatment*Time + (1|Hive), data = Study2)
anova(Bees, test="F")
###

#Number of spores
Study2$Spores <- Study2$SporesRaw * 50000
Study2$Spores.t <- Study2$Spores + 1
Study2$Spores.L <- log10(Study2$Spores.t)
Spores <- lmer(Spores.L ~ Treatment*Time + (1|Hive), data = Study2)
Anova(Spores)

sim_bmond <- simulateResiduals(fittedModel = Spores, n = 250)
hist(sim_bmond)
plot(sim_bmond)
#Same as the previous experiment, not normally distributed. We will used the Poisson distribution

Spores <- mixed(Spores.t ~ Treatment*Time + (1|Hive), family=gaussian(link = "log"), data = Study2, method="LRT")
Spores

lsm<-lsmeans (Spores, list( ~ Time+Ingred))
cld(lsm)

#Split up by time sampled
Time0 <- subset(Study2, Time == 0)
Time1 <- subset(Study2, Time == 1)
Time2 <- subset(Study2, Time == 2)
Time4 <- subset(Study2, Time == 4)

Spores <- lm(Spores.L ~ Treatment, data = Time4)
anova(Spores)
#Significant differences only occur at time 1

lsm<-lsmeans (Spores, list( ~ Treatment))
cld(lsm)

#Control and Nosevit spring had greater spore counts than Fumagilin fall in the first time of sampling.

#Prevalence of Nosema
mixed(Prev ~ Treatment*Time + (1|Hive), family=binomial(link = "logit"), data = Study2, method="LRT")
###



