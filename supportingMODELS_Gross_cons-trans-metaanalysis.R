# R code for conservation translocation meta-analysis supplemental models

# Clear the environment
rm(list = ls())

# Set working directory
#FIXME: Uncomment and set working directory for local computer
#setwd("<<insert working directory here>>")

# load necessary packages
library("MCMCglmm")
library("metafor")


# Import data
load(file = "Gross_cons-trans-metaanalysis.rdata")
# Load previous models
load(file = "mcmc.taxo1.rdata")
load(file = "mcmc.taxo2.rdata")
load(file = "mcmc.taxo3.rdata")


########################################################################
######    #############################################################
#####  #  ############################################################
####  ##  #############      Meta-analsysis validation      #########
###         ########################################################
#######   #########################################################
##################################################################
# Setup MCMC specifications and priors
nsamp <- 3000
BURN <- 3000
THIN <- 150
NITT <- BURN + nsamp*THIN

# Priors
a <- 1000  #<-- scale factor for F-distribution priors
# Taxonomy
prFxdGR.taxo2 <- list(R = list(V = 1e-6, fix = 1),
G = list(G1 = list(V = diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*a),
        G2 = list(V = diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*a),
        G3 = list(V = diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*a),
        G4 = list(V = diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*a),
        G5 = list(V = diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*a),
        G6 = list(V = diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*a),
        G7 = list(V = diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*a),
        G8 = list(V = diag(1), fix = 1)))


## Meta-analysis validation and tests for publication bias
# The following four code blocks adapted from two sources:
# 1. Moatt et al. 2016 Dryad: http://dx.doi.org/10.5061/dryad.3fc02
# 2. Farquharson et al. 2018 suppl. code: https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-018-03500-9/MediaObjects/41467_2018_3500_MOESM5_ESM.txt

# 1. Time-lag bias
mcmc.timelag <- MCMCglmm(yi ~ PrintYear,
  random = ~ class + order + fam + genus + species + studyID + groupID + idh(SE):units,
  data = dat,
  prior = prFxdGR.taxo2,
  thin = THIN, burnin = BURN, nitt = NITT)

# 2. Egger's regression
prediction<-predict(mcmc.taxo1, marginal=~order + fam + genus + species + studyID + groupID + idh(SE):units)
precision<-sqrt(1/dat$vi)
mr <- dat$yi - prediction
z.mr <- mr * precision
lm.egger<-lm(z.mr ~ precision)
summary(lm.egger)

# 3. Funnel plot
graphics.off()
par(mfrow=c(1,2))
plot(dat$yi, precision, xlab = 'Log odds ratio', ylab = 'Precision')
abline(v=posterior.mode(mcmc.taxo1$Sol)[1],lwd=1) #v is the meta-analytic mean
plot(mr, precision, xlab = 'Meta-analytic residuals', ylab = 'Precision')
abline(v=0,lwd=1,lty=2)

# 4. Trim-and-fill analysis
metafor.tf <- rma(yi=mr,sei=1/precision)
summary(metafor.tf)
TFL <- trimfill(metafor.tf, side = "left", estimator = "R0")
TFR <- trimfill(metafor.tf, side = 'right', estimator = 'R0')
summary(TFL)
summary(TFR)
graphics.off()
range(mr)
funnel(TFL, xlab = "Meta-analytic residuals", xlim = c(-6,6))
funnel(TFR, xlab = "Meta-analytic residuals", xlim = c(-6, 6))

## Outlier identification
# based on studentized deleted residuals (Viechtbauer & Cheung 2010; DOI: 10.1002/jrsm.11)
# Code adapted from metafor documentation
metafor.res2 <- rma.mv(yi, vi, mods = cbind(class, scope, sex, strategy),
  random = list(~ 1 | order, ~ 1 | fam, ~ 1 | genus, ~ 1 | species, ~ 1 | studyID, ~ 1 | groupID), data=dat)

#XXX WARNING: the next line might take several HOURS !!!! XXX
rs.resid <- rstudent(metafor.res2)
save("metafor.res2", "rs.resid", file = "outlierID.rdata")
load("outlierID.rdata")

# NAs returned for two estimates' residuals. Not sure why
which(is.na(rs.resid$z))

# 43 estimates identified as outliers, based on z-scores > |1.96|
outliers <- which(rs.resid$z > 1.96 | rs.resid$z < -1.96)






#########################
# Running same model structure as taxo-1, 2, & 3, for meta-analysis validations with outliers omitted


# Prior from previous models:
prFxdGR.taxo <- list(R = list(V = 1e-6, fix = 1),
G = list(G1 = list(V = diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*a),
        G2 = list(V = diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*a),
        G3 = list(V = diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*a),
        G4 = list(V = diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*a),
        G5 = list(V = diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*a),
        G6 = list(V = diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*a),
        G7 = list(V = diag(1), fix = 1)))


# Subset data to remove outliers
dat.out <- dat[-outliers,]
dat.out <- droplevels(dat.out)
# Model 4: Wild vs. translocated comparisons
mcmc.taxo4 <- MCMCglmm(yi ~ class + scope + sex + strategy +
    sex*scope + sex*strategy + scope*strategy,
  random = ~ order + fam + genus + species + studyID + groupID + idh(SE):units,
  data = dat.out,
  prior = prFxdGR.taxo,
  thin = THIN, burnin = BURN, nitt = NITT)


# Model 5: Model 4 subset of studies for which number of generations in captivity is known
dat.out.capt <- dat.out[which(dat.out$capt.present == 1),]
dat.out.capt[] <- droplevels(dat.out.capt)

mcmc.taxo5 <- MCMCglmm(yi ~ class + scope + sex + capt.bin + 
    sex*scope + sex*capt.bin + scope*capt.bin,
  random = ~ order + fam + genus + species + studyID + groupID + idh(SE):units,
  data = dat.out.capt,
  prior = prFxdGR.taxo,
  thin = THIN, burnin = BURN, nitt = NITT)

# Model 6: Model 4 subset of studies involving both enriched and unenriched comparisons
dat.out.enrich <- dat.out[which(dat.out$enrich.present == 1),]
dat.out.enrich[] <- droplevels(dat.out.enrich)

mcmc.taxo6 <- MCMCglmm(yi ~ class + enrich.bin + scope + sex + strategy +
    scope*enrich.bin + sex*enrich.bin + enrich.bin*strategy,
  random = ~ order + fam + genus + species + studyID + groupID + idh(SE):units,
  data = dat.out.enrich,
  prior = prFxdGR.taxo,
  thin = THIN, burnin = BURN, nitt = NITT)

## Model validation

# 1. Time-lag bias
mcmc.timelag.out <- MCMCglmm(yi ~ PrintYear,
  random = ~ class + order + fam + genus + species + studyID + groupID + idh(SE):units,
  data = dat.out,
  prior = prFxdGR.taxo2,
  thin = THIN, burnin = BURN, nitt = NITT)

# 2. Egger's regression
prediction.out<-predict(mcmc.taxo7, marginal=~order + fam + genus + species + studyID + groupID + idh(SE):units)
precision.out<-sqrt(1/dat.out$vi)
mr.out <- dat.out$yi - prediction.out
z.mr.out <- mr.out * precision.out
lm.egger.out<-lm(z.mr.out ~ precision.out)
summary(lm.egger.out)

# 3. Funnel plot
graphics.off()
par(mfrow=c(1,2))
plot(dat.out$yi, precision.out, xlab = 'Log odds ratio', ylab = 'Precision')
abline(v=posterior.mode(mcmc.taxo7$Sol)[1],lwd=1) #v is the meta-analytic mean
plot(mr.out, precision.out, xlab = 'Meta-analytic residuals', ylab = 'Precision')
abline(v=0,lwd=1,lty=2)

# 4. Trim-and-fill analysis
metafor.tf.out <- rma(yi=mr.out,sei=1/precision.out)
summary(metafor.tf.out)
TFL.out <- trimfill(metafor.tf.out, side = "left", estimator = "R0")
TFR.out <- trimfill(metafor.tf.out, side = 'right', estimator = 'R0')
summary(TFL.out)
summary(TFR.out)
graphics.off()
range(mr.out)
funnel(TFL.out, xlab = "Meta-analytic residuals", xlim = c(-6,6))
funnel(TFR.out, xlab = "Meta-analytic residuals", xlim = c(-6, 6))


# Save several models in one object:
## for plotting and summarizing in another file
save(list = c("mcmc.taxo1", "mcmc.taxo2", "mcmc.taxo3",
    "mcmc.taxo4", "mcmc.taxo5", "mcmc.taxo6",
    "rs.resid", 
    "dat.out", "dat.out.capt", "dat.out.enrich"),
  file = "allModels_Gross_cons-trans-metaanalysis.rdata")

 
