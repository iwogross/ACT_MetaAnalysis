# R code for conservation translocation meta-analysis supporting/validation models

# Clear the environment
rm(list = ls())

# Set working directory
#FIXME: Uncomment and set working directory for local computer
#setwd("<<insert working directory here>>")


# load necessary packages
library("MCMCglmm")
library("metafor")
library("emmeans")
library("ggplot2")
library("ggridges")
library("reshape")


source("plot_utils.R")  #<-- necessary to plot 95% HPD intervals in Figures 4 & 5

# Snag some colors from base R's `Okabe-Ito` palette
## suitable palette for color-related vision accessibility
### first have a look
palOI <- palette.colors(NULL, "Okabe-Ito")
pie(rep(1, length(palOI)), col = palOI, labels = names(palOI))
# grab subset of 3 colors for figures
palOI3 <- palOI[c("bluishgreen", "vermillion", "blue")]
  class(palOI3) <- "palette"
  attr(palOI3, "name") <- "Okabe-Ito3"
  
dev.off()  #<-- turn off pie chart


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

#XXX
#if already run below models and saved, then just load instead of re-running
# Go to next section labelled "XXX Post-Model Run"
#load(file = "mcmc.timelag.rdata"); load(file = "mcmc.timelag.out.rdata")  
#load(file = "allModels_Gross_cons-trans-metaanalysis.rdata")


# Setup MCMC specifications and priors
nsamp <- 3000
BURN <- 3000
THIN <- 500
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





#XXX Post-Model Run
# Pull marginals posterior modes for overall model effect
taxo.o <- emmeans(object = mcmc.timelag, specs = ~ 1, data = dat)
df.o <- data.frame(hpd.summary(taxo.o, point.est = posterior.mode))

# Return the number of unique studies and estimates for the overall effect
n.stu.o <- length(unique(dat$studyID))
n.est.o <- length(unique(dat$estimateID))
df.o <- cbind(df.o, n.est.o, n.stu.o)

# Pull marginals posterior modes for class fixed effect
taxo.y <- emmeans(object = mcmc.timelag, specs = ~PrintYear, data = dat)
df.y <- data.frame(hpd.summary(taxo.y, point.est = posterior.mode))

# Return the number of unique studies and estimates across class levels
n.stu.y <- as.vector(by(dat$studyID,list(dat$PrintYear), function(x) length(unique(x))))
n.est.y <- as.vector(table(dat$PrintYear))
df.y <- cbind(df.y, n.est.y, n.stu.y)

# rbind class and overall dfs
df <- data.frame(mapply(c, df.y, df.o, SIMPLIFY = FALSE))
levels(df$PrintYear)[levels(df$PrintYear) == "overall"] <- "Overall"

df$ns <- NA
df$ns <- paste0("(",df$n.stu[], ", ", df$n.est[],")")
df$PrintYear <- factor(df$PrintYear, levels = rev(levels(df$PrintYear)))

# Fancy posterior distribution by by-class
taxo.y.em <- as.data.frame(as.mcmc.emmGrid(taxo.y))
taxo.o.em <- as.data.frame(as.mcmc.emmGrid(taxo.o))
colnames(taxo.y.em) <- gsub("PrintYear ", "", colnames(taxo.y.em))
colnames(taxo.o.em) <- gsub("1 overall", "Overall", colnames(taxo.o.em))
taxo.y.em <- cbind(taxo.y.em, taxo.o.em)
melt.y.em <- melt(taxo.y.em)

  

# Fig. 4
pdf(file = "./Fig4.pdf",
  width = 6, height = 8)
ggplot(melt.y.em, aes(x = value, y = variable)) +
      xlim(-4, 4) +
      scale_y_discrete(limits = unique(rev(melt.y.em$variable))) +
    geom_vline(xintercept = 0, linetype = "51", #<--lty is units on then units off
      linewidth = 1.75, color = "grey10") +
    geom_density_ridges2(rel_min_height = 5.5e-3, scale = 1.5,
        stat = "density_ridges_HPDCrI",
        quantile_lines = TRUE,
        quantiles = c(0.95),
        quantile_fun = HPD_fun,
        fill = palOI3[3],
        alpha = 0.75) +  #<-- transparency (0 more -- 1 less/solid)
        
    theme_classic() +
    ylab("Year") + xlab("ln odds ratio (95% CrI)")

dev.off()



# 2. Egger's regression
prediction <- predict(mcmc.taxo1, marginal = ~order + fam + genus + species + studyID + groupID + idh(SE):units)
precision <- sqrt(1/dat$vi)
mr <- dat$yi - prediction
z.mr <- mr * precision
lm.egger <- lm(z.mr ~ precision)
summary(lm.egger)

# 3. Funnel plot

## Fig. 5
pdf(file = "./Fig5.pdf",
  width = 9, height = 6)

par(mfrow = c(1,2), mar = c(5, 5, 0.5, 1), cex.lab = 1.5, cex.axis = 1.0)
  plot(precision ~ dat$yi, type = "n", axes = FALSE,
    xlab = "ln odds ratio", ylab = "Precision")
  abline(v = posterior.mode(mcmc.taxo1$Sol)[1], lwd = 3) #v is meta-analytic mean
  points(precision ~ dat$yi, pch = 21, lwd = 1.6, col = "grey50")
  axis(1)
  axis(2)
  mtext("A", side = 3, line = -1.9, outer = TRUE, at = 0.02, cex = 2)
        
  plot(precision ~ mr, type = "n", axes = FALSE,
    xlab = "Meta-analytic residuals", ylab = "Precision")
  abline(v = 0, lwd = 3)
  points(precision ~ mr, pch = 21, lwd = 1.6, col = "grey50")
  axis(1)
  axis(2)
  mtext("B", side = 3, line = -1.9, outer = TRUE, at = 0.52, cex = 2)
  
dev.off() 



  
# 4. Trim-and-fill analysis
metafor.tf <- rma(yi = mr, sei = 1/precision)
summary(metafor.tf)
TFL <- trimfill(metafor.tf, side = "left", estimator = "R0")
TFR <- trimfill(metafor.tf, side = "right", estimator = "R0")
summary(TFL)
summary(TFR)

range(mr)

par(mfrow = c(1, 2))
  funnel(TFL, xlab = "Meta-analytic residuals", xlim = c(-6, 6))
  funnel(TFR, xlab = "Meta-analytic residuals", xlim = c(-6, 6))








################################################################################



############################
## Outlier identification ##
############################
# based on studentized deleted residuals (Viechtbauer & Cheung 2010; DOI: 10.1002/jrsm.11)
# Code adapted from metafor documentation

#XXX
#if already run below models and saved, then just load instead of re-running
# Go to next section labelled "XXX Post-Model Run"

metafor.res2 <- rma.mv(yi, vi, mods = cbind(class, scope, sex, strategy),
  random = list(~ 1 | order, ~ 1 | fam, ~ 1 | genus, ~ 1 | species, ~ 1 | studyID, ~ 1 | groupID), data=dat)

#XXX WARNING: the next line might take 12 HOURS !!!! XXX
## if on Linux/Unix consider: 
### rs.resid <- rstudent(metafor.res2, parallel = "multicore", ncpus = 3, progbar = TRUE)
rs.resid <- rstudent(metafor.res2, progbar = TRUE)
save("metafor.res2", "rs.resid", file = "outlierID.rdata")
#load("outlierID.rdata")

# Outliers identified based on z-scores > |1.96|
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

save("mcmc.taxo4", file = "mcmc.taxo4.rdata")


# Model 5: Model 4 subset of studies for which number of generations in captivity is known
dat.out.capt <- dat.out[which(dat.out$capt.present == 1),]
dat.out.capt[] <- droplevels(dat.out.capt)

mcmc.taxo5 <- MCMCglmm(yi ~ class + scope + sex + capt.bin + 
    sex*scope + sex*capt.bin + scope*capt.bin,
  random = ~ order + fam + genus + species + studyID + groupID + idh(SE):units,
  data = dat.out.capt,
  prior = prFxdGR.taxo,
  thin = THIN, burnin = BURN, nitt = NITT)

save("mcmc.taxo5", file = "mcmc.taxo5.rdata")

# Model 6: Model 4 subset of studies involving both enriched and unenriched comparisons
dat.out.enrich <- dat.out[which(dat.out$enrich.present == 1),]
dat.out.enrich[] <- droplevels(dat.out.enrich)

mcmc.taxo6 <- MCMCglmm(yi ~ class + enrich.bin + scope + sex + strategy +
    scope*enrich.bin + sex*enrich.bin + enrich.bin*strategy,
  random = ~ order + fam + genus + species + studyID + groupID + idh(SE):units,
  data = dat.out.enrich,
  prior = prFxdGR.taxo,
  thin = THIN, burnin = BURN, nitt = NITT)

save("mcmc.taxo6", file = "mcmc.taxo6.rdata")



# 1. Time-lag bias
mcmc.timelag.out <- MCMCglmm(yi ~ PrintYear,
  random = ~ class + order + fam + genus + species + studyID + groupID + idh(SE):units,
  data = dat.out,
  prior = prFxdGR.taxo2,
  thin = THIN, burnin = BURN, nitt = NITT)






#XXX Post-Model Run
# Pull marginals posterior modes for overall model effect
taxo.o.out <- emmeans(object = mcmc.timelag.out, specs = ~1, data = dat.out)
df.o.out <- data.frame(hpd.summary(taxo.o.out, point.est = posterior.mode))

# Return the number of unique studies and estimates for the overall effect
n.stu.o.out <- length(unique(dat.out$studyID))
n.est.o.out <- length(unique(dat.out$estimateID))
df.o.out <- cbind(df.o.out, n.est.o.out, n.stu.o.out)

# Pull marginals posterior modes for class fixed effect
taxo.y.out <- emmeans(object = mcmc.timelag.out, specs = ~PrintYear, data = dat.out)
df.y.out <- data.frame(hpd.summary(taxo.y.out, point.est = posterior.mode))

# Return the number of unique studies and estimates across class levels
n.stu.y.out <- as.vector(by(dat.out$studyID,list(dat.out$PrintYear), function(x) length(unique(x))))
n.est.y.out <- as.vector(table(dat.out$PrintYear))
df.y.out <- cbind(df.y.out, n.est.y.out, n.stu.y.out)

# rbind class and overall dfs
df.out <- data.frame(mapply(c, df.y.out, df.o.out, SIMPLIFY = FALSE))
levels(df.out$PrintYear)[levels(df.out$PrintYear) == "overall"] <- "Overall"

df.out$ns <- NA
df.out$ns <- paste0("(",df.out$n.stu[], ", ", df.out$n.est[],")")
df.out$PrintYear <- factor(df.out$PrintYear, levels = rev(levels(df.out$PrintYear)))

# Fancy posterior distribution by by-class
taxo.y.out.em <- as.data.frame(as.mcmc.emmGrid(taxo.y.out))
taxo.o.out.em <- as.data.frame(as.mcmc.emmGrid(taxo.o.out))
colnames(taxo.y.out.em)<-gsub("PrintYear ","",colnames(taxo.y.out.em))
colnames(taxo.o.out.em)<-gsub("1 overall","Overall",colnames(taxo.o.out.em))
taxo.y.out.em <- cbind(taxo.y.out.em, taxo.o.out.em)
melt.y.out.em <- melt(taxo.y.out.em)



# Supplementary Figure 1 (See creation of Fig. 4 above)
pdf(file = "./Supplementary_Fig1.pdf",
  width = 6, height = 8)
ggplot(melt.y.out.em, aes(x = value, y = variable)) +
    xlim(-4, 4) +
    scale_y_discrete(limits = unique(rev(melt.y.out.em$variable))) +
    geom_vline(xintercept = 0, linetype = "51", #<--lty is units on then units off
      size = 1.75, color = "grey10") +
    geom_density_ridges2(rel_min_height = 5.5e-3, scale = 1.5,
    	stat = "density_ridges_HPDCrI",
        quantile_lines = TRUE,
        quantiles = c(0.95),
        quantile_fun = HPD_fun,
        fill = palOI3[3],  #<-- transparency (0 more -- 1 less/solid)
        alpha = 0.75) +
        
    theme_classic() +
    ylab("Year") + xlab("ln odds ratio (95% CrI)")

dev.off()   
    

    
    
    
    
# 2. Egger's regression
prediction.out <- predict(mcmc.taxo4, marginal = ~order + fam + genus + species + studyID + groupID + idh(SE):units)
precision.out<-sqrt(1/dat.out$vi)
mr.out <- dat.out$yi - prediction.out
z.mr.out <- mr.out * precision.out
lm.egger.out <- lm(z.mr.out ~ precision.out)
summary(lm.egger.out)

# 3. Funnel plot
## Supplemental Fig. 2
pdf(file = "./Supplementary_Fig2.pdf",
  width = 9, height = 6)

par(mfrow = c(1,2), mar = c(5, 5, 0.5, 1), cex.lab = 1.5, cex.axis = 1.0)
  plot(precision.out ~ dat.out$yi, type = "n", axes = FALSE,
    xlim = c(-4.5, 4.1),
    xlab = "ln odds ratio", ylab = "Precision")
  abline(v = posterior.mode(mcmc.taxo4$Sol)[1], lwd = 3) #v is meta-analytic mean
  points(precision.out ~ dat.out$yi, pch = 21, lwd = 1.6, col = "grey50")
  axis(1)
  axis(2)
  mtext("A", side = 3, line = -1.9, outer = TRUE, at = 0.02, cex = 2)
      
  plot(precision.out ~ mr.out, type = "n", axes = FALSE,
    xlim = c(-4.5, 4.1),
    xlab = "Meta-analytic residuals", ylab = "Precision")
  abline(v = 0, lwd = 3)
  points(precision.out ~ mr.out, pch = 21, lwd = 1.6, col = "grey50")
  axis(1)
  axis(2)
  mtext("B", side = 3, line = -1.9, outer = TRUE, at = 0.52, cex = 2)
  
dev.off() 



# 4. Trim-and-fill analysis
metafor.tf.out <- rma(yi = mr.out, sei = 1/precision.out)
summary(metafor.tf.out)
TFL.out <- trimfill(metafor.tf.out, side = "left", estimator = "R0")
TFR.out <- trimfill(metafor.tf.out, side = "right", estimator = "R0")
summary(TFL.out)
summary(TFR.out)

range(mr.out)

par(mfrow = c(1, 2))
  funnel(TFL.out, xlab = "Meta-analytic residuals", xlim = c(-6,6))
  funnel(TFR.out, xlab = "Meta-analytic residuals", xlim = c(-6, 6))



################################################################################
# Save several models in one object:
## for plotting and summarizing later (or from code in another file)
save(list = c("mcmc.taxo1", "mcmc.taxo2", "mcmc.taxo3",
    "mcmc.taxo4", "mcmc.taxo5", "mcmc.taxo6",
    "mcmc.timelag", "mcmc.timelag.out",
    "rs.resid", "outliers",
    "dat.out", "dat.out.capt", "dat.out.enrich"),
  file = "allModels_Gross_cons-trans-metaanalysis.rdata")

 
