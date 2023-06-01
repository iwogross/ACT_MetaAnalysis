
########################################################################
######   ##############################################################
####     #############################################################
######   #########           Data preparation           #############
######   ###########################################################
####       ########################################################
##################################################################
# R code for conservation translocation meta-analysis

# Clear the environment
rm(list = ls())

# Set working directory
#FIXME: Uncomment and set working directory for local computer
#setwd("<<insert working directory here>>")


# load necessary packages
library("MCMCglmm")
library("metafor")


# Import data
dat <- read.table(file = "DATA_Gross_cons-trans-metaanalysis.txt",
  header = TRUE, sep = "\t")


# Explicitly format all variables 
dat <- within(dat, {
  estimateID <- as.factor(estimateID)
  estimator <- as.factor(estimator)
  studyID <- as.factor(studyID)
  groupID <- as.factor(groupID)
  strategy <- as.factor(strategy)
  scope <- factor(scope,
    levels = c("Movement", "Constitution", "Growth", "Survival", "Reproduction"))
  enrich <- as.factor(enrich)
  enrich.present <- as.factor(enrich.present)
  capt.present <- as.factor(capt.present)
  enrich.bin <- as.factor(enrich.bin)
  capt.gens <- as.factor(capt.gens)
  capt.bin <- as.factor(capt.bin)
  class <- factor(class,
    levels = c("Bivalvia", "Gastropoda", "Malacostraca", "Echinoidea",
      "Actinopterygii", "Amphibia", "Reptilia", "Aves", "Mammalia"))
  phylum <- as.factor(phylum)
  order <- as.factor(order)
  fam <- as.factor(fam)
  genus <- as.factor(genus)
  species <- as.factor(species)
  sex <- factor(sex, levels = c("Female", "Mixed", "Male"))
  n1i <- as.numeric(n1i)
  m1i <- as.numeric(m1i)
  sd1i <- as.numeric(sd1i)
  n2i <- as.numeric(n2i)
  m2i <- as.numeric(m2i)
  sd2i <- as.numeric(sd2i)
  A <- as.integer(A)
  B <- as.integer(B)
  C <- as.integer(C)
  D <- as.integer(D)
  PrintYear <- as.factor(PrintYear)
})




########################################################################
######     ############################################################
####   ##   ##########################################################
########   ##########         Calculate effect sizes         ########
######   ###########        with R package 'metafor'        ########
####        #######################################################
##################################################################

# Parse Log odds ratio estimates to allow separate analyses on one dataset
LORrows <- which(dat$estimator == 'LOR')

# Calculate standardized mean difference estimates
dat <- escalc(measure="D2ORN", m1i = m2i, sd1i = sd2i, n1i = n2i, m2i = m1i, sd2i = sd1i, n2i = n1i, data=dat)

# Calculate log odds ratios and convert to estimates of standardized mean difference
dat[LORrows, ] <- escalc(measure="OR", ai = C, bi = D, ci = A, di = B, data=dat[LORrows, ])

# Calculate standard errors on the calculated effect sizes
dat$SE <- sqrt(dat$vi)


# Create data subsets
## subset of studies for which number of generations in captivity is known (model 2)
dat.capt <- dat[which(dat$capt.present == 1),]
dat.capt <- droplevels(dat.capt)

## subset of studies involving both enriched and unenriched comparisons (model 3)
dat.enrich <- dat[which(dat$enrich.present == 1),]
dat.enrich <- droplevels(dat.enrich)


# Now save as R data object to keep above variable assignments and new variables
## will use this object in other files of this repository\
save("dat", "dat.capt", "dat.enrich",
  file = "Gross_cons-trans-metaanalysis.rdata")



########################################################################
#####       ###########################################################
########   ###########################################################
######    #############       Conduct meta-analysis         #########
########   ###########      with R package 'MCMCglmm'      #########
####      #########################################################
##################################################################


# Setup MCMC specifications and priors
nsamp <- 3000
BURN <- 3000
THIN <- 500
NITT <- BURN + nsamp*THIN

# Priors
a <- 1000  #<-- scale factor for F-distribution priors
# Taxonomy
prFxdGR.taxo <- list(R = list(V = 1e-6, fix = 1),
G = list(G1 = list(V = diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*a),
        G2 = list(V = diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*a),
        G3 = list(V = diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*a),
        G4 = list(V = diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*a),
        G5 = list(V = diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*a),
        G6 = list(V = diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*a),
        G7 = list(V = diag(1), fix = 1)))

# Model 1: Wild vs. translocated comparisons
mcmc.taxo1 <- MCMCglmm(yi ~ class + scope + sex + strategy +
    sex*scope + sex*strategy + scope*strategy,
  random = ~ order + fam + genus + species + studyID + groupID + idh(SE):units,
  data = dat,
  prior = prFxdGR.taxo,
  thin = THIN, burnin = BURN, nitt = NITT)

save("mcmc.taxo1", file = "mcmc.taxo1.rdata")

# Model 2: M1 subset of studies for which number of generations in captivity is known
mcmc.taxo2 <- MCMCglmm(yi ~ class + scope + sex + capt.bin + 
    sex*scope + sex*capt.bin + scope*capt.bin,
  random = ~ order + fam + genus + species + studyID + groupID + idh(SE):units,
  data = dat.capt,
  prior = prFxdGR.taxo,
  thin = THIN, burnin = BURN, nitt = NITT)

save("mcmc.taxo2", file = "mcmc.taxo2.rdata")


# Model 3: M1 subset of studies involving both enriched and unenriched comparisons

mcmc.taxo3 <- MCMCglmm(yi ~ class + enrich.bin + scope + sex + strategy +
    scope*enrich.bin + sex*enrich.bin + enrich.bin*strategy,
  random = ~ order + fam + genus + species + studyID + groupID + idh(SE):units,
  data = dat.enrich,
  prior = prFxdGR.taxo,
  thin = THIN, burnin = BURN, nitt = NITT)

save("mcmc.taxo3", file = "mcmc.taxo3.rdata")



