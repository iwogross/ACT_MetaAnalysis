# R code to generate results and plots for conservation translocation meta-analysis
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

# Import data
load(file = "Gross_cons-trans-metaanalysis.rdata")

# Import data
load("meta_mcmc_dat.rdata")
load("mcmc.taxo1.rdata")
load("mcmc.taxo2.rdata")
load("mcmc.taxo3.rdata")
load("mcmc.taxo7.rdata")
load("mcmc.taxo8.rdata")
load("mcmc.taxo9.rdata")
load("outlierID.rdata")
load("Egger_prediction.rdata")
load("Egger_prediction.out.rdata")
load("mcmc.timelag.rdata")
load("mcmc.timelag.out.rdata")
load("trim-and-fill.rdata")
load("trim-and-fill_out.rdata")



########################################
## Calculation of model heterogeneity ##
########################################
# Source: Kamiya et al. 2014 (https://datadryad.org/stash/dataset/doi:10.5061/dryad.288j7)

# Weight and mean variance
Weight<-1/dat$vi
MV<-sum(Weight*(length(Weight)-1))/(sum(Weight)^2-sum(Weight^2))

## I2 order
I2.order<-100*(mcmc.taxo1$VCV[,"order"])/(mcmc.taxo1$VCV[,"order"]+mcmc.taxo1$VCV[,"fam"]+mcmc.taxo1$VCV[,"genus"]+mcmc.taxo1$VCV[,"species"]+mcmc.taxo1$VCV[,"studyID"]+mcmc.taxo1$VCV[,"groupID"]+mcmc.taxo1$VCV[,"units"]+MV)
summary(I2.order)
posterior.mode(I2.order)
HPDinterval(I2.order)

# I2 family
I2.fam<-100*(mcmc.taxo1$VCV[,"fam"])/(mcmc.taxo1$VCV[,"order"]+mcmc.taxo1$VCV[,"fam"]+mcmc.taxo1$VCV[,"genus"]+mcmc.taxo1$VCV[,"species"]+mcmc.taxo1$VCV[,"studyID"]+mcmc.taxo1$VCV[,"groupID"]+mcmc.taxo1$VCV[,"units"]+MV)
summary(I2.fam)
posterior.mode(I2.fam)
HPDinterval(I2.fam)

# I2 genus
I2.genus<-100*(mcmc.taxo1$VCV[,"genus"])/(mcmc.taxo1$VCV[,"order"]+mcmc.taxo1$VCV[,"fam"]+mcmc.taxo1$VCV[,"genus"]+mcmc.taxo1$VCV[,"species"]+mcmc.taxo1$VCV[,"studyID"]+mcmc.taxo1$VCV[,"groupID"]+mcmc.taxo1$VCV[,"units"]+MV)
summary(I2.genus)
posterior.mode(I2.genus)
HPDinterval(I2.genus)

# I2 species
I2.species<-100*(mcmc.taxo1$VCV[,"species"])/(mcmc.taxo1$VCV[,"order"]+mcmc.taxo1$VCV[,"fam"]+mcmc.taxo1$VCV[,"genus"]+mcmc.taxo1$VCV[,"species"]+mcmc.taxo1$VCV[,"studyID"]+mcmc.taxo1$VCV[,"groupID"]+mcmc.taxo1$VCV[,"units"]+MV)
summary(I2.species)
posterior.mode(I2.species)
HPDinterval(I2.species)

# I2 studyID
I2.study<-100*(mcmc.taxo1$VCV[,"studyID"])/(mcmc.taxo1$VCV[,"order"]+mcmc.taxo1$VCV[,"fam"]+mcmc.taxo1$VCV[,"genus"]+mcmc.taxo1$VCV[,"species"]+mcmc.taxo1$VCV[,"studyID"]+mcmc.taxo1$VCV[,"groupID"]+mcmc.taxo1$VCV[,"units"]+MV)
summary(I2.study)
posterior.mode(I2.study)
HPDinterval(I2.study)

# I2 groupID
I2.group<-100*(mcmc.taxo1$VCV[,"groupID"])/(mcmc.taxo1$VCV[,"order"]+mcmc.taxo1$VCV[,"fam"]+mcmc.taxo1$VCV[,"genus"]+mcmc.taxo1$VCV[,"species"]+mcmc.taxo1$VCV[,"studyID"]+mcmc.taxo1$VCV[,"groupID"]+mcmc.taxo1$VCV[,"units"]+MV)
summary(I2.group)
posterior.mode(I2.group)
HPDinterval(I2.group)

# I2 residual
I2.res<-100*(mcmc.taxo1$VCV[,"units"])/(mcmc.taxo1$VCV[,"order"]+mcmc.taxo1$VCV[,"fam"]+mcmc.taxo1$VCV[,"genus"]+mcmc.taxo1$VCV[,"species"]+mcmc.taxo1$VCV[,"studyID"]+mcmc.taxo1$VCV[,"groupID"]+mcmc.taxo1$VCV[,"units"]+MV)
summary(I2.res)
posterior.mode(I2.res)
HPDinterval(I2.res)

# I2 total
I2.total<-100*(mcmc.taxo1$VCV[,"units"]+mcmc.taxo1$VCV[,"order"]+mcmc.taxo1$VCV[,"fam"]+mcmc.taxo1$VCV[,"genus"]+mcmc.taxo1$VCV[,"species"]+mcmc.taxo1$VCV[,"studyID"]+mcmc.taxo1$VCV[,"groupID"])/(mcmc.taxo1$VCV[,"order"]+mcmc.taxo1$VCV[,"fam"]+mcmc.taxo1$VCV[,"genus"]+mcmc.taxo1$VCV[,"species"]+mcmc.taxo1$VCV[,"studyID"]+mcmc.taxo1$VCV[,"groupID"]+mcmc.taxo1$VCV[,"units"]+MV)
summary(I2.total)
posterior.mode(I2.total)
HPDinterval(I2.total)

# I2 taxonomic signal
I2.taxo<-100*(mcmc.taxo1$VCV[,"order"]+mcmc.taxo1$VCV[,"fam"]+mcmc.taxo1$VCV[,"genus"]+mcmc.taxo1$VCV[,"species"])/(mcmc.taxo1$VCV[,"order"]+mcmc.taxo1$VCV[,"fam"]+mcmc.taxo1$VCV[,"genus"]+mcmc.taxo1$VCV[,"species"]+mcmc.taxo1$VCV[,"studyID"]+mcmc.taxo1$VCV[,"groupID"]+mcmc.taxo1$VCV[,"units"]+MV)
summary(I2.taxo)
posterior.mode(I2.taxo)
HPDinterval(I2.taxo)



addmargins(table(dat$class, dat$scope))
by(dat$studyID,list(dat$class), function(x) length(unique(x)))

# Return the number of unique studies and estimates across class levels
n.stu <- as.vector(by(dat$studyID, list(dat$class), function(x) length(unique(x))))
n.est <- as.vector(table(dat$class))
df.c <- cbind(df.c,n.stu, n.est)
df.c[,3:5] <- round(df.c[,3:5],2)
df.c


#FOR LOG ODDS PLOTS
bonk <- emmeans(object = mcmc.taxo1, specs = ~class, data = dat)
df.c <- data.frame(hpd.summary(bonk, point.est = posterior.mode))
n.est <- length(unique(dat$estimateID))
df.c <- cbind(df.c, n.est)
df.c[,3:5] <- round(df.c[,3:5],2)

bonk <- emmeans(object = mcmc.taxo1, specs = ~1, data = dat)

df.c <- data.frame(hpd.summary(as.mcmc.emmGrid(bonk), point.est = posterior.mode))

# FOR ODDS RATIO TABLES
bonk <- emmeans(object = mcmc.taxo2, specs = ~capt.bin*scope, data = dat.capt)
df.c <- data.frame(hpd.summary(bonk, point.est = posterior.mode))
n.est <- length(unique(dat$estimateID))
df.c <- cbind(df.c, n.est)
df.c[,3:5] <- round(df.c[,3:5],2)

# Overall effect table
taxo.c <- emmeans(object = mcmc.taxo1, specs = ~1, data = dat)
df.c <- data.frame(hpd.summary(taxo.c, point.est = posterior.mode))
n.stu <- length(unique(dat$studyID))
n.est <- length(unique(dat$estimateID))
df.c <- cbind(df.c,n.stu, n.est)
df.c[,2:4] <- round(df.c[,2:4],2)
df.c

#One-way effect table
taxo.c <- emmeans(object = mcmc.taxo1, specs = ~class, data = dat)
df.o <- data.frame(hpd.summary(taxo.c, point.est = posterior.mode))
n.stu <- as.vector(by(dat$studyID, list(dat$class), function(x) length(unique(x))))
n.est <- as.vector(table(dat$class))
df.o <- cbind(df.o, n.stu, n.est)
df.o[,2:4] <- round(df.o[,2:4],2)
df.o

# Interaction effect tables
taxo.c <- emmeans(object = mcmc.taxo1, specs = ~1, data = dat)
df.c <- data.frame(hpd.summary(taxo.c, point.est = posterior.mode))
n.stu <- as.vector(by(dat$studyID, list(dat$scope, dat$sex), function(x) length(unique(x))))
n.est <- as.vector(table(dat$scope, dat$sex))
df.c <- cbind(df.c,n.stu, n.est)
df.c[,3:5] <- round(df.c[,3:5],2)
df.c

#Log odds ratio
bonk <- emmeans(object = mcmc.taxo2, specs = ~scope*capt.bin, data = dat.capt)
hpd.summary(bonk, point.est = posterior.mode)
bonk.or <- exp(as.mcmc.emmGrid(bonk))

#mean(as.mcmc.emmGrid(bonk)[,1] >= 1)

# Odds ratio
posterior.mode(bonk.or)
HPDinterval(bonk.or)
# "decreased odds"
1-posterior.mode(bonk.or)
mean(bonk.or[5] >= 1)




################################################################################
# DUMP from "...-plots.R"
#### Figure 1: Overall effects across classes ####

# Pull marginals posterior modes for overall model effect
taxo.o <- emmeans(object = mcmc.taxo1, specs = ~1, data = dat)
df.o <- data.frame(hpd.summary(taxo.o, point.est = posterior.mode))

# Return the number of unique studies and estimates across class levels
n.stu <- length(unique(dat$studyID))
n.est <- length(unique(dat$estimateID))
df.o <- cbind(df.o, n.est, n.stu)

# Pull marginals posterior modes for class fixed effect
taxo.c <- emmeans(object = mcmc.taxo1, specs = ~class, data = dat)
df.c <- data.frame(hpd.summary(taxo.c, point.est = posterior.mode))

# Return the number of unique studies and estimates across class levels
n.stu <- as.vector(by(dat$studyID,list(dat$class), function(x) length(unique(x))))
n.est <- as.vector(table(dat$class))
df.c <- cbind(df.c, n.est, n.stu)

# rbind class and overall dfs
df <- data.frame(mapply(c, df.c, df.o, SIMPLIFY = FALSE))
levels(df$class)[levels(df$class) == "overall"] <- "Overall"

# Tack on parentheticals stating n studies and n estimates
df$ns <- NA
df$ns <- paste0("(",df$n.stu[], ", ", df$n.est[],")")
df$class <- factor(df$class, levels = rev(levels(df$class)))

# Bind columns of posterior estimates, melt into proper format for plotting
taxo.c.em <- as.data.frame(as.mcmc.emmGrid(taxo.c))
taxo.o.em <- as.data.frame(as.mcmc.emmGrid(taxo.o))
colnames(taxo.c.em)<-gsub("class ","",colnames(taxo.c.em))
colnames(taxo.o.em)<-gsub("1 overall","Overall",colnames(taxo.o.em))
taxo.c.em <- cbind(taxo.c.em, taxo.o.em)
melt.c.em <- melt(taxo.c.em)


## Figure 1 plot
ridges.1 <- ggplot(melt.c.em, aes(x=value, y = variable)) +
     coord_cartesian(ylim=c(1.4, 10.2)) +
     xlim(-9, 7) +
     scale_y_discrete(limits = unique(rev(melt.c.em$variable))) +
         geom_density_ridges2(rel_min_height = 5e-8, scale =1.45,
         quantile_lines = TRUE,
         calc_ecdf = TRUE,
         quantiles = c(0.025, 0.975),
         fill = wes_palette("Darjeeling1")[4],
         alpha = 0.7) +
     theme_classic() +
     geom_vline(xintercept = 0, linetype = 3, size = 1.2) +
     #geom_point(data = aggregate(value ~ variable, melt.c.em, posterior.mode), # My best attempt at visualizing the posterior mode
      #  aes(x = value, y = variable),
      #  color = "black", shape = 17,
      #  size = 3,
      #  inherit.aes = FALSE) +
     ylab("") + xlab("ln odds ratio (95% CrI)") +
     geom_text(data = df,aes(y = class, x = upper.HPD),
         label = df$ns,
         nudge_y = .33,
         nudge_x = 2)

# Posterior mode lines
temp.1 <- ggplot_build(ridges.1)
build.1 <- as.data.frame(temp.1$data[1])

res <- by(build.1, build.1$group, function(i) i[which(i$density == max(i$density)),])
density_lines.1 <- do.call(rbind, res)

ridges.1 + geom_segment(data = density_lines.1, 
               aes(x = x, y = ymin, xend = x, 
                   yend = ymin+density*scale*iscale))

### Figure 2: Class x enrichment ####

# Pull marginals posterior modes for overall model effect
taxo.o <- emmeans(object = mcmc.taxo3, specs = ~enrich.bin, data = dat.enrich)
df.o <- data.frame(class = as.factor(rep("Overall", 2)), hpd.summary(taxo.o, point.est = posterior.mode))

# Return the number of unique studies and estimates across class levels
n.stu <- length(unique(dat.enrich$studyID))
n.est <- length(unique(dat.enrich$estimateID))
df.o <- cbind(df.o, n.est, n.stu)

# Pull marginals posterior modes for class fixed effect
taxo.c <- emmeans(object = mcmc.taxo3, specs = ~class*enrich.bin, data = dat.enrich)
df.c <- data.frame(hpd.summary(taxo.c, point.est = posterior.mode))

# Return the number of unique studies and estimates across class levels
n.stu <- as.vector(by(dat.enrich$studyID,list(dat.enrich$class), function(x) length(unique(x))))
n.est <- as.vector(table(dat.enrich$class))
df.c <- cbind(df.c, n.est, n.stu)

# rbind class and overall dfs
df <- rbind(df.c[1:5,], df.o[1,],df.c[6:10,],df.o[2,])

# Tack on parentheticals stating n studies and n estimates overall and by class
df$ns <- NA
df$ns <- paste0("(",df$n.stu[], ", ", df$n.est[],")")
df$class <- factor(df$class, levels = rev(levels(df$class)))

# Bind columns of posterior estimates, melt into proper format for plotting
taxo.c.em <- as.data.frame(as.mcmc.emmGrid(taxo.c))
taxo.o.em <- as.data.frame(as.mcmc.emmGrid(taxo.o))
taxo.em <- cbind(taxo.c.em[,1:5], taxo.o.em[,1, drop = FALSE], taxo.c.em[,6:10], taxo.o.em[,2, drop = FALSE])
melt.em <- melt(taxo.em)
melt.em <- cbind(melt.em, enrich = rep(0:1, each = 18000), class = rep(as.character(df$class), each = 3000))


# Figure 2 plot
ridges.2 <- ggplot(melt.em, aes(x=value, y = class, fill = as.factor(enrich))) +
      coord_cartesian(ylim=c(1.5, 6.3)) +
      xlim(-13,13) +
      scale_y_discrete(limits = unique(rev(melt.em$class))) +
      scale_fill_manual(values = wes_palette("Darjeeling1"),
          labels=c('Enrichment absent', 'Enrichment present')) +
      geom_density_ridges2(rel_min_height = 5e-8,
          scale = 1.15,
          quantile_lines = TRUE,
          calc_ecdf = TRUE,
          quantiles = c(0.025, 0.975),
          alpha = 0.7) +
      theme_classic() +
      theme(legend.position = "bottom",
         legend.title = element_blank()) +
      geom_vline(xintercept = 0, linetype = 3, size = 1.2) +
      ylab("") + xlab("ln odds ratio (95% CrI)") +
      geom_text(data = df[1:6,], aes(y = class, x = upper.HPD),
         inherit.aes=FALSE, parse=FALSE,
         label = df[1:6,]$ns, 
         nudge_y =0.25, nudge_x = 5)

# Posterior mode lines
temp.2 <- ggplot_build(ridges.2)
build.2 <- as.data.frame(temp.2$data[1])

res <- by(build.2, build.2$group, function(i) i[which(i$density == max(i$density)),])
density_lines.2 <- do.call(rbind, res)

density_lines.2 <- density_lines.2[order(density_lines.2$ymin),]
density_lines.2$group <- rep(1:6, each = 2)
density_lines.2$group <- as.factor(density_lines.2$group)
density_lines.2$enrich <- rep(0:1, 6)
density_lines.2$enrich <- as.factor(rep(0:1, 6))

ridges.2 + geom_segment(data = density_lines.2, 
               aes(x = x, y = ymin, xend = x, 
                   yend = ymin+density*scale*iscale))




############################
#### Findings reporting ####
############################


## Line 182-183: translocated organisms have 67% decreased odds of out-performing their wild-resident counterparts.
## Below I've included a progression of calculations based on Todd Steury's notes regarding logistic regression output
## Did I exponentiate in the right place??

est.1 <- emmeans(object = mcmc.taxo1, specs = ~1, data = dat)

## "Wild-residents have 167% increased odd of out-performing translocates"
posterior.mode(exp(as.mcmc.emmGrid(est.1)))
#1 overall 
#2.670096

## "Translocates have 37% the odds of out-performing wild-residents"
1/posterior.mode(exp(as.mcmc.emmGrid(est.1)))
#1 overall 
#0.3745183

## The final product
## "Translocates have 63% decreased odds of out-performing wild-residents"
1-(1/posterior.mode(exp(as.mcmc.emmGrid(est.1))))
#1 overall 
#0.6254817




## Lines 210-215
## Did I do these correctly?

## "the posterior mode of our overall model’s marginal log odds ratio estimate (LOR: 0.98; CrI: 0.40 to 1.68)"
posterior.mode(as.mcmc.emmGrid(est.1))
#1 overall 
#0.9803469

HPDinterval(as.mcmc.emmGrid(est.1))
#              lower    upper
#1 overall 0.4046175 1.684843
#attr(,"Probability")
#[1] 0.95

## "Thus, we find that wild-residents have 167% increased odds of out-performing their translocated conspecifics (OR: 2.67; CrI: 1.32 to 4.95)"
posterior.mode(exp(as.mcmc.emmGrid(est.1)))
#1 overall 
#2.670096

HPDinterval(exp(as.mcmc.emmGrid(est.1)))
#             lower    upper
#1 overall 1.325463 4.950923
#attr(,"Probability")
#[1] 0.95

## "... we can infer that wild-residents have a 0.9997 probability of out-performing wild-resident conspecifics following release"
mean(as.mcmc.emmGrid(est.1) > 0)
#[1] 0.9996667




