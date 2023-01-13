# R code to generate results and plots for conservation translocation meta-analysis
# Clear the environment
rm(list = ls())

# Set working directory
#FIXME: Uncomment and set working directory for local computer
#setwd("<<insert working directory here>>")
setwd("/Users/iwo/ACT_MetaAnalysis")

# load necessary packages
library("MCMCglmm")
library("metafor")
library("emmeans")
library("ggplot2")
library("ggridges")
library("reshape")
library("wesanderson")

source("plot_utils.R")  #<-- necessary to plot 95% HPD intervals in Figures 1 & 2


# Import data
load(file = "Gross_cons-trans-metaanalysis.rdata")
load(file = "allModels_Gross_cons-trans-metaanalysis.rdata")

#######################
# Table 1: Comparison sample sizes organized by taxonomic class and the performance metric measured.
addmargins(table(dat$class, dat$scope))

#######################
# Findings reporting

# Demonstrating findings reporting for an excerpt from the second paragraph of the Results
# "The model’s overall log odds ratio... 
est.lor <- emmeans(object = mcmc.taxo1, specs = ~1, data = dat)
hpd.summary(est.lor, point.est = posterior.mode)

# ...was converted to an odds ratio by exponentiating across the posterior distribution...
est.or <- exp(as.mcmc.emmGrid(est.lor))
posterior.mode(est.or)
HPDinterval(est.or)

# ...Indicating that translocated organisms have 69.9% decreased odds of out-performing their wild-resident counterparts.
1-posterior.mode(est.or)

# By quantifying the proportion of posterior samples with LOR estimates greater than or equal to zero,
# we infer that translocated organisms had an overall 0.002 probability of out-performing wild-resident conspecifics following release."
mean(as.mcmc.emmGrid(est.lor) >= 0)
mean(est.or >= 1)

############################################
# Figure 1: Overall effects across classes 

# Pull marginals posterior modes for overall model effect
taxo.o.1 <- emmeans(object = mcmc.taxo1, specs = ~1, data = dat)
df.o.1 <- data.frame(hpd.summary(taxo.o.1, point.est = posterior.mode))

# Return the number of unique studies and estimates across class levels
n.stu.1 <- length(unique(dat$studyID))
n.est.1 <- length(unique(dat$estimateID))
df.o.1 <- cbind(df.o.1, n.est.1, n.stu.1)

# Pull marginals posterior modes for class fixed effect
taxo.c.1 <- emmeans(object = mcmc.taxo1, specs = ~class, data = dat)
df.c.1 <- data.frame(hpd.summary(taxo.c.1, point.est = posterior.mode))

# Return the number of unique studies and estimates across class levels
n.stu.c.1 <- as.vector(by(dat$studyID,list(dat$class), function(x) length(unique(x))))
n.est.c.1 <- as.vector(table(dat$class))
df.c.1 <- cbind(df.c.1, n.est.c.1, n.stu.c.1)

# rbind class and overall dfs
df.1 <- data.frame(mapply(c, df.c.1, df.o.1, SIMPLIFY = FALSE))
levels(df.1$class)[levels(df.1$class) == "overall"] <- "Overall"

# Tack on parentheticals stating n studies and n estimates
df.1$ns <- NA
df.1$ns <- paste0("(",df.1$n.stu[], ", ", df.1$n.est[],")")
df.1$class <- factor(df.1$class, levels = rev(levels(df.1$class)))

# Bind columns of posterior estimates, melt into proper format for plotting
taxo.c.em.1 <- as.data.frame(as.mcmc.emmGrid(taxo.c.1))
taxo.o.em.1 <- as.data.frame(as.mcmc.emmGrid(taxo.o.1))
colnames(taxo.c.em.1)<-gsub("class ","",colnames(taxo.c.em.1))
colnames(taxo.o.em.1)<-gsub("1 overall","Overall",colnames(taxo.o.em.1))
taxo.c.em.1 <- cbind(taxo.c.em.1, taxo.o.em.1)
melt.c.em.1 <- melt(taxo.c.em.1)


## Figure 1 plot
ridges.1 <- ggplot(melt.c.em.1, aes(x=value, y = variable)) +
     coord_cartesian(ylim=c(1.4, 10.2)) +
     xlim(-9, 7) +
     scale_y_discrete(limits = unique(rev(melt.c.em.1$variable))) +
         geom_density_ridges2(rel_min_height = 5e-8, scale =1.45,
         quantile_lines = TRUE,
         calc_ecdf = TRUE,
         quantiles = c(0.025, 0.975),
         fill = wes_palette("Darjeeling1")[4],
         alpha = 0.7) +
     theme_classic() +
     geom_vline(xintercept = 0, linetype = 3, size = 1.2) +
     ylab("") + xlab("ln odds ratio (95% CrI)") +
     geom_text(data = df.1,aes(y = class, x = upper.HPD),
         label = df.1$ns,
         nudge_y = .33,
         nudge_x = 2)

# Posterior mode lines
temp.1 <- ggplot_build(ridges.1)
build.1 <- as.data.frame(temp.1$data[1])

res.1 <- by(build.1, build.1$group, function(i) i[which(i$density == max(i$density)),])
density_lines.1 <- do.call(rbind, res.1)

# Add 95% credible intervals to plot
## calculate for each group
hpd_data <- aggregate(value ~ variable, data = melt.c.em.1,
  FUN = function(x) HPDinterval(as.mcmc(x), prob = 0.95))
  ## add "group" column to correspond to `class` convention
  hpd_data$group <- as.factor(as.integer(hpd_data$variable))
## Add Credible interval limits to plot data
build.1[, c("lcl", "ucl")] <- hpd_data[match(build.1$group, hpd_data$group),
    "value"]

# replot, but this time add 95% credible interval bars along each plot/ridges
## x-axis BEFORE plotting the density    
ggplot(melt.c.em.1, aes(x=value, y = variable)) +
     coord_cartesian(ylim=c(1.4, 10.2)) +
     xlim(-9, 7) +
     scale_y_discrete(limits = unique(rev(melt.c.em.1$variable))) +
    # 95% credible intervals as thick lines along bottom of each density
    ## plot before the densities so 95% CrI bars not over top density below
    geom_segment(data = build.1, 
      aes(x = lcl, xend = ucl, y = ymin, yend = ymin), size = 3) +         
    geom_density_ridges2(rel_min_height = 5e-8, scale =1.45,
         stat = "density_ridges_HPDCrI", 
         quantile_lines = TRUE,
         calc_ecdf = TRUE,
#         quantiles = c(0.025, 0.975),
	 quantiles = c(0.95),
	 quantile_fun = HPD_fun,
         fill = wes_palette("Darjeeling1")[4],
         alpha = 0.7) +
    # posterior mode lines
    geom_segment(data = density_lines.1, 
               aes(x = x, y = ymin, xend = x, 
                   yend = ymin+density*scale*iscale)) +
                   
     theme_classic() +
     geom_vline(xintercept = 0, linetype = 3, size = 1.2) +
     ylab("") + xlab("ln odds ratio (95% CrI)") +
     geom_text(data = df.1,aes(y = class, x = upper.HPD),
         label = df.1$ns,
         nudge_y = .33,
         nudge_x = 2)



################################
# Figure 2: Influence of enrichment across classes

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

# Add 95% credible intervals to plot
## calculate for each group and enrichment level
hpd_data2 <- aggregate(value ~ class + enrich, data = melt.em,
  FUN = function(x) HPDinterval(as.mcmc(x), prob = 0.95))
  ## add "group" column to correspond to `class` convention
  ### first, recreate "variable" column
  hpd_data2$variable <- with(hpd_data2,
    paste0("class ", class, ", enrich.bin ", enrich))
    ## go back by hand and change "Overall" to drop all but enrich part
    ### clunky, but to match how it comes out in melt.em
    hpd_data2$variable <- gsub(pattern = "class Overall, ", replacement = "",
      x = hpd_data2$variable)
  hpd_data2$group <- match(hpd_data2$variable, levels(melt.em$variable))
## Add Credible interval limits to plot data
build.2[, c("lcl", "ucl")] <- hpd_data2[match(build.2$group, hpd_data2$group),
    "value"]
## also need to add enrichment column for "fill" argument
build.2$enrich <- hpd_data2[match(build.2$group, hpd_data2$group), "enrich"] 

# replot, but this time add 95% credible interval bars along each plot/ridges
## x-axis BEFORE plotting the density    
ggplot(melt.em, aes(x=value, y = class, fill = as.factor(enrich))) +
      coord_cartesian(ylim=c(1.5, 6.3)) +
      xlim(-13,13) +
      scale_y_discrete(limits = unique(rev(melt.em$class))) +
      scale_fill_manual(values = wes_palette("Darjeeling1"),
          labels=c('Enrichment absent', 'Enrichment present')) +
      geom_density_ridges2(rel_min_height = 5e-8,
          scale = 1.15,
          stat = "density_ridges_HPDCrI", 
          quantile_lines = TRUE,
          calc_ecdf = TRUE,
# 	  quantiles = c(0.025, 0.975),         
          quantiles = 0.95,
          quantile_fun = HPD_fun,
          alpha = 0.7) +
   # posterior mode lines
   geom_segment(data = density_lines.2, 
            aes(x = x, y = ymin, xend = x, 
                yend = ymin+density*scale*iscale)) +

      theme_classic() +
      theme(legend.position = "bottom",
         legend.title = element_blank()) +
      geom_vline(xintercept = -0.1, linetype = 3, size = 1.2) +
      ylab("") + xlab("ln odds ratio (95% CrI)") +
      geom_text(data = df[1:6,], aes(y = class, x = upper.HPD),
         inherit.aes=FALSE, parse=FALSE,
         label = df[1:6,]$ns, 
         nudge_y =0.25, nudge_x = 5)



#####################################
# Calculation of model heterogeneity 

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
