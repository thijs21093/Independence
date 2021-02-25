library(coda)
library(lattice)
library(MCMCpack)
library(cowplot)
library(bayesplot)
library(ggplot2)
library(mcmcse)
library(polycor)
library(dplyr)
library(extrafont)
library(MCMCvis)
library(gridExtra)
library(grid)
library(readxl)
library(tidyverse)


## Graphical settings
color_scheme_set("darkgray")
par(mar=c(1,1,1,1)) # Margins
options(scipen=999)

## Before starting
setwd("C:/Users/Thijs/surfdrive/Project A - Accountability landscape/Data/QUANTitative datasets/Complete datasets/Documentation unobtrustive/Independence/Independence - git") # Working directory

## Data
data <- read_excel("independence - complete -221220.xlsx", na = "99")
data <- data %>% mutate_if(is.numeric, as.factor)
data <- data %>%  column_to_rownames(var = "agency")

## Correlation plot
cor <- hetcor(data[c(1:8,11:17, 19)], use = "pairwise.complete.obs")
cor$correlations <- format(cor$correlations, digits = 1)
grid.table(cor$correlations)


## Chain 
post.samp <- MCMCordfactanal(~ ah.term + ah.selection + ah.quorum + ah.dismissal +
            ah.reappointment + mb.reappointment + ah.independence + ah.requirement +
            mb.fixed + mb.term + mb.reappointment + mb.independence + mb.requirement +
              agency.independence + agency.report + agency.program + agency.appeal +
              agency.forum + agency.composition + agency.discharge,
                             seed = 1090, factors = 1, data = data,
                             lambda.constraints = list(ah.independence = list(2,"+"),
                                                       agency.independence = list(2,"+")),
                             burnin = 100000, mcmc = 1000000, thin = 200,
                             verbose = TRUE,
                             store.lambda = TRUE, store.scores = TRUE, tune = 1)

save.image(file = "independence - 221220.RData")

## General info
varnames(post.samp) # Variable names
nchain(post.samp) # Number of chains
mcpar(post.samp) # Start iteration, end iteration and thinning interval

## Subset parameters
names <- attr(post.samp, which = "dimnames")
names <- names[[2]]
names

lambda.1 <- names[c(1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37)] # Lambda 1 (negative item difficulty parameter)
lambda.2 <- names[c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38)] # Lambda 2 (factor loadings)
gamma <- names[39:40] # Gamma 
agencies <- names[51:95] # Agencies

## Summary statistics
summary(post.samp[,1:38])

## Thin and burn-in

  # Acceptance rate
round(attr(post.samp, which = "accept.rates"), 2) # ???

  # Auto-correlation
acf1 <- acfplot(post.samp[,1:38],  thin = 1, ylab = "", xlab = "No thinning", par.settings = simpleTheme(col.line = "blue"),  par.strip.text = list(cex = 1),  aspect = 0.2, layout = c(5, 8))
acf5 <- acfplot(post.samp[,1:38],  thin = 5, ylab = "", xlab = "thin = 5", par.settings = simpleTheme(col.line = "red"), par.strip.text = list(cex = 1), aspect = 0.2, layout = c(5, 8))
acf15 <- acfplot(post.samp[,1:38], thin = 15, ylab = "", xlab = "thin = 15", par.settings = simpleTheme(col.line = "orange"), par.strip.text = list(cex = 1), aspect = 0.2, layout = c(5, 8))

acfplot <- plot_grid(acf1, acf5, acf15, nrow = 3)
acfplot

  # Raftery and Lewis
raftery <- raftery.diag(post.samp)  # I: Values of I much greater than 1 (I > 5) indicate high within-chain correlations and likely convergence failure and reparametrization is advised.
                                    # Nmin: required numbers of iterations if the samples were independent
                                    # N: required numbers of iterations for each variable
                                    # Minimal burn-in

sink("Raftery MCMC 23-12.txt")
print(raftery)
sink()

# Check burnin and thin
post.check <- mcmc(post.samp, start = 10000, thin = 20)
mcpar(post.check) # other

acfplot(post.check, ylab = "", par.settings = simpleTheme(col.line = "brown"), par.strip.text = list(cex = 0.5), aspect = 0.2, layout = c(10, 10))
plot(post.check, trace = TRUE, density = FALSE, smooth = FALSE)
plot(post.check, trace = FALSE, density = TRUE, smooth = TRUE)

# New run
post.thin2 <- MCMCordfactanal(~ ah.term + ah.selection + ah.quorum + ah.dismissal +
                               ah.reappointment + mb.reappointment +
                                ah.independence + ah.requirement +
                               mb.fixed + mb.independence + mb.requirement +
                               agency.independence + agency.report + agency.program + agency.appeal +
                               agency.forum + agency.discharge,
                             seed = 1090, factors = 1, data = data,
                             lambda.constraints = list(ah.independence = list(2,"+"),
                                                       agency.independence = list(2,"+")),
                             burnin = 10000, mcmc = 1000000, thin = 20,
                             verbose = TRUE,
                             store.lambda = TRUE, store.scores = TRUE, tune = 0.95)

post.thin3 <- MCMCordfactanal(~ ah.term + ah.selection + ah.quorum + ah.dismissal +
                                ah.reappointment + mb.reappointment + ah.independence + ah.requirement +
                                mb.fixed + mb.independence + mb.requirement +
                                agency.independence + agency.report + agency.program + agency.appeal +
                                agency.forum + agency.discharge,
                              seed = 1090, factors = 1, data = data,
                              lambda.constraints = list(ah.independence = list(2,"+"),
                                                        agency.independence = list(2,"+")),
                              burnin = 100000, mcmc = 1000000, thin = 200,
                              verbose = TRUE,
                              store.lambda = TRUE, store.scores = TRUE, tune = 0.85)

save.image(file = "independence - 241220.RData")

  # ACF
acfplot(post.thin, ylab = "", par.settings = simpleTheme(col.line = "brown"), par.strip.text = list(cex = 0.5), aspect = 0.2, layout = c(10, 10))

  # Trace
plot(post.thin, trace = TRUE, density = FALSE, smooth = FALSE)

 # Density
plot(post.thin, trace = FALSE, density = TRUE, smooth = TRUE)

## Diagnostics

  # Geweke 
geweke <- geweke.diag(post.thin) # Values of Z which fall in the extreme tails of N(0, 1) indicate that the chain has not yet converged

sink("Geweke MCMC 23-12.txt")
print(geweke)
sink()
geweke.plot(post.thin)

    # Raftery
raftery.thin <- raftery.diag(post.thin)   # I: Values of I much greater than 1 (I > 5) indicate high within-chain correlations and likely convergence failure and reparametrization is advised.
                                          # Nmin: required numbers of iterations if the samples were independent
                                          # N: required numbers of iterations for each variable
                                          # Minimal burn-in

sink("Raftery MCMC thin 23-12.txt")
print(raftery.thin)
sink()

  # Effect size
single.ess <- effectiveSize(post.thin) # For each variable

sink("ESS MCMC thin 24-12.txt")
print(single.ess)
sink()

multi.ess <- multiESS(post.thin)
min.ess <- minESS(p = 95, alpha = .05, eps = .05)
multi.ess - min.ess # Accounting for dependencies

  # Acceptance rate
attr(post.thin, which = "accept.rates") # ???
acceptance <- 1 - rejectionRate(post.samp)
acceptance
  
  # Heidelberger and Welch's convergence diagnostic
heidel <- heidel.diag(post.samp)  # Stationarity test
sink("Heidelberger MCMC thin 24-12.txt")
print(heidel)
sink() 

# Latent variable
post.list <- coda::as.mcmc.list(post.thin) # Create MCMC list

MCMCplot(post.list, 
         params = agencies, 
         rank = TRUE,
         xlab = 'Independence',
         guide_lines = TRUE,
         xlim = c(-3, 3),
         sz_labels = 1.5,
         sz_med = 2,
         sz_thick = 3,
         sz_thin = 1,
         sz_ax = 2,
         mar = c(1,1,1,1))

l1.plot <- MCMCsummary(post.list, 
            params = lambda.1,
            Rhat = FALSE,
            n.eff = FALSE,
            round = 2)

l2.plot <- MCMCsummary(post.list, 
                       params = lambda.2,
                       Rhat = FALSE,
                       n.eff = FALSE,
                       round = 2)

gamma.plot <- MCMCsummary(post.list, 
                          params = gamma,
                          Rhat = FALSE,
                          n.eff = FALSE,
                          round = 2)

tt3 <- ttheme_default(core = list(fg_params = list(hjust = 0, x = 0.1)),
                      rowhead = list(fg_params = list(hjust = 0, x = 0)))

grid.table(l1.plot,  theme = tt3)
grid.table(l2.plot,  theme = tt3)
grid.table(gamma.plot,  theme = tt3) # Plot? -> Koop & Hanretty
