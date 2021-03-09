## Before starting

  # Loading packages
library(coda)
library(lattice)
library(MCMCpack)
library(ggplot2)
library(mcmcse)
library(polycor)
library(dplyr)
library(gridExtra)
library(readxl)
library(ggmcmc)
library(stringr)
library(tibble)


  # Graphical settings
par(mar=c(1,1,1,1)) # Margins
options(scipen=999)
theme_update(text = element_text(size = 12, face = "bold"))

  # Set working directory
setwd("C:/Users/Thijs/surfdrive/Project A - Accountability landscape/Data/QUANTitative datasets/Complete datasets/Documentation unobtrustive/Independence")

  # Prepare data
data <- read_excel("independence - data.xlsx")
data <- data %>% mutate_if(is.numeric, as.ordered)
data <- data %>%  column_to_rownames(var = "agency")

## Descriptives and correlations

  # Check distribution
counts <- data %>% 
  lapply(table) %>% 
  lapply(as.data.frame)

counts %>%
  print()

  # Correlation matrix
cor <- hetcor(data, use = "pairwise.complete.obs")
grid.table(round(cor$correlations, 2))

## Note: The variables 'agency.discharge' and 'agency.forum' have strong
## negative correlations with the other variables and are therefore
## not included in the chain.

## Running the chain

  # Priors
l0.prior <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                  -0.08,
                  0.4,
                  0,
                  0.1,
                  0.66,
                  0.73,
                  0,
                  0.73,
                  0.79,
                  0,
                  0.77,
                  0.08,
                  0,
                  0) 
dim(l0.prior) <- c(14,2) # Adapted from from Hanretty &  Koop (2012).

  # Chain 
sample.post <- MCMCordfactanal(~ ah.selection +
                                    ah.term +
                                    ah.quorum +  
                                    ah.dismissal +  
                                    ah.independence +
                                    ah.reappointment +
                                    ah.requirement +     
                                    mb.fixed +          
                                    mb.independence +  
                                    mb.requirement +
                                    agency.independence + 
                                    agency.report +      
                                    agency.program +
                                    agency.appeal,
                                  lambda.constraints = 
                                    list(mb.fixed = list(2,"+")),
                                  seed = 1090,
                                  factors = 1,
                                  data = data,
                                  burnin = 10000,
                                  mcmc = 1000000,
                                  thin = 5,
                                  verbose = 10000,
                                  L0 = 1,
                                  store.lambda = TRUE,
                                  store.scores = TRUE,
                                  tune = 4,
                                  l0.prior = l0.prior) # mcmc: DO NOT CHANGE
save.image()

## Preparing for analysis: extracting labels for later

  # Dimension names
dimnames <- attr(sample.post, which = "dimnames")
names <- dimnames[[2]]
names 

  # Filter
lambda <-  grep("^Lambda", names, value = TRUE) 
lambda.1 <- grep("1$", lambda, value = TRUE)  # Lambda 1 (negative item difficulty parameter)
lambda.2 <- grep("2$", lambda, value = TRUE) # Lambda 2 (factor loadings or the item discrimination parameters)
gamma <- grep("gamma", names, value = TRUE)  # Gamma 
phi <-  grep("phi", names, value = TRUE)

  # Create labels
lambda.1.label <- lambda.1 %>% str_sub(start = 7, end = -3)%>% paste0(" (lambda 1)")
lambda.2.label <- lambda.2 %>% str_sub(start = 7, end = -3) %>% paste0(" (lambda 2)")
gamma.label <- gamma  %>% str_sub(start = 8) %>% paste0(" (gamma)")
agencies.label  <- phi %>% str_sub(start = 5, end = -3)

  # Create dataframe
lambda.1.df <- data.frame(
  Parameter = lambda.1,
  Label = lambda.1.label)

lambda.2.df <- data.frame(
  Parameter = lambda.2,
  Label = lambda.2.label)

gamma.df <- data.frame(
  Parameter = gamma,
  Label = gamma.label)

agencies.df <- data.frame(
  Parameter = agencies,
  Label = agencies.label)

total.df <- rbind(lambda.1.df, lambda.2.df, gamma.df, agencies.df)

 # Convert mcmc object to tibble 
sample.post.converted <- ggs(as.mcmc.list(sample.post), par_labels = total.df)
sample.post.agencies <- ggs(as.mcmc.list(sample.post), par_labels = agencies.df, family = "^phi")

## Key information of the chain

  # Summary
summary <- summary(sample.post)
scores.df <- summary[["statistics"]] %>%
  as.data.frame() 

scores <- scores.df %>%
  mutate(Agency = row.names(scores.df)) %>%
  filter(str_detect(Agency, "^phi"))  %>%
   mutate(Agency = str_remove(Agency, ".2"),
          Agency = str_remove(Agency, "phi."))

write.csv(scores, "scores independence.csv", row.names = FALSE)

# Acceptance
attr(sample.post, which = "accept.rates")

  # Trace
plot(sample.post, trace = TRUE, density = FALSE, smooth = FALSE)

  # Density
plot(sample.post, trace = FALSE, density = TRUE, smooth = TRUE)

  # ACF
acfplot(sample.post,
        ylab = "",
        par.settings = simpleTheme(col.line = "brown"),
        par.strip.text = list(cex = 0.4),
        aspect = 0.5,
        layout = c(11, 7))

autocorr.diag(sample.post)

autocorr.plot(sample.post, lag.max = 10)

  # Running means
ggs_running(sample.post.converted) +
  facet_wrap(~ Parameter, ncol = 11, scales = "free") +
  geom_line(size = 2, colour = "black") + 
  theme(legend.position = "none")

## Diagnostics

  # Raftery
raftery <- raftery.diag(sample.post)
raftery

sink("Raftery  Independence.txt")
print(raftery)
sink()

  # Geweke 
geweke <- geweke.diag(sample.post) 
geweke[["z"]] %>% as_tibble() %>% filter(value < -1.96 | value > 1.96) # Parameters outside expected range

sink("Geweke  Independence.txt")
print(geweke)
sink()

geweke.plot(sample.post)

  # Effective size
single.ess <- effectiveSize(sample.post) # For each variable
single.ess

sink("ESS  Independence.txt")
print(single.ess) # The larger the better. ESSs < 100 is really bad and > 200 sufficient. On the other hand chasing ESSs > 10000 may be a waste of computational resources
sink()

multi.ess <- multiESS(sample.post)
min.ess <- minESS(p = 95, alpha = .05, eps = .05)
multi.ess - min.ess # Accounting for dependencies

  # Heidelberger and Welch's convergence diagnostic
heidel <- heidel.diag(sample.post)  # Stationarity test
sink("Heidelberger Independence.txt.txt")
print(heidel)
sink() 

## Latent variable

  # Plot
ggs_caterpillar(sample.post.agencies) +
  scale_x_continuous(name =  "θ Formal independence",
                     breaks = seq(-3, 3, 1),  
                     limits = c(-3.25, 3.25)) +
  labs(y = "") +
  theme(axis.title.x  = element_text(size = 20, face = "italic"),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "grey50", size = 0.1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black",
                                    fill = NA,
                                    size = 0.4))

## References
  # Hanretty, C. & Koop, C, (2012). Measuring the formal independence of regulatory agencies,
  # Journal of European Public Policy, 19:2, 198-216, DOI: 10.1080/13501763.2011.607357
