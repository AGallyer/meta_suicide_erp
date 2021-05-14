library(tidyverse)
library(metafor)
library(readxl)
library(compute.es)
library(PRISMAstatement)
library(metaviz)
library(metameta)
library(ggthemes)
library(cowplot)
#PRISMA statement
prisma(found = 298,
       found_other = 4,
       no_dupes = 219, 
       screened = 219, 
       screen_exclusions = 163, 
       full_text = 56,
       full_text_exclusions = 29, 
       qualitative = 27, 
       quantitative = 27)

#Import data
data <- read_excel("meta_spreadsheet.xlsx", sheet = "effect_size_coding")

meta_data <- data %>% 
  dplyr::select(Study_ID, effect_id, Description, 
                type, g, var_g, everything()) %>% 
  mutate("abs_g" = abs(g), 
         "lbci" = g - (1.96*sqrt(var_g)), 
         "ubci" = g + (1.96*sqrt(var_g)))


meta_data$scoring <- recode_factor(as.factor(meta_data$scoring), raw = 0, 
                                   difference = 1, slope = 1)


#Separate ERP Studies ---------------------------------------------------
# Create separate data for each type and outcome
# LDAEP
ldaep_si_data <- meta_data %>%
  filter(type == "LDAEP" & Description == "ideation")# 3 effect sizes 2 studies. Meta analyze

ldaep_sa_data <- meta_data %>%
  filter(type == "LDAEP" & Description == "attempt")# 12 effect sizes, 5 studies. Meta analyze

ldaep_risk_data <- meta_data %>%
  filter(type == "LDAEP" & Description == "risk")# 3 effect sizes, 3 studies. Meta analyze

ldaep_sasi_data <- meta_data %>% 
  filter(type == "LDAEP" & suicide_group == 'SA' & control_group == 'SI')# 11 effect sizes, 4 studies meta analyze

## Analyses
### LDAEP SI Analysis
ldaep_estimate <- rma.mv(yi = g ~ 1,
                         V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID),
                         data = ldaep_si_data,
                         method = "REML",  
                         slab = paste(paste(paste(Authors, Year, sep=", ("), ")", 
                                                              sep = ""), effect_id, sep = " "))

ldaep_robust <- robust(ldaep_estimate, cluster = ldaep_si_data$Study_ID,
                       adjust = TRUE)
summary(ldaep_robust)# g = - 0.2720, p = 0.6803

#### Compute I^2 for this model
W <- diag(1/ldaep_si_data$var_g)
X <- model.matrix(ldaep_robust)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(ldaep_robust$sigma2) / (sum(ldaep_robust$sigma2) + (ldaep_robust$k-ldaep_robust$p)/sum(diag(P)))

#### Create forest plot
fig2a <- viz_forest(ldaep_robust, type = "study_only", variant = "classic", study_labels = ldaep_estimate[["slab"]], 
                    col = "#1e89c6", 
           text_size = 5, x_limit = c(-5.0, 4), x_breaks = c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4), 
           xlab = expression("Hedges' " * italic("g"))) + ggtitle("LDAEP")

summarydata <- data.frame("x.diamond" = c(ldaep_robust$ci.lb,
                                         coef(ldaep_robust),
                                         ldaep_robust$ci.ub,
                                          coef(ldaep_robust)),
                          "y.diamond" = c(nrow(ldaep_si_data) + 1,
                                          nrow(ldaep_si_data) + 1.3,
                                          nrow(ldaep_si_data) + 1,
                                          nrow(ldaep_si_data) + 0.7),
                          "diamond_group" = rep(1, times = 4))

ylabels <- c("Summary", ldaep_robust$slab)

fig2a <- fig2a + geom_polygon(data = summarydata, aes(x = x.diamond, y = y.diamond, group = diamond_group), 
                     color= "black", fill = "#1e89c6", 
                     size = 0.1) + scale_y_continuous(name = "",
                                                      breaks = c(4, 3, 2, 1), labels = ylabels)
### LDAEP SA Analysis
ldaep_estimate <- rma.mv(yi = g ~ 1,
                         V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID),
                         data = ldaep_sa_data,
                         method = "REML", 
                         slab = paste(paste(paste(Authors, Year, sep=", ("), ")", 
                                            sep = ""), effect_id, sep = " "))

ldaep_robust <- robust(ldaep_estimate, cluster = ldaep_sa_data$Study_ID,
                       adjust = TRUE)
summary(ldaep_robust)# g = 0.021, p = 0.946

#### Compute I^2 for this model
W <- diag(1/ldaep_sa_data$var_g)
X <- model.matrix(ldaep_robust)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(ldaep_robust$sigma2) / (sum(ldaep_robust$sigma2) + (ldaep_robust$k-ldaep_robust$p)/sum(diag(P)))

#### Create forest plot
fig3a <- viz_forest(ldaep_robust, type = "study_only", variant = "classic", study_labels = ldaep_estimate[["slab"]], 
                    col = "#1e89c6", xlab = expression("Hedges' " * italic("g")), 
                    text_size = 5) + ggtitle("LDAEP")

summarydata <- data.frame("x.diamond" = c(ldaep_robust$ci.lb,
                                          coef(ldaep_robust),
                                          ldaep_robust$ci.ub,
                                          coef(ldaep_robust)),
                          "y.diamond" = c(nrow(ldaep_sa_data) + 1,
                                          nrow(ldaep_sa_data) + 1.3,
                                          nrow(ldaep_sa_data) + 1,
                                          nrow(ldaep_sa_data) + 0.7),
                          "diamond_group" = rep(1, times = 4))

ylabels <- c("Summary", ldaep_robust$slab)
nrows <- nrow(ldaep_sa_data)
ybreaklim <- nrows + 1
ybreaks <- c(sort(1:ybreaklim, decreasing = TRUE))

fig3a <- fig3a + geom_polygon(data = summarydata, aes(x = x.diamond, y = y.diamond, group = diamond_group), 
                              color= "black", fill = "#1e89c6", 
                              size = 0.1) + scale_y_continuous(name = "",
                                                               breaks = ybreaks, labels = ylabels)



fig3a

#### LDAEP Risk Publication Bias
ldaep_weights <- weights(ldaep_robust)

ldaep_sa_data$weights <- ldaep_weights

ldaep_sa_data <- mutate(ldaep_sa_data, "sqrt_weights" = sqrt(weights))

ldaep_estimate <- rma.mv(yi = abs_g ~ sqrt_weights, 
                         V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID), 
                         data = ldaep_sa_data,
                         method = "REML")

ldaep_eggers <- robust(ldaep_estimate, cluster = ldaep_sa_data$Study_ID, 
                       adjust = TRUE)

summary(ldaep_eggers)# not significant

### LDAEP Risk Analysis
ldaep_estimate <- rma.mv(yi = g ~ 1,
                         V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID),
                         data = ldaep_risk_data,
                         method = "REML", 
                         slab = paste(paste(paste(Authors, Year, sep=", ("), ")", 
                                            sep = ""), effect_id, sep = " "))

ldaep_robust <- robust(ldaep_estimate, cluster = ldaep_risk_data$Study_ID,
                       adjust = TRUE)
summary(ldaep_robust)# g = 0.175, p = 0.552

#### Compute I^2 for this model
W <- diag(1/ldaep_risk_data$var_g)
X <- model.matrix(ldaep_robust)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(ldaep_robust$sigma2) / (sum(ldaep_robust$sigma2) + (ldaep_robust$k-ldaep_robust$p)/sum(diag(P)))

#### Create forest plot
fig4a <- viz_forest(ldaep_robust, type = "study_only", variant = "classic", study_labels = ldaep_estimate[["slab"]], 
                    col = "#1e89c6", xlab = expression("Hedges' " * italic("g")), 
                    text_size = 5, x_limit = c(-2.0, 2.0)) + ggtitle("LDAEP")

summarydata <- data.frame("x.diamond" = c(ldaep_robust$ci.lb,
                                          coef(ldaep_robust),
                                          ldaep_robust$ci.ub,
                                          coef(ldaep_robust)),
                          "y.diamond" = c(nrow(ldaep_risk_data) + 1,
                                          nrow(ldaep_risk_data) + 1.3,
                                          nrow(ldaep_risk_data) + 1,
                                          nrow(ldaep_risk_data) + 0.7),
                          "diamond_group" = rep(1, times = 4))

ylabels <- c("Summary", ldaep_robust$slab)
nrows <- nrow(ldaep_risk_data)
ybreaklim <- nrows + 1
ybreaks <- c(sort(1:ybreaklim, decreasing = TRUE))

fig4a <- fig4a + geom_polygon(data = summarydata, aes(x = x.diamond, y = y.diamond, group = diamond_group), 
                              color= "black", fill = "#1e89c6", 
                              size = 0.1) + scale_y_continuous(name = "",
                                                               breaks = ybreaks, labels = ylabels)

fig4a

### LDAEP SA vs. SI analysis
ldaep_estimate <- rma.mv(yi = g ~ 1,
                         V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID),
                         data = ldaep_sasi_data,
                         method = "REML",  
                         slab = paste(paste(paste(Authors, Year, sep=", ("), ")", 
                                            sep = ""), effect_id, sep = " "))

ldaep_robust <- robust(ldaep_estimate, cluster = ldaep_sasi_data$Study_ID,
                       adjust = TRUE)
summary(ldaep_robust)# g = 0.168, p = 0.648

#### Compute I^2 for this model
W <- diag(1/ldaep_sasi_data$var_g)
X <- model.matrix(ldaep_robust)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(ldaep_robust$sigma2) / (sum(ldaep_robust$sigma2) + (ldaep_robust$k-ldaep_robust$p)/sum(diag(P)))
#69.17%
#### Create forest plot
fig5a <- viz_forest(ldaep_robust, type = "study_only", variant = "classic", study_labels = ldaep_estimate[["slab"]], 
                    col = "#1e89c6", xlab = expression("Hedges' " * italic("g")), 
                    text_size = 5) + ggtitle("LDAEP")

summarydata <- data.frame("x.diamond" = c(ldaep_robust$ci.lb,
                                          coef(ldaep_robust),
                                          ldaep_robust$ci.ub,
                                          coef(ldaep_robust)),
                          "y.diamond" = c(nrow(ldaep_sasi_data) + 1,
                                          nrow(ldaep_sasi_data) + 1.3,
                                          nrow(ldaep_sasi_data) + 1,
                                          nrow(ldaep_sasi_data) + 0.7),
                          "diamond_group" = rep(1, times = 4))

ylabels <- c("Summary", ldaep_robust$slab)
nrows <- nrow(ldaep_sasi_data)
ybreaklim <- nrows + 1
ybreaks <- c(sort(1:ybreaklim, decreasing = TRUE))

fig5a <- fig5a + geom_polygon(data = summarydata, aes(x = x.diamond, y = y.diamond, group = diamond_group), 
                              color= "black", fill = "#1e89c6", 
                              size = 0.1) + scale_y_continuous(name = "",
                                                               breaks = ybreaks, labels = ylabels)



fig5a


#### LDAEP SA vs. SI Publication Bias
ldaep_weights <- weights(ldaep_robust)

ldaep_sasi_data$weights <- ldaep_weights

ldaep_sasi_data <- mutate(ldaep_sasi_data, "sqrt_weights" = sqrt(weights))

ldaep_estimate <- rma.mv(yi = abs_g ~ sqrt_weights, 
                         V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID), 
                         data = ldaep_sasi_data,
                         method = "REML")

ldaep_eggers <- robust(ldaep_estimate, cluster = ldaep_sasi_data$Study_ID, 
                       adjust = TRUE)

summary(ldaep_eggers)# not significant

# LPP
lpp_si_data <- meta_data %>%
  filter(type == "LPP" & Description == "ideation")# 34 effect sizes, 3 studies. Meta analyze

lpp_sa_data <- meta_data %>%
  filter(type == "LPP" & Description == "attempt")# 1 study, 3 effect sizes. Meta analyze

lpp_risk_data <- meta_data %>%
  filter(type == "LPP" & Description == "risk")# 9 effect sizes, 2 studies. Meta analyze

lpp_sasi_data <- meta_data %>% 
  filter(type == "LPP" & suicide_group == 'SA' & control_group == 'SI')# Same effects as SA LPP dataset


## Analyses
###LPP SI Analysis
lpp_estimate <- rma.mv(yi = g ~ 1,
                       V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID),
                       data = lpp_si_data,
                       method = "REML", slab = paste(paste(paste(Authors, Year, sep=", ("), ")", 
                                                           sep = ""), effect_id, sep = " "))

lpp_robust <- robust(lpp_estimate, cluster = lpp_si_data$Study_ID,
                     adjust = TRUE)
summary(lpp_robust)# g = -0.053, p = 0.407

#### Compute I^2 for this model
W <- diag(1/lpp_si_data$var_g)
X <- model.matrix(lpp_robust)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(lpp_robust$sigma2) / (sum(lpp_robust$sigma2) + (lpp_robust$k-lpp_robust$p)/sum(diag(P)))
#10.73%
#### Create forest plot
fig2b <- viz_forest(lpp_robust, type = "study_only", variant = "classic", study_labels = ldaep_estimate[["slab"]], 
                    col = "#1e89c6", xlab = expression("Hedges' " * italic("g")), 
                    text_size = 5) + ggtitle("LPP")

summarydata <- data.frame("x.diamond" = c(lpp_robust$ci.lb,
                                          coef(lpp_robust),
                                          lpp_robust$ci.ub,
                                          coef(lpp_robust)),
                          "y.diamond" = c(nrow(lpp_si_data) + 1,
                                          nrow(lpp_si_data) + 1.3,
                                          nrow(lpp_si_data) + 1,
                                          nrow(lpp_si_data) + 0.7),
                          "diamond_group" = rep(1, times = 4))

ylabels <- c("Summary", lpp_robust$slab)
nrows <- nrow(lpp_si_data)
ybreaklim <- nrows + 1
ybreaks <- c(sort(1:ybreaklim, decreasing = TRUE))

fig2b <- fig2b + geom_polygon(data = summarydata, aes(x = x.diamond, y = y.diamond, group = diamond_group), 
                              color= "black", fill = "#1e89c6", 
                              size = 0.1) + scale_y_continuous(name = "",
                                                               breaks = ybreaks, labels = ylabels)



fig2b



#### LPP SI publication bias
lpp_weights <- weights(lpp_robust)

lpp_si_data$weights <- lpp_weights

lpp_si_data <- mutate(lpp_si_data, "sqrt_weights" = sqrt(weights))

lpp_estimate <- rma.mv(yi = abs_g ~ sqrt_weights, 
                         V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID), 
                         data = lpp_si_data,
                         method = "REML")

lpp_eggers <- robust(lpp_estimate, cluster = lpp_si_data$Study_ID, 
                       adjust = TRUE)

summary(lpp_eggers)# not significant

### LPP Risk Analysis
lpp_estimate <- rma.mv(yi = g ~ 1,
                       V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID),
                       data = lpp_risk_data,
                       method = "REML", slab = paste(paste(paste(Authors, Year, sep=", ("), ")", 
                                                           sep = ""), effect_id, sep = " "))


lpp_robust <- robust(lpp_estimate, cluster = lpp_risk_data$Study_ID,
                     adjust = TRUE)
summary(lpp_robust)# g = -0.290, p = 0.157

W <- diag(1/lpp_risk_data$var_g)
X <- model.matrix(lpp_robust)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(lpp_robust$sigma2) / (sum(lpp_robust$sigma2) + (lpp_robust$k-lpp_robust$p)/sum(diag(P)))
#9.15%

#### Create Forest Plot
fig4b <- viz_forest(lpp_robust, type = "study_only", variant = "classic", study_labels = lpp_estimate[["slab"]], 
                    col = "#1e89c6", xlab = expression("Hedges' " * italic("g")), 
                    text_size = 5, x_limit = c(-1.5, 1.0)) + ggtitle("LPP")

summarydata <- data.frame("x.diamond" = c(lpp_robust$ci.lb,
                                          coef(lpp_robust),
                                          lpp_robust$ci.ub,
                                          coef(lpp_robust)),
                          "y.diamond" = c(nrow(lpp_risk_data) + 1,
                                          nrow(lpp_risk_data) + 1.3,
                                          nrow(lpp_risk_data) + 1,
                                          nrow(lpp_risk_data) + 0.7),
                          "diamond_group" = rep(1, times = 4))

ylabels <- c("Summary", lpp_robust$slab)
nrows <- nrow(lpp_risk_data)
ybreaklim <- nrows + 1
ybreaks <- c(sort(1:ybreaklim, decreasing = TRUE))

fig4b <- fig4b + geom_polygon(data = summarydata, aes(x = x.diamond, y = y.diamond, group = diamond_group), 
                              color= "black", fill = "#1e89c6", 
                              size = 0.1) + scale_y_continuous(name = "",
                                                               breaks = ybreaks, labels = ylabels)



fig4b

# P3
p3_si_data <- meta_data %>%
  filter(type == "P3" & Description == "ideation")# 53 effect sizes, 3 studies. meta analyze

p3_sa_data <- meta_data %>%
  filter(type == "P3" & Description == "attempt")# 15 effect sizes, 6 studies. Meta analyze

p3_risk_data <- meta_data %>%
  filter(type == "P3" & Description == "risk")# 28 effect sizes, 3 studies. meta analyze

p3_sasi_data <- meta_data %>% 
  filter(type == "P3" & suicide_group == 'SA' & control_group == 'SI')# 10 effect sizes, 3 studies. meta analyze


## Analyses
###P3 SI Analysis
p3_estimate <- rma.mv(yi = g ~ 1,
                      V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID),
                      data = p3_si_data,
                      method = "REML", slab = paste(paste(paste(Authors, Year, sep=", ("), ")", 
                                                          sep = ""), effect_id, sep = " "))

p3_robust <- robust(p3_estimate, cluster = p3_si_data$Study_ID,
                    adjust = TRUE)
summary(p3_robust)# g = 0.067, p = 0.827

W <- diag(1/p3_si_data$var_g)
X <- model.matrix(p3_robust)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(p3_robust$sigma2) / (sum(p3_robust$sigma2) + (p3_robust$k-p3_robust$p)/sum(diag(P)))
#79.75

#### Create Forest Plot
fig2c <- viz_forest(p3_robust, type = "study_only", variant = "classic", study_labels = p3_estimate[["slab"]], 
                    col = "#1e89c6", xlab = expression("Hedges' " * italic("g")), 
                    text_size = 3.5, x_limit = c(-2, 2)) + ggtitle("P300")

summarydata <- data.frame("x.diamond" = c(p3_robust$ci.lb,
                                          coef(p3_robust),
                                          p3_robust$ci.ub,
                                          coef(p3_robust)),
                          "y.diamond" = c(nrow(p3_si_data) + 1,
                                          nrow(p3_si_data) + 1.3,
                                          nrow(p3_si_data) + 1,
                                          nrow(p3_si_data) + 0.7),
                          "diamond_group" = rep(1, times = 4))

ylabels <- c("Summary", p3_robust$slab)
nrows <- nrow(p3_si_data)
ybreaklim <- nrows + 1
ybreaks <- c(sort(1:ybreaklim, decreasing = TRUE))

fig2c <- fig2c + geom_polygon(data = summarydata, aes(x = x.diamond, y = y.diamond, group = diamond_group), 
                              color= "black", fill = "#1e89c6", 
                              size = 0.1) + scale_y_continuous(name = "",
                                                               breaks = ybreaks, labels = ylabels)



fig2c


#### P300 SI publication bias
p3_weights <- weights(p3_robust)

p3_si_data$weights <- p3_weights

p3_si_data <- mutate(p3_si_data, "sqrt_weights" = sqrt(weights))

p3_estimate <- rma.mv(yi = abs_g ~ sqrt_weights, 
                       V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID), 
                       data = p3_si_data,
                       method = "REML")

p3_eggers <- robust(p3_estimate, cluster = p3_si_data$Study_ID, 
                     adjust = TRUE)

summary(p3_eggers)# not significant


###P3 SA Analysis
p3_estimate <- rma.mv(yi = g ~ 1,
                      V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID),
                      data = p3_sa_data,
                      method = "REML", slab = paste(paste(paste(Authors, Year, sep=", ("), ")", 
                                                          sep = ""), effect_id, sep = " "))

p3_robust <- robust(p3_estimate, cluster = p3_sa_data$Study_ID,
                    adjust = TRUE)
summary(p3_robust)# g = 0.067, p = 0.827

W <- diag(1/p3_sa_data$var_g)
X <- model.matrix(p3_robust)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(p3_robust$sigma2) / (sum(p3_robust$sigma2) + (p3_robust$k-p3_robust$p)/sum(diag(P)))
#93.97

#### Create Forest Plot
fig3b <- viz_forest(p3_robust, type = "study_only", variant = "classic", study_labels = p3_estimate[["slab"]], 
                    col = "#1e89c6", xlab = expression("Hedges' " * italic("g")), 
                    text_size = 5) + ggtitle("P300")

summarydata <- data.frame("x.diamond" = c(p3_robust$ci.lb,
                                          coef(p3_robust),
                                          p3_robust$ci.ub,
                                          coef(p3_robust)),
                          "y.diamond" = c(nrow(p3_sa_data) + 1,
                                          nrow(p3_sa_data) + 1.3,
                                          nrow(p3_sa_data) + 1,
                                          nrow(p3_sa_data) + 0.7),
                          "diamond_group" = rep(1, times = 4))

ylabels <- c("Summary", p3_robust$slab)
nrows <- nrow(p3_sa_data)
ybreaklim <- nrows + 1
ybreaks <- c(sort(1:ybreaklim, decreasing = TRUE))

fig3b <- fig3b + geom_polygon(data = summarydata, aes(x = x.diamond, y = y.diamond, group = diamond_group), 
                              color= "black", fill = "#1e89c6", 
                              size = 0.1) + scale_y_continuous(name = "",
                                                               breaks = ybreaks, labels = ylabels)



fig3b


#### P300 SA publication bias
p3_weights <- weights(p3_robust)

p3_sa_data$weights <- p3_weights

p3_sa_data <- mutate(p3_sa_data, "sqrt_weights" = sqrt(weights))

p3_estimate <- rma.mv(yi = abs_g ~ sqrt_weights, 
                      V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID), 
                      data = p3_sa_data,
                      method = "REML")

p3_eggers <- robust(p3_estimate, cluster = p3_sa_data$Study_ID, 
                    adjust = TRUE)

summary(p3_eggers)# not significant

### P3 Risk Analysis
p3_estimate <- rma.mv(yi = g ~ 1,
                      V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID),
                      data = p3_risk_data,
                      method = "REML", slab = paste(paste(paste(Authors, Year, sep=", ("), ")", 
                                                          sep = ""), effect_id, sep = " "))

p3_robust <- robust(p3_estimate, cluster = p3_risk_data$Study_ID,
                    adjust = TRUE)
summary(p3_robust)# g = 0.2497, p = .3218

W <- diag(1/p3_risk_data$var_g)
X <- model.matrix(p3_robust)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(p3_robust$sigma2) / (sum(p3_robust$sigma2) + (p3_robust$k-p3_robust$p)/sum(diag(P)))
#65.09

#### Create Forest Plot
fig4c <- viz_forest(p3_robust, type = "study_only", variant = "classic", study_labels = p3_estimate[["slab"]], 
                    col = "#1e89c6", xlab = expression("Hedges' " * italic("g")), 
                    text_size = 5) + ggtitle("P300")

summarydata <- data.frame("x.diamond" = c(p3_robust$ci.lb,
                                          coef(p3_robust),
                                          p3_robust$ci.ub,
                                          coef(p3_robust)),
                          "y.diamond" = c(nrow(p3_risk_data) + 1,
                                          nrow(p3_risk_data) + 1.3,
                                          nrow(p3_risk_data) + 1,
                                          nrow(p3_risk_data) + 0.7),
                          "diamond_group" = rep(1, times = 4))

ylabels <- c("Summary", p3_robust$slab)
nrows <- nrow(p3_risk_data)
ybreaklim <- nrows + 1
ybreaks <- c(sort(1:ybreaklim, decreasing = TRUE))

fig4c <- fig4c + geom_polygon(data = summarydata, aes(x = x.diamond, y = y.diamond, group = diamond_group), 
                              color= "black", fill = "#1e89c6", 
                              size = 0.1) + scale_y_continuous(name = "",
                                                               breaks = ybreaks, labels = ylabels)



fig4c


#### P300 risk publication bias
p3_weights <- weights(p3_robust)

p3_risk_data$weights <- p3_weights

p3_risk_data <- mutate(p3_risk_data, "sqrt_weights" = sqrt(weights))

p3_estimate <- rma.mv(yi = abs_g ~ sqrt_weights, 
                      V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID), 
                      data = p3_risk_data,
                      method = "REML")

p3_eggers <- robust(p3_estimate, cluster = p3_risk_data$Study_ID, 
                    adjust = TRUE)

summary(p3_eggers)# not significant


###P3 SA vs. SI Analysis
p3_estimate <- rma.mv(yi = g ~ 1,
                      V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID),
                      data = p3_sasi_data,
                      method = "REML", slab = paste(paste(paste(Authors, Year, sep=", ("), ")", 
                                                          sep = ""), effect_id, sep = " "))

p3_robust <- robust(p3_estimate, cluster = p3_sasi_data$Study_ID,
                    adjust = TRUE)
summary(p3_robust)# g = -0.40, p = 0.6422

W <- diag(1/p3_sasi_data$var_g)
X <- model.matrix(p3_robust)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(p3_robust$sigma2) / (sum(p3_robust$sigma2) + (p3_robust$k-p3_robust$p)/sum(diag(P)))
#93.06

#### Create Forest Plot
fig5b <- viz_forest(p3_robust, type = "study_only", variant = "classic", study_labels = p3_estimate[["slab"]], 
                    col = "#1e89c6", xlab = expression("Hedges' " * italic("g")), 
                    text_size = 5, x_limit = c(-4, 3)) + ggtitle("P300")

summarydata <- data.frame("x.diamond" = c(p3_robust$ci.lb,
                                          coef(p3_robust),
                                          p3_robust$ci.ub,
                                          coef(p3_robust)),
                          "y.diamond" = c(nrow(p3_sasi_data) + 1,
                                          nrow(p3_sasi_data) + 1.3,
                                          nrow(p3_sasi_data) + 1,
                                          nrow(p3_sasi_data) + 0.7),
                          "diamond_group" = rep(1, times = 4))

ylabels <- c("Summary", p3_robust$slab)
nrows <- nrow(p3_sasi_data)
ybreaklim <- nrows + 1
ybreaks <- c(sort(1:ybreaklim, decreasing = TRUE))

fig5b <- fig5b + geom_polygon(data = summarydata, aes(x = x.diamond, y = y.diamond, group = diamond_group), 
                              color= "black", fill = "#1e89c6", 
                              size = 0.1) + scale_y_continuous(name = "",
                                                               breaks = ybreaks, labels = ylabels)



fig5b


#### P300 SA vs. SI -only publication bias
p3_weights <- weights(p3_robust)

p3_sasi_data$weights <- p3_weights

p3_sasi_data <- mutate(p3_sasi_data, "sqrt_weights" = sqrt(weights))

p3_estimate <- rma.mv(yi = abs_g ~ sqrt_weights, 
                      V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID), 
                      data = p3_sasi_data,
                      method = "REML")

p3_eggers <- robust(p3_estimate, cluster = p3_sasi_data$Study_ID, 
                    adjust = TRUE)

summary(p3_eggers)# Significant publication bias


# RewP
rewp_si_data <- meta_data %>%
  filter(type == "RewP" & Description == "ideation")# 12 effect sizes, 3 studies. Meta analyze

rewp_sa_data <- meta_data %>%
  filter(type == "RewP" & Description == "attempt")# 5 effect sizes, 2 studies. Meta analyze

rewp_risk_data <- meta_data %>%
  filter(type == "RewP" & Description == "risk")# 3 effect sizes, 1 study. Narrative review only 

rewp_sasi_data <- meta_data %>% 
  filter(type == "RewP" & suicide_group == 'SA' & control_group == 'SI')# 0 effect sizes


## Analyses
### RewP SI Analysis
rewp_estimate <- rma.mv(yi = g ~ 1,
                      V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID),
                      data = rewp_si_data,
                      method = "REML", slab = paste(paste(paste(Authors, Year, sep=", ("), ")", 
                                                          sep = ""), effect_id, sep = " "))

rewp_robust <- robust(rewp_estimate, cluster = rewp_si_data$Study_ID,
                    adjust = TRUE)
summary(rewp_robust)# g = -0.06, p = .0406

W <- diag(1/rewp_si_data$var_g)
X <- model.matrix(rewp_robust)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(rewp_robust$sigma2) / (sum(rewp_robust$sigma2) + (rewp_robust$k-rewp_robust$p)/sum(diag(P)))
#4.29%

#### Create Forest Plot
fig2d <- viz_forest(rewp_robust, type = "study_only", variant = "classic", study_labels = rewp_estimate[["slab"]], 
                    col = "#1e89c6", xlab = expression("Hedges' " * italic("g")), 
                    text_size = 5) + ggtitle("RewP")

summarydata <- data.frame("x.diamond" = c(rewp_robust$ci.lb,
                                          coef(rewp_robust),
                                          rewp_robust$ci.ub,
                                          coef(rewp_robust)),
                          "y.diamond" = c(nrow(rewp_si_data) + 1,
                                          nrow(rewp_si_data) + 1.3,
                                          nrow(rewp_si_data) + 1,
                                          nrow(rewp_si_data) + 0.7),
                          "diamond_group" = rep(1, times = 4))

ylabels <- c("Summary", rewp_robust$slab)
nrows <- nrow(rewp_si_data)
ybreaklim <- nrows + 1
ybreaks <- c(sort(1:ybreaklim, decreasing = TRUE))

fig2d <- fig2d + geom_polygon(data = summarydata, aes(x = x.diamond, y = y.diamond, group = diamond_group), 
                              color= "black", fill = "#1e89c6", 
                              size = 0.1) + scale_y_continuous(name = "",
                                                               breaks = ybreaks, labels = ylabels)



fig2d


#### RewP SI publication bias
rewp_weights <- weights(rewp_robust)

rewp_si_data$weights <- rewp_weights

rewp_si_data <- mutate(rewp_si_data, "sqrt_weights" = sqrt(weights))

rewp_estimate <- rma.mv(yi = abs_g ~ sqrt_weights, 
                      V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID), 
                      data = rewp_si_data,
                      method = "REML")

rewp_eggers <- robust(rewp_estimate, cluster = rewp_si_data$Study_ID, 
                    adjust = TRUE)

summary(rewp_eggers)# not significant

### RewP SA Analysis
rewp_estimate <- rma.mv(yi = g ~ 1,
                        V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID),
                        data = rewp_sa_data,
                        method = "REML", slab = paste(paste(paste(Authors, Year, sep=", ("), ")", 
                                                            sep = ""), effect_id, sep = " "))

rewp_robust <- robust(rewp_estimate, cluster = rewp_sa_data$Study_ID,
                      adjust = TRUE)
summary(rewp_robust)# g = -0.12, p = .5449

W <- diag(1/rewp_sa_data$var_g)
X <- model.matrix(rewp_robust)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(rewp_robust$sigma2) / (sum(rewp_robust$sigma2) + (rewp_robust$k-rewp_robust$p)/sum(diag(P)))
#12.33%

#### Create Forest Plot
fig3c <- viz_forest(rewp_robust, type = "study_only", variant = "classic", study_labels = rewp_estimate[["slab"]], 
                    col = "#1e89c6", xlab = expression("Hedges' " * italic("g")), 
                    text_size = 5, x_limit = c(-2.25, 2)) + ggtitle("RewP")

summarydata <- data.frame("x.diamond" = c(rewp_robust$ci.lb,
                                          coef(rewp_robust),
                                          rewp_robust$ci.ub,
                                          coef(rewp_robust)),
                          "y.diamond" = c(nrow(rewp_sa_data) + 1,
                                          nrow(rewp_sa_data) + 1.3,
                                          nrow(rewp_sa_data) + 1,
                                          nrow(rewp_sa_data) + 0.7),
                          "diamond_group" = rep(1, times = 4))

ylabels <- c("Summary", rewp_robust$slab)
nrows <- nrow(rewp_sa_data)
ybreaklim <- nrows + 1
ybreaks <- c(sort(1:ybreaklim, decreasing = TRUE))

fig3c <- fig3c + geom_polygon(data = summarydata, aes(x = x.diamond, y = y.diamond, group = diamond_group), 
                              color= "black", fill = "#1e89c6", 
                              size = 0.1) + scale_y_continuous(name = "",
                                                               breaks = ybreaks, labels = ylabels)



fig3c

#CNV
cnv_si_data <- meta_data %>%
  filter(type == "CNV" & Description == "ideation")# 3 effect sizes, 1 study. review only 

cnv_sa_data <- meta_data %>%
  filter(type == "CNV" & Description == "attempt")# 4 effect sizes, 4 studies. Meta analyze

cnv_risk_data <- meta_data %>%
  filter(type == "CNV" & Description == "risk")# no effects

cnv_sasi_data <- meta_data %>% 
  filter(type == "CNV" & suicide_group == 'SA' & control_group == 'SI')# 1 effect sizes


## Analyses
### CNV SA Analysis
cnv_estimate <- rma.mv(yi = g ~ 1,
                        V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID),
                        data = cnv_sa_data,
                        method = "REML", slab = paste(paste(paste(Authors, Year, sep=", ("), ")", 
                                                            sep = ""), effect_id, sep = " "))

cnv_robust <- robust(cnv_estimate, cluster = cnv_sa_data$Study_ID,
                      adjust = TRUE)
summary(cnv_robust)# g = 0.44, p = 0.4157

W <- diag(1/cnv_sa_data$var_g)
X <- model.matrix(cnv_robust)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(cnv_robust$sigma2) / (sum(cnv_robust$sigma2) + (cnv_robust$k-cnv_robust$p)/sum(diag(P)))
#84.45%

#### Create Forest Plot
fig3d <- viz_forest(cnv_robust, type = "study_only", variant = "classic", study_labels = cnv_estimate[["slab"]], 
                    col = "#1e89c6", xlab = expression("Hedges' " * italic("g")), 
                    text_size = 5) + ggtitle("CNV")

summarydata <- data.frame("x.diamond" = c(cnv_robust$ci.lb,
                                          coef(cnv_robust),
                                          cnv_robust$ci.ub,
                                          coef(cnv_robust)),
                          "y.diamond" = c(nrow(cnv_sa_data) + 1,
                                          nrow(cnv_sa_data) + 1.3,
                                          nrow(cnv_sa_data) + 1,
                                          nrow(cnv_sa_data) + 0.7),
                          "diamond_group" = rep(1, times = 4))

ylabels <- c("Summary", cnv_robust$slab)
nrows <- nrow(cnv_sa_data)
ybreaklim <- nrows + 1
ybreaks <- c(sort(1:ybreaklim, decreasing = TRUE))

fig3d <- fig3d + geom_polygon(data = summarydata, aes(x = x.diamond, y = y.diamond, group = diamond_group), 
                              color= "black", fill = "#1e89c6", 
                              size = 0.1) + scale_y_continuous(name = "",
                                                               breaks = ybreaks, labels = ylabels)



fig3d
# N2
n2_si_data <- meta_data %>%
  filter(type == "N2" & Description == "ideation")# no effects

n2_sa_data <- meta_data %>%
  filter(type == "N2" & Description == "attempt")# 6 effect sizes, 1 study. narrative review

n2_risk_data <- meta_data %>%
  filter(type == "N2" & Description == "risk")# 2 effects, 1 study. Narrative Review

# P2
p2_si_data <- meta_data %>%
  filter(type == "P2" & Description == "ideation")# 3 effects, 1 study. narrative

p2_sa_data <- meta_data %>%
  filter(type == "P2" & Description == "attempt")# 1 effect, 1 study. Narrative Review

p2_risk_data <- meta_data %>%
  filter(type == "P2" & Description == "risk")# no studies



# Statistical Power -------------------------------------------------------
## Format all studies needed for mapower_ul
study1_data <- meta_data %>% 
  filter(Study_ID == "1") %>%
  mutate("study" = paste(paste(Authors, Year, sep=", ("), ")", 
                         sep = "")) %>% 
  select(study, "yi" = g, "lower" = lbci, "upper" = ubci)
 
study2_data <- meta_data %>% 
  filter(Study_ID == "2") %>%
  mutate("study" = paste(paste(Authors, Year, sep=", ("), ")", 
                         sep = "")) %>% 
  select(study, "yi" = g, "lower" = lbci, "upper" = ubci)

study3_data <- meta_data %>% 
  filter(Study_ID == "3") %>%
  mutate("study" = paste(paste(Authors, Year, sep=", ("), ")", 
                         sep = "")) %>% 
  select(study, "yi" = g, "lower" = lbci, "upper" = ubci)

study4_data <- meta_data %>% 
  filter(Study_ID == "4") %>%
  mutate("study" = paste(paste(Authors, Year, sep=", ("), ")", 
                         sep = "")) %>% 
  select(study, "yi" = g, "lower" = lbci, "upper" = ubci)

study5_data <- meta_data %>% 
  filter(Study_ID == "5") %>%
  mutate("study" = paste(paste(Authors, Year, sep=", ("), ")", 
                         sep = "")) %>% 
  select(study, "yi" = g, "lower" = lbci, "upper" = ubci)


study6_data <- meta_data %>% 
  filter(Study_ID == "6") %>%
  mutate("study" = paste(paste(Authors, Year, sep=", ("), ")", 
                         sep = "")) %>% 
  select(study, "yi" = g, "lower" = lbci, "upper" = ubci)

study7_data <- meta_data %>% 
  filter(Study_ID == "7") %>%
  mutate("study" = paste(paste(Authors, Year, sep=", ("), ")", 
                         sep = "")) %>% 
  select(study, "yi" = g, "lower" = lbci, "upper" = ubci)

study8_data <- meta_data %>% 
  filter(Study_ID == "8") %>%
  mutate("study" = paste(paste(Authors, Year, sep=", ("), ")", 
                         sep = "")) %>% 
  select(study, "yi" = g, "lower" = lbci, "upper" = ubci)

study9_data <- meta_data %>% 
  filter(Study_ID == "9") %>%
  mutate("study" = paste(paste(Authors, Year, sep=", ("), ")", 
                         sep = "")) %>% 
  select(study, "yi" = g, "lower" = lbci, "upper" = ubci)

study10_data <- meta_data %>% 
  filter(Study_ID == "10") %>%
  mutate("study" = paste(paste(Authors, Year, sep=", ("), ")", 
                         sep = "")) %>% 
  select(study, "yi" = g, "lower" = lbci, "upper" = ubci)

study11_data <- meta_data %>% 
  filter(Study_ID == "11") %>%
  mutate("study" = paste(paste(Authors, Year, sep=", ("), ")", 
                         sep = "")) %>% 
  select(study, "yi" = g, "lower" = lbci, "upper" = ubci)

study12_data <- meta_data %>% 
  filter(Study_ID == "12") %>%
  mutate("study" = paste(paste(Authors, Year, sep=", ("), ")", 
                         sep = "")) %>% 
  select(study, "yi" = g, "lower" = lbci, "upper" = ubci)

study13_data <- meta_data %>% 
  filter(Study_ID == "13") %>%
  mutate("study" = paste(paste(Authors, Year, sep=", ("), ")", 
                         sep = "")) %>% 
  select(study, "yi" = g, "lower" = lbci, "upper" = ubci)

study14_data <- meta_data %>% 
  filter(Study_ID == "14") %>%
  mutate("study" = paste(paste(Authors, Year, sep=", ("), ")", 
                         sep = "")) %>% 
  select(study, "yi" = g, "lower" = lbci, "upper" = ubci)

study15_data <- meta_data %>% 
  filter(Study_ID == "15") %>%
  mutate("study" = paste(paste(Authors, Year, sep=", ("), ")", 
                         sep = "")) %>% 
  select(study, "yi" = g, "lower" = lbci, "upper" = ubci)

study16_data <- meta_data %>% 
  filter(Study_ID == "16") %>%
  mutate("study" = paste(paste(Authors, Year, sep=", ("), ")", 
                         sep = "")) %>% 
  select(study, "yi" = g, "lower" = lbci, "upper" = ubci)

study17_data <- meta_data %>% 
  filter(Study_ID == "17") %>%
  mutate("study" = paste(paste(Authors, Year, sep=", ("), ")", 
                         sep = "")) %>% 
  select(study, "yi" = g, "lower" = lbci, "upper" = ubci)

study18_data <- meta_data %>% 
  filter(Study_ID == "18") %>%
  mutate("study" = paste(paste(Authors, Year, sep=", ("), ")", 
                         sep = "")) %>% 
  select(study, "yi" = g, "lower" = lbci, "upper" = ubci)

study19_data <- meta_data %>% 
  filter(Study_ID == "19") %>%
  mutate("study" = paste(paste(Authors, Year, sep=", ("), ")", 
                         sep = "")) %>% 
  select(study, "yi" = g, "lower" = lbci, "upper" = ubci)

study20_data <- meta_data %>% 
  filter(Study_ID == "20") %>%
  mutate("study" = paste(paste(Authors, Year, sep=", ("), ")", 
                         sep = "")) %>% 
  select(study, "yi" = g, "lower" = lbci, "upper" = ubci)

study21_data <- meta_data %>% 
  filter(Study_ID == "21") %>%
  mutate("study" = paste(paste(Authors, Year, sep=", ("), ")", 
                         sep = "")) %>% 
  select(study, "yi" = g, "lower" = lbci, "upper" = ubci)

study22_data <- meta_data %>% 
  filter(Study_ID == "22") %>%
  mutate("study" = paste(paste(Authors, Year, sep=", ("), ")", 
                         sep = "")) %>% 
  select(study, "yi" = g, "lower" = lbci, "upper" = ubci)

study23_data <- meta_data %>% 
  filter(Study_ID == "23") %>%
  mutate("study" = paste(paste(Authors, Year, sep=", ("), ")", 
                         sep = "")) %>% 
  select(study, "yi" = g, "lower" = lbci, "upper" = ubci)

study24_data <- meta_data %>% 
  filter(Study_ID == "24") %>%
  mutate("study" = paste(paste(Authors, Year, sep=", ("), ")", 
                         sep = "")) %>% 
  select(study, "yi" = g, "lower" = lbci, "upper" = ubci)

study25_data <- meta_data %>% 
  filter(Study_ID == "25") %>%
  mutate("study" = paste(paste(Authors, Year, sep=", ("), ")", 
                         sep = "")) %>% 
  select(study, "yi" = g, "lower" = lbci, "upper" = ubci)

study26_data <- meta_data %>% 
  filter(Study_ID == "26") %>%
  mutate("study" = paste(paste(Authors, Year, sep=", ("), ")", 
                         sep = "")) %>% 
  select(study, "yi" = g, "lower" = lbci, "upper" = ubci)

study27_data <- meta_data %>% 
  filter(Study_ID == "27") %>%
  mutate("study" = paste(paste(Authors, Year, sep=", ("), ")", 
                         sep = "")) %>% 
  select(study, "yi" = g, "lower" = lbci, "upper" = ubci)

## Calculate statistical power given range of effect sizes for each study effect. Put in bogus true effect
power_1 <- mapower_ul(dat = study1_data, 
                          observed_es = 0.00, 
                          name = study1_data$study[1])

power1_data <- power_1$power_median_dat

power_2 <- mapower_ul(dat = study2_data, 
                      observed_es = 0.00, 
                      name = study2_data$study[1])
power2_data <- power_2$power_median_dat

power_3 <- mapower_ul(dat = study3_data, 
                      observed_es = 0.00, 
                      name = study3_data$study[1])
power3_data <- power_3$power_median_dat

power_4 <- mapower_ul(dat = study4_data, 
                      observed_es = 0.00, 
                      name = study4_data$study[1])
power4_data <- power_4$power_median_dat

power_5 <- mapower_ul(dat = study5_data, 
                      observed_es = 0.00, 
                      name = study5_data$study[1])
power5_data <- power_5$power_median_dat

power_6 <- mapower_ul(dat = study6_data, 
                      observed_es = 0.00, 
                      name = study6_data$study[1])
power6_data <- power_6$power_median_dat

power_7 <- mapower_ul(dat = study7_data, 
                      observed_es = 0.00, 
                      name = study7_data$study[1])
power7_data <- power_7$power_median_dat

power_8 <- mapower_ul(dat = study8_data, 
                      observed_es = 0.00, 
                      name = study8_data$study[1])
power8_data <- power_8$power_median_dat

power_9 <- mapower_ul(dat = study9_data, 
                      observed_es = 0.00, 
                      name = study9_data$study[1])
power9_data <- power_9$power_median_dat

power_10 <- mapower_ul(dat = study10_data, 
                      observed_es = 0.00, 
                      name = study10_data$study[1])
power10_data <- power_10$power_median_dat

power_11 <- mapower_ul(dat = study11_data, 
                      observed_es = 0.00, 
                      name = study11_data$study[1])
power11_data <- power_11$power_median_dat

power_12 <- mapower_ul(dat = study12_data, 
                      observed_es = 0.00, 
                      name = study12_data$study[1])
power12_data <- power_12$power_median_dat

power_13 <- mapower_ul(dat = study13_data, 
                      observed_es = 0.00, 
                      name = study13_data$study[1])
power13_data <- power_13$power_median_dat

power_14 <- mapower_ul(dat = study14_data, 
                      observed_es = 0.00, 
                      name = study14_data$study[1])
power14_data <- power_14$power_median_dat

power_15 <- mapower_ul(dat = study15_data, 
                      observed_es = 0.00, 
                      name = study15_data$study[1])
power15_data <- power_15$power_median_dat

power_16 <- mapower_ul(dat = study16_data, 
                      observed_es = 0.00, 
                      name = study16_data$study[1])
power16_data <- power_16$power_median_dat

power_17 <- mapower_ul(dat = study17_data, 
                      observed_es = 0.00, 
                      name = study17_data$study[1])
power17_data <- power_17$power_median_dat

power_18 <- mapower_ul(dat = study18_data, 
                      observed_es = 0.00, 
                      name = study18_data$study[1])
power18_data <- power_18$power_median_dat

power_19 <- mapower_ul(dat = study19_data, 
                      observed_es = 0.00, 
                      name = study19_data$study[1])
power19_data <- power_19$power_median_dat

power_20 <- mapower_ul(dat = study20_data, 
                      observed_es = 0.00, 
                      name = study20_data$study[1])
power20_data <- power_20$power_median_dat

power_21 <- mapower_ul(dat = study21_data, 
                      observed_es = 0.00, 
                      name = study21_data$study[1])
power21_data <- power_21$power_median_dat

power_22 <- mapower_ul(dat = study22_data, 
                      observed_es = 0.00, 
                      name = study22_data$study[1])
power22_data <- power_22$power_median_dat

power_23 <- mapower_ul(dat = study23_data, 
                      observed_es = 0.00, 
                      name = study23_data$study[1])
power23_data <- power_23$power_median_dat

power_24 <- mapower_ul(dat = study24_data, 
                      observed_es = 0.00, 
                      name = study24_data$study[1])
power24_data <- power_24$power_median_dat

power_25 <- mapower_ul(dat = study25_data, 
                      observed_es = 0.00, 
                      name = study25_data$study[1])
power25_data <- power_25$power_median_dat

power_26 <- mapower_ul(dat = study26_data, 
                      observed_es = 0.00, 
                      name = study26_data$study[1])
power26_data <- power_26$power_median_dat

power_27 <- mapower_ul(dat = study27_data, 
                      observed_es = 0.00, 
                      name = study27_data$study[1])
power27_data <- power_27$power_median_dat

powerall_data <- rbind(power1_data, power2_data, power3_data, power4_data, 
                       power5_data, power6_data, power7_data, power8_data, power9_data,
                       power10_data, power11_data, power12_data, power13_data, 
                       power14_data, power15_data, power16_data, power17_data, 
                       power18_data, power19_data, power20_data, power21_data, 
                       power22_data, power23_data, power24_data, power25_data, 
                       power26_data, power27_data
                       )

powerall_data <- powerall_data %>% 
  select("study_name" = meta_analysis_name, es01:es1)

min_max_power <- (min(powerall_data[,11])) * 100

# Calculate maximum observed median effect size power
max_max_power <- (max(powerall_data[,11])) * 100

# Calculate mean observed median effect size power
mean_max_power <- (mean(powerall_data[,11])) * 100

# Calculate median observed median effect size power
median_max_power <- (median(powerall_data[,11])) * 100

long_power_data <- pivot_longer(powerall_data, es01:es1, names_to = "effect_size")


sufficient_power <- long_power_data %>% 
  filter(value >= 0.8)

mean(.4, .8, .8, .8, .4, .9, .4, 1.0, .9, .8, 1.0, .9, .9, .9, .9, 
       .7, .8, .8, .7, .8, .5) 
# Create Figures ----------------------------------------------------------

#Figure 2
figure2 <- plot_grid(fig2a, fig2b, fig2c, fig2d, labels = c('A', 'B', 'C', 'D'))

save_plot("fig2.png", figure2, ncol = 2, base_height = 12, base_width = 6, dpi = 600)

figure2

#Figure 3
figure3 <- plot_grid(fig3a, fig3b, fig3c, fig3d, labels = c('A', 'B', 'C', 'D'))

save_plot("fig3.png", figure3, ncol = 2, base_height = 12, base_width = 6, dpi = 600)

figure3


#Figure 4
figure4 <- plot_grid(fig4a, fig4b, fig4c, labels = c('A', 'B', 'C'))

save_plot("fig4.png", figure4, ncol = 2, base_height = 12, base_width = 6, dpi = 600)

figure4



#Figure 5
figure5 <- plot_grid(fig5a, fig5b, labels = c('A', 'B', 'C'))

save_plot("fig5.png", figure5, ncol = 2, base_height = 12, base_width = 6, dpi = 600)

figure5



#Figure 6
firepower_plot <- ggplot(data = long_power_data) +
  geom_tile(aes(
    x = effect_size,
    y = reorder(study_name, desc(study_name)),
    fill = value
  )) +
  theme(aspect.ratio = 0.3) +
  scale_fill_gradient2(
    name = "Power",
    midpoint = 0.5,
    low = "white",
    mid = "orange",
    high = "red"
  )+ theme_tufte(base_family = "Helvetica")

firepower_plot <-
  firepower_plot + theme(strip.text.x = element_blank())

firepower_plot <- firepower_plot + labs(x ="Effect size", y = "")

firepower_plot <-
  firepower_plot + scale_x_discrete(
    labels = c(
      "es01" = "0.1",
      "es02" = "0.2",
      "es03" = "0.3",
      "es04" = "0.4",
      "es05" = "0.5",
      "es06" = "0.6",
      "es07" = "0.7",
      "es08" = "0.8",
      "es09" = "0.9",
      "es1" = "1"
    ))
firepower_plot <- firepower_plot + theme(text = element_text(size = 20))

tiff(file = "figure_6.tiff", width = 12, height = 10, units = "in", 
     res = 800, compression = "lzw")
firepower_plot
dev.off()


# Supplemental Analyses ---------------------------------------------------
#Subset data for overall analyses
ideation_data <- meta_data %>% 
  filter(Description == "ideation" & abs_g < 2.9)#One outlier effect size 



ideation_data$type <- recode_factor(as.factor(ideation_data$type), LDAEP = 0, 
                                    LPP = 1, RewP = 2, P3 = 3, CNV = 4, P2 = 5)

attempt_data <- meta_data %>% 
  filter(Description == "attempt" & abs_g < 1.9)#three outliers 

attempt_data$type <- recode_factor(as.factor(attempt_data$type), LDAEP = 0, 
                                   LPP = 1, RewP = 2, P3 = 3, CNV = 4, N2 = 5, 
                                   N1 = 6, PINV = 6, SPN = 6)

risk_data <- meta_data %>% 
  filter(Description == "risk")

risk_data$type <- recode_factor(as.factor(risk_data$type), LDAEP = 0, LPP = 1, 
                                   RewP = 2, P3 = 3, N2 = 4)

sasi_data <- meta_data %>% 
  filter(suicide_group == 'SA' & control_group == 'SI' & abs_g < 1.90)# Two outliers

sasi_data$type <- recode_factor(as.factor(sasi_data$type), LDAEP = 0, LPP = 1, 
                 P3 = 2, N2 = 3, CNV = 4, N1 = 4)

#Ideation meta-analysis##############
ideation_estimate <- rma.mv(yi = abs_g ~ 1, 
                    V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID), 
                    data = ideation_data,
                    method = "REML", 
                    slab = paste(paste(paste(Authors, Year, sep=", ("), ")", 
                                       sep = ""), effect_id, sep = " "))

robust_ideation <- robust(ideation_estimate, cluster = ideation_data$Study_ID,
                          adjust = TRUE)

summary(robust_ideation)

#Compute I^2 for this model
W <- diag(1/ideation_data$var_g)
X <- model.matrix(robust_ideation)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(robust_ideation$sigma2) / (sum(robust_ideation$sigma2) + (robust_ideation$k-robust_ideation$p)/sum(diag(P)))
# I2 = 18.46%
#Funnel Plot ideation
### create contour enhanced funnel plot (with funnel centered at 0)
funnel_1 <- viz_sunset(robust_ideation, 
           contours = TRUE,
           power_contours = "continuous", 
           text_size = 5, power_stats = FALSE, 
           x_breaks = c(-1.5, -1.0, -0.5, 0, 0.5, 1, 1.5), 
           x_limit = c(-1.75, 1.75), 
           xlab = expression("Hedges' " * italic("g"))) + theme(legend.position = "none")

funnel_power <- metafor::funnel(robust_ideation)# Used later to calculate power

### SI publication bias
ideation_weights <- weights(robust_ideation)

ideation_data$weights <- ideation_weights

ideation_data <- mutate(ideation_data, "sqrt_weights" = sqrt(weights))

ideation_estimate <- rma.mv(yi = abs_g ~ sqrt_weights, 
                            V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID), 
                            data = ideation_data,
                            method = "REML")

ideation_eggers <- robust(ideation_estimate, cluster = ideation_data$Study_ID, 
                          adjust = TRUE)

summary(ideation_eggers)# significant small sample bias

average_ideation <- ideation_data %>% #aggregate effect sizes within studies
      group_by(Study_ID) %>% 
      summarize("mean_abs_g" = mean(abs_g), 
                "mean_var_g" = mean(var_g))

weightfunct(effect = average_ideation$mean_abs_g, 
            v = average_ideation$mean_var_g)

### SI statistical power
ideation_data$power <- s_power(se = funnel_power$y, 
                     true_effect = 0.24, 
                     sig_level = 0.05)

mean(ideation_data$power)# 0.209
min(ideation_data$power)# 0.094
max(ideation_data$power)# 0.398

#ERP type moderation ideation
ideation_estimate <- rma.mv(yi = abs_g ~ type, 
                                 V = var_g, 
                            random = list(~ 1 | effect_id, ~ 1 | Study_ID), 
                                 data = ideation_data,
                                 method = "REML")

ideation.type.mod <- robust(ideation_estimate, cluster = ideation_data$Study_ID, 
                            adjust = TRUE)

summary(ideation.type.mod)# Overall effect of type, no estimates significant

#ERP Scoring method moderation ideation
ideation_estimate <- rma.mv(yi = abs_g ~ scoring, 
                                       V = var_g, 
                            random = list(~ 1 | effect_id, ~ 1 | Study_ID), 
                                       data = ideation_data,
                                       method = "REML")

ideation.scoring.mod <- robust(ideation_estimate, cluster = ideation_data$Study_ID, 
                               adjust = TRUE)

summary(ideation.scoring.mod) #No effect of scoring method

#Age moderation ideation
ideation_estimate <- rma.mv(yi = abs_g ~ mean_age, 
                              V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID), 
                              data = ideation_data,
                              method = "REML")

ideation.age.mod <- robust(ideation_estimate, cluster = ideation_data$Study_ID, adjust = TRUE)

summary(ideation.age.mod)#No evidence that age matters for strength of effect

#Plot figure 4
fig.4 <- ggplot(data = ideation_data, aes(x = sqrt_weights, y = abs_g, color = as.factor(Study_ID))) + geom_point(size = 5) + theme_classic(base_size = 20)+ geom_abline(slope = ideation_eggers$b[2], intercept = ideation_eggers$b[1], color = "black", size = 1.0)#Supplying custom regressin line
tiff(file = "figure_4_eggers.tiff", width = 10, height = 12, units = "in", 
     res = 800, compression = "lzw")
fig.4 + labs(x = expression(sqrt("W")["i"]), y = expression("Hedges' " * italic("g")), color = "Study") + scale_color_stata(labels = c("Albanese et al. (2019b)", "Song et al. (2019)", "Tsypes, Owens, & Gibb (2019)", "Baik et al. (2018)", "Kudinova et al. (2016)", "Lee et al. (2014)", "Marsic (2012)", "Albanese & Schmidt (unpublished, 2020)", "Gallyer et al. (2020)")) + 
  scale_y_continuous(breaks = c(0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00)) 
dev.off()

#Attempt meta-analysis############
attempt_estimate <- rma.mv(yi = abs_g ~ 1, 
                            V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID), 
                            data = attempt_data,
                            method = "REML", 
                           slab = paste(paste(paste(Authors, Year, sep=", ("), ")", sep = ""), effect_id, sep = " "))

robust_attempt <- robust(attempt_estimate, cluster = attempt_data$Study_ID,
                          adjust = TRUE)

summary(robust_attempt)

W <- diag(1/attempt_data$var_g)
X <- model.matrix(robust_attempt)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(robust_attempt$sigma2) / (sum(robust_attempt$sigma2) + (robust_attempt$k-robust_attempt$p)/sum(diag(P)))

jpeg(file = "forest1.jpg", width = 10, height = 12, units = "in", 
     res = 800)

forest(robust_attempt, xlab = expression("Hedges' " * italic("g")), 
       mlab = "Robust RE Model", cex = 0.95)

par(cex=1.25, font=2)

### add column headings to the plot
text(-4.65, 48, "Author(s) (Year) Effect ID",  pos=4)
text(6.75, 48, expression(italic("g ") * "[95% CI]"), pos=2)

dev.off()

#Funnel Plot attempt
funnel_2 <- viz_sunset(robust_attempt, 
           contours = TRUE,
           power_contours =  "continuous", 
           text_size = 5, 
           power_stats = FALSE, x_breaks = c(-1.5, -1.0, -0.5, 0, 0.5, 1, 1.5), 
           x_limit = c(-1.75, 1.75), 
           xlab = expression("Hedges' " * italic("g"))) + theme(legend.position = "none")

funnel_power <- metafor::funnel(robust_attempt)

### SA publication bias
attempt_weights <- weights(robust_attempt)

attempt_data$weights <- attempt_weights

attempt_data <- mutate(attempt_data, "sqrt_weights" = sqrt(weights))

attempt_estimate <- rma.mv(yi = abs_g ~ sqrt_weights, 
                            V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID), 
                            data = attempt_data,
                            method = "REML")

attempt_eggers <- robust(attempt_estimate, cluster = attempt_data$Study_ID, adjust = TRUE)

summary(attempt_eggers)# significant small sample bias

average_attempt <- attempt_data %>% #aggregate effect sizes within studies
      group_by(Study_ID) %>% 
      summarize("mean_abs_g" = mean(abs_g), 
                "mean_var_g" = mean(var_g))

weightfunct(effect = average_attempt$mean_abs_g, v = average_attempt$mean_var_g)

### SA statistical power
attempt_data$power <- s_power(se = funnel_power$y, 
                               true_effect = 0.40, 
                               sig_level = 0.05)

mean(attempt_data$power)# 0.303
min(attempt_data$power)# 0.103
max(attempt_data$power)# 0.810


#ERP Type Moderation attempt
attempt_estimate <- rma.mv(yi = abs_g ~ type, 
       V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID), 
       data = attempt_data,
       method = "REML")

attempt.type.mod <- robust(attempt_estimate, cluster = attempt_data$Study_ID, adjust = TRUE)

summary(attempt.type.mod) #Significant effect with RewP and LPP being significantly smaller effects than LDAEP

#ERP Scoring Method Moderation attempt
attempt_estimate <- rma.mv(yi = abs_g ~ scoring, 
                                    V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID), 
                                    data = attempt_data,
                                    method = "REML")

attempt.scoring.mod <- robust(attempt_estimate, cluster = attempt_data$Study_ID, adjust = TRUE)

summary(attempt.scoring.mod)# Not significant, scoring does not seem to account for heterogeneity of effect sizes
#Sample age moderation attempt
attempt_estimate <- rma.mv(yi = abs_g ~ mean_age, 
                            V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID), 
                            data = attempt_data,
                            method = "REML", verbose = TRUE, control=list(optimizer="optim", optmethod="Nelder-Mead")) #Didn't converge first time, switched optimizer and opt method

attempt.age.mod <- robust(attempt_estimate, cluster = attempt_data$Study_ID, adjust = TRUE)

summary(attempt.age.mod)# No evidece that age matters

#Risk meta-analysis############
risk_estimate <- rma.mv(yi = abs_g ~ 1, 
                           V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID), 
                           data = risk_data,
                           method = "REML", 
                        slab = paste(paste(paste(Authors, Year, sep=", ("), ")", sep = ""), effect_id, sep = " "))

robust_risk <- robust(risk_estimate, cluster = risk_data$Study_ID,
                         adjust = TRUE)

summary(robust_risk)

W <- diag(1/risk_data$var_g)
X <- model.matrix(robust_risk)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(robust_risk$sigma2) / (sum(robust_risk$sigma2) + (robust_risk$k-robust_risk$p)/sum(diag(P)))


jpeg(file = "forest2.jpg", width = 10, height = 12, units = "in", 
     res = 800)

forest(robust_risk, xlab = expression("Hedges' " * italic("g")), 
       mlab = "Robust RE Model", cex = 0.95)

par(cex=1.25, font=2)

### add column headings to the plot
text(-3.80, 46.5, "Author(s) (Year) Effect ID",  pos=4)
text(5.20, 46.5, expression(italic("g ") * "[95% CI]"), pos=2)

dev.off()
#Funnel plot risk

funnel_3 <- viz_sunset(robust_risk, 
           contours = TRUE,
           power_contours =  "continuous", 
           text_size = 5, 
           power_stats = FALSE, x_breaks = c(-1.5, -1.0, -0.5, 0, 0.5, 1, 1.5), 
           x_limit = c(-1.75, 1.75), 
           xlab = expression("Hedges' " * italic("g"))) + theme(legend.position = "none")


funnel_power <- metafor::funnel(robust_risk)


### SR publication bias
risk_weights <- weights(robust_risk)

risk_data$weights <- risk_weights

risk_data <- mutate(risk_data, "sqrt_weights" = sqrt(weights))

risk_estimate <- rma.mv(yi = abs_g ~ sqrt_weights, 
                           V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID), 
                           data = risk_data,
                           method = "REML")

risk_eggers <- robust(risk_estimate, cluster = risk_data$Study_ID, adjust = TRUE)

summary(risk_eggers)# significant small sample bias

average_risk <- risk_data %>% #aggregate effect sizes within studies
      group_by(Study_ID) %>% 
      summarize("mean_abs_g" = mean(abs_g), 
                "mean_var_g" = mean(var_g))

weightfunct(effect = average_risk$mean_abs_g, v = average_risk$mean_var_g)

### SR statistical power
risk_data$power <- s_power(se = funnel_power$y, 
                              true_effect = 0.33, 
                              sig_level = 0.05)

mean(risk_data$power)# 0.286
min(risk_data$power)# 0.119
max(risk_data$power)# 0.649

#ERP type moderation risk
risk_estimate <- rma.mv(yi = abs_g ~ type, 
                        V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID), 
                        data = risk_data,
                        method = "REML")

risk.type.mod <- robust(risk_estimate, cluster = risk_data$Study_ID,
                      adjust = TRUE)

summary(risk.type.mod)#Not significant

#ERP scoring method moderation risk
risk_estimate <- rma.mv(yi = abs_g ~ scoring, 
                        V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID), 
                        data = risk_data,
                        method = "REML")

risk.scoring.mod <- robust(risk_estimate, cluster = risk_data$Study_ID,
                        adjust = TRUE)

summary(risk.scoring.mod) #No significant differences
#Sample age moderation risk
risk_estimate <- rma.mv(yi = abs_g ~ mean_age, 
                            V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID), 
                            data = risk_data,
                            method = "REML", control=list(iter.max=10000, rel.tol=1e-8))

risk.age.mod <- robust(risk_estimate, cluster = risk_data$Study_ID, adjust = TRUE)

summary(risk.age.mod)#No evidence that age matters for strength of effect, though it is close

#SA vs. SI meta-analysis############
sasi_estimate <- rma.mv(yi = abs_g ~ 1, 
                        V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID), 
                        data = sasi_data,
                        method = "REML", 
                        slab = paste(paste(paste(Authors, Year, sep=", ("), ")", sep = ""), effect_id, sep = " "))

robust_sasi <- robust(sasi_estimate, cluster = sasi_data$Study_ID,
                      adjust = TRUE)

summary(robust_sasi)
# Calculate I^2
W <- diag(1/sasi_data$var_g)
X <- model.matrix(robust_sasi)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(robust_sasi$sigma2) / (sum(robust_sasi$sigma2) + (robust_sasi$k-robust_sasi$p)/sum(diag(P)))

# SA vs. SI forest plot
jpeg(file = "forest3.jpg", width = 10, height = 12, units = "in", 
     res = 800)
forest(robust_sasi, xlab = expression("Hedges' " * italic("g")), 
       mlab = "Robust RE Model", cex = 0.95)
par(cex=1.25, font=2)

### add column headings to the plot
text(-4.50, 33, "Author(s) (Year) Effect ID",  pos=4)
text(6.45, 33, expression(italic("g ") * "[95% CI]"), pos=2)

dev.off()
# SA vs. SI funnel plot
funnel_4 <- viz_sunset(robust_sasi, 
           contours = TRUE,
           power_contours =  "continuous", 
           text_size = 5, 
           power_stats = FALSE, x_breaks = c(-1.5, -1.0, -0.5, 0, 0.5, 1, 1.5), 
           x_limit = c(-1.75, 1.75), 
           xlab = expression("Hedges' " * italic("g"))) + theme(legend.position = "none")

funnel_power <- metafor::funnel(robust_sasi)


### SR publication bias
sasi_weights <- weights(robust_sasi)

sasi_data$weights <- sasi_weights

sasi_data <- mutate(sasi_data, "sqrt_weights" = sqrt(weights))

sasi_estimate <- rma.mv(yi = abs_g ~ sqrt_weights, 
                        V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID), 
                        data = sasi_data,
                        method = "REML")

sasi_eggers <- robust(sasi_estimate, cluster = sasi_data$Study_ID, adjust = TRUE)

summary(sasi_eggers)# significant small sample bias

average_sasi <- sasi_data %>% #aggregate effect sizes within studies
      group_by(Study_ID) %>% 
      summarize("mean_abs_g" = mean(abs_g), 
                "mean_var_g" = mean(var_g))

weightfunct(effect = average_sasi$mean_abs_g, v = average_sasi$mean_var_g)

### SR statistical power
sasi_data$power <- s_power(se = funnel_power$y, 
                           true_effect = 0.33, 
                           sig_level = 0.05)

mean(sasi_data$power)# 0.231
min(sasi_data$power)# 0.086
max(sasi_data$power)# 0.648
#ERP type moderation risk
sasi_estimate <- rma.mv(yi = abs_g ~ type,
                        V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID), 
                        data = sasi_data,
                        method = "REML", control=list(maxiter=100000), verbose = TRUE)

sasi.type.mod <- robust(sasi_estimate, cluster = sasi_data$Study_ID,
                        adjust = TRUE)

summary(sasi.type.mod)#Not significant

#ERP scoring method moderation risk
sasi_estimate <- rma.mv(yi = abs_g ~ scoring, 
                        V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID), 
                        data = sasi_data,
                        method = "REML")

sasi.scoring.mod <- robust(sasi_estimate, cluster = sasi_data$Study_ID,
                           adjust = TRUE)

summary(sasi.scoring.mod) #No significant differences
#Sample age moderation risk
sasi_estimate <- rma.mv(yi = abs_g ~ mean_age, 
                        V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID), 
                        data = sasi_data,
                        method = "REML", control=list(iter.max=10000, rel.tol=1e-8))

sasi.age.mod <- robust(sasi_estimate, cluster = sasi_data$Study_ID, adjust = TRUE)

summary(sasi.age.mod)#Evidence that age matters

#Create Figure 8 showing age effect
figure8_data <- sasi_data %>% 
  filter(!is.na(mean_age))
fig.8 <- ggplot(data = figure8_data, aes(x = mean_age, y = abs_g, color = as.factor(Study_ID))) + geom_point(size = 5) + theme_classic(base_size = 20)+ geom_abline(slope = sasi.age.mod$b[2], intercept = sasi.age.mod$b[1], color = "black", size = 1.0)#Supplying custom regressin line

tiff(file = "figure_8_agemod.tiff", width = 10, height = 12, units = "in", 
     res = 800, compression = "lzw")

fig.8 + labs(x = "Average Age", y = expression("Hedges' " * italic("g")), color = "Study") + scale_color_stata(labels = c("Albanese et al. (2019a)","Weinberg et al. (2017)", "Kim & Park (2013)", "Min et al. (2012)", "Chen et al. (2005)", "Hansenne et al. (1996)")) + 
  scale_x_continuous(breaks = c(15, 20, 25, 30, 35, 40, 45, 50))

dev.off()

# Sensitivity Analyses############
#Ideation sensitivity analyses
#Analysis if we didn't remove outliers
ideation_data <- meta_data %>% 
      filter(Description == "ideation")

ideation_estimate <- rma.mv(yi = abs_g ~ 1, 
                            V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID), 
                            data = ideation_data,
                            method = "REML", slab = paste(paste(paste(Authors, Year, sep=", ("), ")", sep = ""), effect_id, sep = " "))

robust_ideation <- robust(ideation_estimate, cluster = ideation_data$Study_ID,
                          adjust = TRUE)
summary(robust_ideation)
#Analysis with only published studies
ideation_data <- meta_data %>% 
      filter(Description == "ideation" & Published == 1 & abs_g < 2.9)

ideation_estimate <- rma.mv(yi = abs_g ~ 1, 
                            V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID), 
                            data = ideation_data,
                            method = "REML", slab = paste(paste(paste(Authors, Year, sep=", ("), ")", sep = ""), effect_id, sep = " "))

robust_ideation <- robust(ideation_estimate, cluster = ideation_data$Study_ID,
                          adjust = TRUE)
summary(robust_ideation)

#Attempt sensitivity analyses
#Analysis if we didn't remove outliers
attempt_data <- meta_data %>% 
      filter(Description == "attempt")

attempt_estimate <- rma.mv(yi = abs_g ~ 1, 
                           V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID), 
                           data = attempt_data,
                           method = "REML", 
                           slab = paste(paste(paste(Authors, Year, sep=", ("), ")", sep = ""), effect_id, sep = " "))

robust_attempt <- robust(attempt_estimate, cluster = attempt_data$Study_ID,
                         adjust = TRUE)
summary(robust_attempt)
#Analysis using only published studies
attempt_data <- meta_data %>% 
      filter(Description == "attempt" & Published == 1 & abs_g < 1.9)

attempt_estimate <- rma.mv(yi = abs_g ~ 1, 
                           V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID), 
                           data = attempt_data,
                           method = "REML", 
                           slab = paste(paste(paste(Authors, Year, sep=", ("), ")", sep = ""), effect_id, sep = " "))

robust_attempt <- robust(attempt_estimate, cluster = attempt_data$Study_ID,
                         adjust = TRUE)
summary(robust_attempt)

#Suicide risk sensitivity analysis
risk_data <- meta_data %>% 
      filter(Description == "risk" & Published == 1)

risk_estimate <- rma.mv(yi = abs_g ~ 1, 
                        V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID), 
                        data = risk_data,
                        method = "REML", 
                        slab = paste(paste(paste(Authors, Year, sep=", ("), ")", sep = ""), effect_id, sep = " "))

robust_risk <- robust(risk_estimate, cluster = risk_data$Study_ID,
                      adjust = TRUE)
summary(robust_risk)

#SA vs. SI-only sensitivity analysis
sasi_data <- meta_data %>% 
      filter(suicide_group == 'SA' & control_group == 'SI')

sasi_estimate <- rma.mv(yi = abs_g ~ 1, 
                        V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID), 
                        data = sasi_data,
                        method = "REML", 
                        slab = paste(paste(paste(Authors, Year, sep=", ("), ")", sep = ""), effect_id, sep = " "))

robust_sasi <- robust(sasi_estimate, cluster = sasi_data$Study_ID,
                      adjust = TRUE)

summary(robust_sasi)

sasi_data <- meta_data %>% 
      filter(suicide_group == 'SA' & control_group == 'SI' & Published == 1 & abs_g < 1.90)

sasi_estimate <- rma.mv(yi = abs_g ~ 1, 
                        V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID), 
                        data = sasi_data,
                        method = "REML", 
                        slab = paste(paste(paste(Authors, Year, sep=", ("), ")", sep = ""), effect_id, sep = " "))

robust_sasi <- robust(sasi_estimate, cluster = sasi_data$Study_ID,
                      adjust = TRUE)

summary(robust_sasi)


#Check mean sample size across ERP type in attempt data############
attempt_data <- meta_data %>% 
      filter(Description == "attempt" & abs_g < 1.9)

attempt_erp <- attempt_data %>% 
      group_by(type) %>% 

### Create Figures
## Figure 3
figure3 <- (funnel_1 + funnel_2)/(funnel_3 + funnel_4)
tiff(file = "figure3_funnel.tiff", width = 10, height = 12, units = "in", 
     res = 800, compression = "lzw")

figure2 + plot_annotation(tag_levels = "A")

dev.off()
# Figure 2
year_data <- read_excel("meta_spreadsheet.xlsx", sheet = "Full Text Review")

year_data <- year_data %>% 
      filter(Study_ID > 0)

year_data$Year[year_data$Year == "unpublished"] <- "2020"

summary_year <- year_data %>% 
      group_by(Year) %>% 
      summarize(n = n())

empty_years <- data.frame(seq(1994, 2020, 1), rep(0, length(seq(1994, 2020, 1))))

empty_years <- rename(empty_years, "Year" = seq.1994..2020..1., "n" = rep.0..length.seq.1994..2020..1...)

summary_year$Year <- as.numeric(summary_year$Year)

missing_years <- anti_join(empty_years, summary_year, by = "Year")

summary_year$n <- as.numeric(summary_year$n)

all_years <- union(summary_year, missing_years, by = "Year") %>% 
      arrange(Year)

jpeg(file = "studies_time.jpg", width = 10, height = 12, units = "in", 
     res = 800)
ggplot(data = all_years, aes(x = Year, y = n, group = 1)) + geom_line(color = "red", size = 2) + geom_point(size = 3) + 
      theme_classic(base_size = 20) + labs(y = "# of Studies Included") + scale_x_continuous(limits = c(1994, 2020), breaks = seq(1994,2020,1)) + scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, 2)) + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
 

