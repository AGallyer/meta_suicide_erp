library(tidyverse)
library(metafor)
library(readxl)
library(compute.es)
library(PRISMAstatement)
library(metaviz)
library(ggthemes)
library(patchwork)
# Function
s_power <- function(se, true_effect, sig_level) {
      
      (1 - stats::pnorm(stats::qnorm(1 - sig_level/2) * 
                              se, abs(true_effect), se)) + 
            stats::pnorm(stats::qnorm(sig_level/2) * 
                               se, abs(true_effect), se)
}
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
  mutate("abs_g" = abs(g))

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

#Compute I^2 for this model
W <- diag(1/ldaep_si_data$var_g)
X <- model.matrix(ldaep_robust)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(ldaep_robust$sigma2) / (sum(ldaep_robust$sigma2) + (ldaep_robust$k-ldaep_robust$p)/sum(diag(P)))

#Create forest plot
fig2a <- viz_forest(ldaep_robust, variant = "classic", study_labels = ldaep_estimate[["slab"]], 
                    col = "#1e89c6", xlab = expression("Hedges' " * italic("g")), 
           text_size = 5)

### LDAEP SA Analysis
ldaep_estimate <- rma.mv(yi = g ~ 1,
                         V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID),
                         data = ldaep_sa_data,
                         method = "REML")

ldaep_robust <- robust(ldaep_estimate, cluster = ldaep_sa_data$Study_ID,
                       adjust = TRUE)
summary(ldaep_robust)# g = 0.021, p = 0.946

W <- diag(1/ldaep_sa_data$var_g)
X <- model.matrix(ldaep_robust)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(ldaep_robust$sigma2) / (sum(ldaep_robust$sigma2) + (ldaep_robust$k-ldaep_robust$p)/sum(diag(P)))

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

summary(ldaep_eggers)# significant small sample bias


### LDAEP Risk Analysis
ldaep_estimate <- rma.mv(yi = g ~ 1,
                         V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID),
                         data = ldaep_risk_data,
                         method = "REML")

ldaep_robust <- robust(ldaep_estimate, cluster = ldaep_risk_data$Study_ID,
                       adjust = TRUE)
summary(ldaep_robust)# g = 0.175, p = 0.552

W <- diag(1/ldaep_risk_data$var_g)
X <- model.matrix(ldaep_robust)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(ldaep_robust$sigma2) / (sum(ldaep_robust$sigma2) + (ldaep_robust$k-ldaep_robust$p)/sum(diag(P)))

# LPP
lpp_si_data <- meta_data %>%
  filter(type == "LPP" & Description == "ideation")# 34 effect sizes, 3 studies. Meta analyze

lpp_sa_data <- meta_data %>%
  filter(type == "LPP" & Description == "attempt")# 1 study, 3 effect sizes. Meta analyze

lpp_risk_data <- meta_data %>%
  filter(type == "LPP" & Description == "risk")# 9 effect sizes, 2 studies. Meta analyze

## Analyses
###LPP SI Analysis
lpp_estimate <- rma.mv(yi = g ~ 1,
                       V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID),
                       data = lpp_si_data,
                       method = "REML")

lpp_robust <- robust(lpp_estimate, cluster = lpp_si_data$Study_ID,
                     adjust = TRUE)
summary(lpp_robust)# g = -0.053, p = 0.407

### LPP Risk Analysis
lpp_estimate <- rma.mv(yi = g ~ 1,
                       V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID),
                       data = lpp_risk_data,
                       method = "REML")

lpp_robust <- robust(lpp_estimate, cluster = lpp_risk_data$Study_ID,
                     adjust = TRUE)
summary(lpp_robust)# g = -0.290, p = 0.157
# P3
p3_si_data <- meta_data %>%
  filter(type == "P3" & Description == "ideation")# 53 effect sizes, 3 studies. meta analyze

p3_sa_data <- meta_data %>%
  filter(type == "P3" & Description == "attempt")# 15 effect sizes, 6 studies. Meta analyze

p3_risk_data <- meta_data %>%
  filter(type == "P3" & Description == "risk")# 28 effect sizes, 3 studies. meta analyze

## Analyses
###P3 SI Analysis
p3_estimate <- rma.mv(yi = g ~ 1,
                      V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID),
                      data = p3_si_data,
                      method = "REML")

p3_robust <- robust(p3_estimate, cluster = p3_si_data$Study_ID,
                    adjust = TRUE)
summary(p3_robust)# g = 0.067, p = 0.827

###P3 SA Analysis
p3_estimate <- rma.mv(yi = g ~ 1,
                      V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID),
                      data = p3_sa_data,
                      method = "REML")

p3_robust <- robust(p3_estimate, cluster = p3_sa_data$Study_ID,
                    adjust = TRUE)
summary(p3_robust)# g = -0.118, p = 0.799

### P3 Risk Analysis
p3_estimate <- rma.mv(yi = g ~ 1,
                      V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID),
                      data = p3_risk_data,
                      method = "REML")

p3_robust <- robust(p3_estimate, cluster = p3_risk_data$Study_ID,
                    adjust = TRUE)

summary(p3_robust)# g = 0.250, p = 0.322

# RewP
rewp_si_data <- meta_data %>%
  filter(type == "RewP" & Description == "ideation")# 12 effect sizes, 3 studies. Meta analyze

rewp_sa_data <- meta_data %>%
  filter(type == "RewP" & Description == "attempt")# 5 effect sizes, 2 studies. Meta analyze

rewp_risk_data <- meta_data %>%
  filter(type == "RewP" & Description == "risk")# 3 effect sizes, 1 study. Meta analyze

## Analyses
### RewP SI Analysis
rewp_estimate <- rma.mv(yi = g ~ 1,
                        V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID),
                        data = rewp_si_data,
                        method = "REML")

rewp_robust <- robust(rewp_estimate, cluster = rewp_si_data$Study_ID,
                      adjust = TRUE)

summary(rewp_robust)# g = 0.007, p = 0.885

### RewP SA Analysis
rewp_estimate <- rma.mv(yi = g ~ 1,
                        V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID),
                        data = rewp_sa_data,
                        method = "REML")

rewp_robust <- robust(rewp_estimate, cluster = rewp_sa_data$Study_ID,
                      adjust = TRUE)

summary(rewp_robust)# g = -0.122, p = 0.545

#CNV
cnv_si_data <- meta_data %>%
  filter(type == "CNV" & Description == "ideation")# 3 effect sizes, 1 study. Meta analyze

cnv_sa_data <- meta_data %>%
  filter(type == "CNV" & Description == "attempt")# 4 effect sizes, 4 studies. Meta analyze

cnv_risk_data <- meta_data %>%
  filter(type == "CNV" & Description == "risk")# no effects

## Analyses
### CNV SA Analysis
cnv_estimate <- rma.mv(yi = g ~ 1,
                       V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID),
                       data = cnv_sa_data,
                       method = "REML")

cnv_robust <- robust(cnv_estimate, cluster = cnv_sa_data$Study_ID,
                     adjust = TRUE)

summary(cnv_robust)# g = 0.413, p = 0.416
# N2
n2_si_data <- meta_data %>%
  filter(type == "N2" & Description == "ideation")# no effects

n2_sa_data <- meta_data %>%
  filter(type == "N2" & Description == "attempt")# 6 effect sizes, 1 study. Meta analyze

n2_risk_data <- meta_data %>%
  filter(type == "N2" & Description == "risk")# 2 effects, 1 study. Narrative Review

# P2
p2_si_data <- meta_data %>%
  filter(type == "P2" & Description == "ideation")# 3 effects, 1 study. Meta analyze

p2_sa_data <- meta_data %>%
  filter(type == "P2" & Description == "attempt")# 1 effect, 1 study. Narrative Review

p2_risk_data <- meta_data %>%
  filter(type == "P2" & Description == "risk")# no studies

### Supplemental Analyses
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
fig.4 + labs(x = expression(sqrt("W")["i"]), y = expression("Hedges' " * italic("g")), color = "Study") + scale_color_stata(labels = c("Albanese et al. (2019b)", "Song et al. (2019)", "Tsypes, Owens, & Gibb (2019)", "Baik et al. (2018)", "Kudinova et al. (2016)", "Lee et al. (2014)", "Marsic (2012)", "Albanese (unpublished, 2020)", "Gallyer et al. (2020)")) + 
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
 

