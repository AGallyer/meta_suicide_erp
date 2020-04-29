library(tidyverse)
library(metafor)
library(readxl)
library(compute.es)
library(PRISMAstatement)
library(multcomp)

#PRISMA statement
prisma(found = 270,
       found_other = 5,
       no_dupes = 200, 
       screened = 200, 
       screen_exclusions = 149, 
       full_text = 51,
       full_text_exclusions = 24, 
       qualitative = 27, 
       quantitative = 24,
       width = 800, height = 800)

#Import data
data <- read_excel("meta_spreadsheet.xlsx", sheet = "effect_size_coding")
meta_data <- data %>% 
  dplyr::select(Study_ID, effect_id, Description, type, g, var_g, everything()) %>% 
  mutate("abs_g" = abs(g))
meta_data$scoring <- recode_factor(as.factor(meta_data$scoring), raw = 0, 
                                   difference = 1, slope = 1)
meta_data$type <- recode_factor(as.factor(meta_data$type), LDAEP = 0, LPP = 1, 
                                RewP = 2, P3 = 3, PSW = 4,
                                N1 = 5, CNV = 5, N2 = 5, P2 = 5, PINV = 5, SPN = 5)
#Subset data
ideation_data <- meta_data %>% 
  filter(Description == "ideation" & abs_g < 2.9)#One outlier effect size 

attempt_data <- meta_data %>% 
  filter(Description == "attempt" & abs_g < 1.9)#three outliers 

risk_data <- meta_data %>% 
  filter(Description == "risk")

#Ideation meta-analysis##############
ideation_estimate <- rma.mv(yi = abs_g ~ 1, 
                    V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID), 
                    data = ideation_data,
                    method = "REML", slab = paste(paste(paste(Authors, Year, sep=", ("), ")", sep = ""), effect_id, sep = " "))

robust_ideation <- robust(ideation_estimate, cluster = ideation_data$Study_ID,
                          adjust = TRUE)

#Compute I^2 for this model
W <- diag(1/ideation_data$var_g)
X <- model.matrix(robust_ideation)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(robust_ideation$sigma2) / (sum(robust_ideation$sigma2) + (robust_ideation$k-robust_ideation$p)/sum(diag(P)))

#Funnel Plot ideation
### create contour enhanced funnel plot (with funnel centered at 0)
tiff(file = "ideation_funnel.tiff", width = 10, height = 12, units = "in", 
     res = 800, compression = "lzw")
par(mar=c(5,4,1,2))
funnel(robust_ideation, xlab = expression("Hedges' " * italic("g")))
dev.off()
#ERP type moderation ideation
ideation_estimate <- rma.mv(yi = abs_g ~ type, 
                                 V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID), 
                                 data = ideation_data,
                                 method = "REML")

ideation.type.mod <- robust(ideation_estimate, cluster = ideation_data$Study_ID, adjust = TRUE)

summary(ideation.type.mod)# Overall effect of type, no estimates significantly stronger or weaker than LDAEP

# Exploratory moderation to test which ERP has the strongest average effect size
# ideation_data$type <- fct_expand(ideation_data$type, "6")
# 
# ideation_data$type <- fct_recode(ideation_data$type, "6" = "0")
# 
# 
# ideation_estimate <- rma.mv(yi = abs_g ~ type - 1,
#                             V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID),
#                             data = ideation_data,
#                             method = "REML")
# 
# ideation.type.mod2 <- robust(ideation_estimate, cluster = ideation_data$Study_ID, adjust = TRUE)
# 
# 
# summary(glht(ideation.type.mod2, linfct=rbind(c(1, -1, 0, 0, 0, 0), c(1, 0, -1, 0, 0, 0),
#                                               c(1, 0, 0, -1, 0, 0), c(1, 0, 0, 0, -1, 0),
#                                               c(1, 0, 0, 0, 0, -1), c(0, 1, -1, 0, 0, 0),
#                                               c(0, 1, 0, -1, 0, 0), c(0, 1, 0, 0, -1, 0),
#                                               c(0, 1, 0, 0, 0, -1), c(0, 0, 1, -1, 0, 0),
#                                               c(0, 0, 1, 0, -1, 0), c(0, 0, 1, 0, 0, -1),
#                                               c(0, 0, 0, 1, -1, 0), c(0, 0, 0, 1, 0, -1),
#                                               c(0, 0, 0, 0, 1, -1))), test = adjusted("holm"))
#ERP Scoring method moderation ideation
ideation_estimate <- rma.mv(yi = abs_g ~ scoring, 
                                       V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID), 
                                       data = ideation_data,
                                       method = "REML")

ideation.scoring.mod <- robust(ideation_estimate, cluster = ideation_data$Study_ID, adjust = TRUE)

summary(ideation.scoring.mod) #No effect of scoring method, surprising result to me
#Age moderation ideation
ideation_estimate <- rma.mv(yi = abs_g ~ mean_age, 
                              V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID), 
                              data = ideation_data,
                              method = "REML")

ideation.age.mod <- robust(ideation_estimate, cluster = ideation_data$Study_ID, adjust = TRUE)

summary(ideation.age.mod)#No evidence that age matters for strength of effect
#Sample size moderation ideation
ideation_data <- ideation_data %>% 
  mutate("cent_n" = N - min(N))
ideation_estimate <- rma.mv(yi = abs_g ~ cent_n, 
                              V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID), 
                              data = ideation_data,
                              method = "REML")

ideation.N.mod <- robust(ideation_estimate, cluster = ideation_data$Study_ID, adjust = TRUE)

summary(ideation.N.mod) #Effect diminishes when N increases. Suggests to me that there is an inflating of effect sizes for the small studies

#Plot figure 3
fig.3 <- ggplot(data = ideation_data, aes(x = N, y = abs_g, color = as.factor(Study_ID))) + geom_point(size = 3) + theme_classic()+ geom_abline(slope = ideation.N.mod$b[2], intercept = ideation.N.mod$b[1], color = "black", size = 1.0)#Supplying custom regressin line
tiff(file = "figure_3_N.tiff", width = 10, height = 12, units = "in", 
     res = 800, compression = "lzw")
fig.3 + labs(x = "Effect Sample Size (N)", y = expression("Hedges' " * italic("g")), color = "Study") + scale_color_brewer(palette = "Dark2", labels = c("Albanese et al. (2019b)", "Song et al. (2019)", "Tsypes, Owens, & Gibb (2019)", "Baik et al. (2018)", "Kudinova et al. (2016)", "Lee et al. (2014)", "Marsic (2012)", "Albanese (unpublished, 2020)")) + 
  scale_y_continuous(breaks = c(0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00)) + scale_x_continuous(breaks = c(50, 75, 100, 125, 150, 175, 200, 225, 250))
dev.off()

#Attempt meta-analysis############
attempt_estimate <- rma.mv(yi = abs_g ~ 1, 
                            V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID), 
                            data = attempt_data,
                            method = "REML", 
                           slab = paste(paste(paste(Authors, Year, sep=", ("), ")", sep = ""), effect_id, sep = " "))

robust_attempt <- robust(attempt_estimate, cluster = attempt_data$Study_ID,
                          adjust = TRUE)

W <- diag(1/attempt_data$var_g)
X <- model.matrix(robust_attempt)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(robust_attempt$sigma2) / (sum(robust_attempt$sigma2) + (robust_attempt$k-robust_attempt$p)/sum(diag(P)))


tiff(file = "attempt_forest.tiff", width = 10, height = 12, units = "in", 
     res = 800, compression = "lzw")

forest(robust_attempt, xlab = expression("Hedges' " * italic("g")), 
       mlab = "Robust RE Model")
par(cex=1.25, font=2)

### add column headings to the plot
text(-4.65, 48, "Author(s) (Year) Effect ID",  pos=4)
text(6.75, 48, expression(italic("g ") * "[95% CI]"), pos=2)
dev.off()

#Funnel Plot attempt
tiff(file = "attempt_funnel.tiff", width = 10, height = 12, units = "in", 
     res = 800, compression = "lzw")
par(mar=c(5,4,1,2))
funnel(robust_attempt, xlab = expression("Hedges' " * italic("g")))
dev.off()
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
                            method = "REML")

attempt.age.mod <- robust(attempt_estimate, cluster = attempt_data$Study_ID, adjust = TRUE)

summary(attempt.age.mod)#Evidence that age matters, with the effect getting SMALLER as age increases
#Sample size moderation attempt
attempt_data <- attempt_data %>% 
  mutate("cent_n" = N - min(N))

attempt_estimate <- rma.mv(yi = abs_g ~ cent_n, 
                              V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID), 
                              data = attempt_data,
                              method = "REML")

attempt.N.mod <- robust(attempt_estimate, cluster = attempt_data$Study_ID, adjust = TRUE)

summary(attempt.N.mod) #Not significant, but could be due to sample size (i.e., 12 studies for 46 effects)

attempt_data <- attempt_data %>% 
  filter(!is.na(mean_age))
#Plot significant moderation of effect by sample age
fig.6 <- ggplot(data = attempt_data, aes(x = mean_age, y = abs_g, color = as.factor(Study_ID))) + geom_point(size = 3) + theme_classic()+ geom_abline(slope = attempt.age.mod$b[2], intercept = attempt.age.mod$b[1], color = "black", size = 1.0)#Supplying custom regressin line
tiff(file = "figure_6_age.tiff", width = 10, height = 12, units = "in", 
     res = 800, compression = "lzw")
fig.6 + labs(x = "Average Age", y = expression("Hedges' " * italic("g")), color = "Study") + scale_color_brewer(palette = "Dark2", labels = c("Albanese et al. (2019a)","Weinberg et al. (2017)", "Kim & Park (2013)", "Min et al. (2012)", "Jandl, Steyer, & Kaschka (2010)", "Chen et al. (2005)", "Hansenne et al. (1996)", "Ashton et al. (1994)")) + 
  scale_x_continuous(breaks = c(15, 20, 25, 30, 35, 40, 45, 50))
dev.off()

#Risk meta-analysis############
risk_estimate <- rma.mv(yi = abs_g ~ 1, 
                           V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID), 
                           data = risk_data,
                           method = "REML", 
                        slab = paste(paste(paste(Authors, Year, sep=", ("), ")", sep = ""), effect_id, sep = " "))

robust_risk <- robust(risk_estimate, cluster = risk_data$Study_ID,
                         adjust = TRUE)

W <- diag(1/risk_data$var_g)
X <- model.matrix(robust_risk)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(robust_risk$sigma2) / (sum(robust_risk$sigma2) + (robust_risk$k-robust_risk$p)/sum(diag(P)))


tiff(file = "risk_forest.tiff", width = 10, height = 12, units = "in", 
     res = 800, compression = "lzw")

forest(robust_risk, xlab = expression("Hedges' " * italic("g")), 
       mlab = "Robust RE Model")
par(cex=1.25, font=2)

### add column headings to the plot
text(-3.80, 39.5, "Author(s) (Year) Effect ID",  pos=4)
text(5.35, 39.5, expression(italic("g ") * "[95% CI]"), pos=2)
dev.off()

#Funnel plot risk
tiff(file = "risk_funnel.tiff", width = 10, height = 12, units = "in", 
     res = 800, compression = "lzw")
par(mar=c(5,4,1,2))
funnel(robust_risk, xlab = expression("Hedges' " * italic("g")))
dev.off()

#ERP type moderation risk
risk_estimate <- rma.mv(yi = abs_g ~ type, 
                        V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID), 
                        data = risk_data,
                        method = "REML")

risk.type.mod <- robust(risk_estimate, cluster = risk_data$Study_ID,
                      adjust = TRUE)

summary(risk.type.mod)#Significant, but no individual estimates different form the LDAEP

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
#Sample size moderation risk
risk_data <- risk_data %>% 
  mutate("cent_n" = N - min(N))

risk_estimate <- rma.mv(yi = abs_g ~ cent_n, 
                        V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID), 
                        data = risk_data,
                        method = "REML")

risk.N.mod <- robust(risk_estimate, cluster = risk_data$Study_ID,
                           adjust = TRUE)

summary(risk.N.mod) #Did not find evidence of N having impact on effect size

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

#Check mean sample size across ERP type in attempt data############
attempt_data <- meta_data %>% 
      filter(Description == "attempt" & abs_g < 1.9)

attempt_erp <- attempt_data %>% 
      group_by(type) %>% 
      summarise("mean_n" = mean(N))

meta_sample <- meta_data %>% 
      group_by(Study_ID) %>% 
      summarise("mean_n" = mean(N))

mean(meta_sample$mean_n)
sd(meta_sample$mean_n)

meta_sample <- meta_data %>% 
      group_by(Description, Study_ID) %>% 
      summarise("mean_n" = mean(N))

fig.9 <- ggplot(meta_sample, aes(x = fct_reorder(Description, mean_n, fun = median, .desc =TRUE), y = mean_n, fill = Description)) + 
      geom_boxplot() + theme_classic() + geom_jitter(position=position_jitter(0.2)) + 
      labs(x = "Meta-analysis STB Group", y = "Study N") + theme(legend.position="none", text = element_text(size = 12)) + 
      scale_fill_brewer(palette = "Dark2") + scale_x_discrete(labels = c('Suicidal Ideation','Suicide Attempt','Suicidality'))
tiff(file = "figure_9.tiff", width = 10, height = 12, units = "in", 
     res = 800, compression = "lzw")
fig.9
dev.off()