library(tidyverse)
library(metafor)
library(readxl)
library(compute.es)
#Import data
data <- read_excel("meta_spreadsheet.xlsx", sheet = "effect_size_coding")
meta_data <- data %>% 
  select(Study_ID, effect_id, Description, type, g, var_g, everything()) %>% 
  mutate("abs_g" = abs(g))
meta_data$scoring <- recode_factor(as.factor(meta_data$scoring), raw = 0, 
                                   difference = 1, slope = 1)
meta_data$type <- recode_factor(as.factor(meta_data$type), LDAEP = 0, LPP = 1, 
                                RewP = 2, P3 = 3, 
                                N1 = 4, CNV = 4, N2 = 4, P2 = 4, PINV = 4, SPN = 4)
#Subset data
ideation_data <- meta_data %>% 
  filter(Description == "ideation" & abs_g < 2.9)#One outlier effect size

attempt_data <- meta_data %>% 
  filter(Description == "attempt" & abs_g < 1.9)#three outliers

risk_data <- meta_data %>% 
  filter(Description == "risk")

#Ideation meta-analysis
ideation_estimate <- rma.mv(yi = abs_g ~ 1, 
                    V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID), 
                    data = ideation_data,
                    method = "REML", slab = paste(paste(paste(Authors, Year, sep=", ("), ")", sep = ""), effect_id, sep = " "))

robust_ideation <- robust(ideation_estimate, cluster = ideation_data$Study_ID,
                          adjust = TRUE)
#Forest plot ideation
tiff(file = "ideation_forest.tiff", width = 10, height = 12, units = "in", 
     res = 800, compression = "lzw")
forest(robust_ideation, xlab = expression("Hedges' " * italic("g")), 
       mlab = "Robust RE Model")
par(cex=1.25, font=2)

### add column headings to the plot
text(-3.5, 54, "Author(s) (Year) Effect ID",  pos=4)
text(4.3, 54, expression(italic("g ") * "[95% CI]"), pos=2)
dev.off()

#Funnel Plot ideation
### create contour enhanced funnel plot (with funnel centered at 0)
tiff(file = "ideation_funnel.tiff", width = 10, height = 12, units = "in", 
     res = 800, compression = "lzw")
par(mar=c(5,4,1,2))
funnel(robust_ideation, level=c(90, 95, 99), shade=c("white", "gray55", "gray75"), refline=0, legend=TRUE,
       xlab = expression("Hedges' " * italic("g")))
dev.off()
#ERP type moderation ideation
ideation_estimate <- rma.mv(yi = abs_g ~ type, 
                                 V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID), 
                                 data = ideation_data,
                                 method = "REML")

ideation.type.mod <- robust(ideation_estimate, cluster = ideation_data$Study_ID, adjust = TRUE)

summary(ideation.type.mod)# Overall effect of type, no estimates significantly stronger or weaker than LDAEP

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

summary(ideation.N.mod) #Effect diminishes from .40 to .1 when N increases. Suggests to me that there is an inflating of effect sizes for the small studies
#Attempt meta-analysis
attempt_estimate <- rma.mv(yi = abs_g ~ 1, 
                            V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID), 
                            data = attempt_data,
                            method = "REML", 
                           slab = paste(paste(paste(Authors, Year, sep=", ("), ")", sep = ""), effect_id, sep = " "))

robust_attempt <- robust(attempt_estimate, cluster = attempt_data$Study_ID,
                          adjust = TRUE)

tiff(file = "attempt_forest.tiff", width = 10, height = 12, units = "in", 
     res = 800, compression = "lzw")

forest(robust_attempt, xlab = expression("Hedges' " * italic("g")), 
       mlab = "Robust RE Model")
par(cex=1.25, font=2)

### add column headings to the plot
text(-4.75, 48, "Author(s) (Year) Effect ID",  pos=4)
text(7.0, 48, expression(italic("g ") * "[95% CI]"), pos=2)
dev.off()

#Funnel Plot attempt
tiff(file = "attempt_funnel.tiff", width = 10, height = 12, units = "in", 
     res = 800, compression = "lzw")
par(mar=c(5,4,1,2))
funnel(robust_attempt, level=c(90, 95, 99), shade=c("white", "gray55", "gray75"), refline=0, legend=TRUE,
       xlab = expression("Hedges' " * italic("g")))
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
#Risk meta-analysis
risk_estimate <- rma.mv(yi = abs_g ~ 1, 
                           V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID), 
                           data = risk_data,
                           method = "REML", 
                        slab = paste(paste(paste(Authors, Year, sep=", ("), ")", sep = ""), effect_id, sep = " "))

robust_risk <- robust(risk_estimate, cluster = risk_data$Study_ID,
                         adjust = TRUE)


tiff(file = "risk_forest.tiff", width = 10, height = 12, units = "in", 
     res = 800, compression = "lzw")

forest(robust_risk, xlab = expression("Hedges' " * italic("g")), 
       mlab = "Robust RE Model")
par(cex=1.25, font=2)

### add column headings to the plot
text(-3.80, 13.5, "Author(s) (Year) Effect ID",  pos=4)
text(5.0, 13.5, expression(italic("g ") * "[95% CI]"), pos=2)
dev.off()

#Funnel plot risk
tiff(file = "risk_funnel.tiff", width = 10, height = 12, units = "in", 
     res = 800, compression = "lzw")
par(mar=c(5,3,1,1))
funnel(robust_risk, level=c(90, 95, 99), shade=c("white", "gray55", "gray75"), refline=0, legend=TRUE,
       xlab = expression("Hedges' " * italic("g")))
dev.off()

#ERP type moderation risk
risk_estimate <- rma.mv(yi = abs_g ~ type, 
                        V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID), 
                        data = risk_data,
                        method = "REML")

risk.type.mod <- robust(risk_estimate, cluster = risk_data$Study_ID,
                      adjust = TRUE)

summary(risk.type.mod)#Not significant, likely due to power but unsure.
#ERP scoring method moderation risk
risk_estimate <- rma.mv(yi = abs_g ~ scoring, 
                        V = var_g, random = list(~ 1 | effect_id, ~ 1 | Study_ID), 
                        data = risk_data,
                        method = "REML")

risk.scoring.mod <- robust(risk_estimate, cluster = risk_data$Study_ID,
                        adjust = TRUE)

summary(risk.scoring.mod) #Again no differences
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