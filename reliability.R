
# packages ----------------------------------------------------------------
library(tidyverse)
library(readxl)
library(irr)
library(psych)

# import data -------------------------------------------------------------
austin_data <- read_excel("meta_spreadsheet_austin.xlsx", sheet = "effect_size_coding")

austin_data <- select(austin_data, effect_id, design, Description, mean_age, n_suicide, n_control, scoring, type, N)

sean_data <- read_excel("meta_spreadsheet_sean.xlsx", sheet = "effect_size_coding")

sean_data <- select(sean_data, effect_id, sean_design, sean_Description, sean_mean_age, 
                    sean_n_suicide, sean_n_control, sean_scoring, sean_type, sean_N)

full_data <- full_join(austin_data, sean_data, by = "effect_id")

design_data <- select(full_data, design, sean_design)

description_data <- select(full_data, Description, sean_Description)

age_data <- select(full_data, mean_age, sean_mean_age)

n_suicide_data <- select(full_data, n_suicide, sean_n_suicide)

n_control_data <- select(full_data, n_control, sean_n_control)

scoring_data <- select(full_data, scoring, sean_scoring)

type_data <- select(full_data, type, sean_type)

N_data <- select(full_data, N, sean_N)


# Agreement rate for categorical ----------------------------------------------------------

agree(design_data) #99.5%

agree(description_data)# 97.9%

agree(scoring_data)# 100%

agree(type_data)# 87.5%


# Kappa for categorical variables -----------------------------------------

kappa2(design_data)# 0.989

kappa2(description_data)# 0.965

kappa2(scoring_data)# 1.0

kappa2(type_data)# 0.833


# Intercoder correlation for continuous -----------------------------------
cor(age_data, use = "pairwise.complete.obs")# .99

cor(n_suicide_data, use = "pairwise.complete.obs")# 1

cor(n_control_data, use = "pairwise.complete.obs")# 1

cor(N_data, use = "pairwise.complete.obs")# .99

# ICC for continuous data -------------------------------------------------
ICC(age_data)# 1.00

ICC(n_suicide_data)# 1.00

ICC(n_control_data)# 1.00

ICC(N_data)# 0.99


# Label discrepancies for discussion --------------------------------------
different <- (design_data$design == design_data$sean_design)

design_data$different <- different

different <- (description_data$Description == description_data$sean_Description)

description_data$different <- different

different <- (scoring_data$scoring == scoring_data$sean_scoring)

scoring_data$different <- different

different <- (type_data$type == type_data$sean_type)

type_data$different <- different

different <- (N_data$N == N_data$sean_N)
N_data$different <- different