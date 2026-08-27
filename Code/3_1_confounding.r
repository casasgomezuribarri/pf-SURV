# This script checks checks for systematic differences
# in gametocytaemia and age at exposure across replicates

# author: ivan casas

# env
##########################################################################################################
library(readxl)
library(lme4)
library(dplyr)
library(tidyr)
library(lmerTest)
library(car)
source("Code/functions.r")


path <- "Data/experiments.xlsx"
data <- read_excel(path, sheet = "feeds")
# view(data)

data$temp_mean <- as.factor(data$temp_mean)

# differences in gametocytaemia and exflagellation?
##########################################################################################################

game_model <- glm(gametocytaemia_inf ~ temp_mean, data = data)
summary(game_model) # use for coef

# p computation type iii analysis - this is a massive overkill here but for consistency
# reencode vars into sum-to-0 contrasts  (needed for type iii analysis)
local_contrasts <- list(
    temp_mean = contr.sum(2)
)

# create model matrix manually with newly encoded vars
game_matrix <- model.matrix(~temp_mean, data = data, contrasts.arg = local_contrasts)[, -1]
data_sumto0 <- cbind(data, as.data.frame(game_matrix)) # bind to data

# fit full model
game_model_sumto0 <- glm(gametocytaemia_inf ~ temp_mean, data = data_sumto0)

# sanity check: these two models will have different results, and sumto0 rresults are not interpretable. Check the models are the same via other metrics
logLik(game_model) # 31.60124 (df=3)
logLik(game_model_sumto0) # 31.60124 (df=3)
AIC(game_model) # -57.20248
AIC(game_model_sumto0) # -57.20248

# fit nested model
no_tm <- glm(gametocytaemia_inf ~ 1, data = data_sumto0)
p_table_game <- rbind(
    calc_lrt(game_model_sumto0, no_tm, "temp_mean")
)

print(p_table_game)


# differences in age at feeding?
##########################################################################################################

# first pull all data into one longer column
datalong <- data %>%
    pivot_longer(
        cols = c(gambiae_age, coluzzii_age),
        names_to = "species",
        values_to = "age"
    ) |>
    mutate(species = factor(species))
# view(datalong)

age_model <- glm(age ~ temp_mean * species, data = datalong)
summary(age_model) # use oinly for coefs

# p computation type iii analysis
# reencode vars into sum-to-0 contrasts  (needed for type iii analysis)
local_contrasts <- list(
    temp_mean = contr.sum(2),
    species = contr.sum(2)
)

# create model matrix manually with newly encoded vars
age_matrix <- model.matrix(~ temp_mean * species, data = datalong, contrasts.arg = local_contrasts)[, -1]
colnames(age_matrix) <- gsub(":", "_", colnames(age_matrix)) # remove colons
datalong_sumto0 <- cbind(datalong, as.data.frame(age_matrix)) # bind to datalong

# fit full model
age_model_sumto0 <- glm(age ~ temp_mean1 + species1 + temp_mean1_species1, data = datalong_sumto0)

# sanity check: these two models will have different results, and sumto0 rresults are not interpretable. Check the models are the same via other metrics
logLik(age_model) # -27.46431 (df=5)
logLik(age_model_sumto0) # -27.46431 (df=5)
AIC(age_model) # 64.92861
AIC(age_model_sumto0) # 64.92861

# fit nested models
no_tm <- glm(age ~ species1 + species1:temp_mean1, data = datalong_sumto0)
no_sp <- glm(age ~ temp_mean1 + species1:temp_mean1, data = datalong_sumto0)
no_tm_sp <- glm(age ~ species1 + temp_mean1, data = datalong_sumto0)

p_table_age <- rbind(
    calc_lrt(age_model_sumto0, no_tm, "temp_mean"),
    calc_lrt(age_model_sumto0, no_sp, "species"),
    calc_lrt(age_model_sumto0, no_tm_sp, "temp_mean:species")
)

print(p_table_age)
