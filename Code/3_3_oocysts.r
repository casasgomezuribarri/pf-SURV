# This script replicates the oocyst analysis and generates the paper's figures
# oocyst prevalence over time
# oocyst prevalences + stats
# oocyst counts + stats

# author: ivan casas gomez-uribarri


##################################################################################################################################
# Environment
##################################################################################################################################
source("Code/compile_data_survival.r") # sets the right wd and prepares the data
source("Code/functions.r")

# load packages
packages <- c(
    "ggplot2",
    "tidyr",
    "emmeans",
    "cowplot",
    "dplyr",
    "tidyverse",
    "betareg",
    "GGally",
    "ggh4x",
    "lme4",
    "grid",
    "glmmTMB",
    "car",
    "legendry"
)
for (i in packages) {
    if (!require(i, character.only = TRUE)) install.packages(i)
    library(i, character.only = TRUE)
}

##################################################################################################################################
# Oocyst prevalence by treatment by day
##################################################################################################################################
# choose colours
selected_colors <- c("#0c36f6", "#5a8cd1", "#da0000", "#e976d6")

# subset data by species
prevalence_gambiae <- prevalence[prevalence$species == "An. gambiae", ]
prevalence_coluzzii <- prevalence[prevalence$species == "An. coluzzii", ]

# expand dataset with all possible combinations so all bars have the same width
prevalence_gambiae <- prevalence_gambiae %>%
    complete(dpi, treatment, fill = list(mean_prevalence = 0, se_prevalence = 0))
prevalence_coluzzii <- prevalence_coluzzii %>%
    complete(dpi, treatment, fill = list(mean_prevalence = 0, se_prevalence = 0))


# plot time series (we ended up choosing not to have smooth lines)
prevalence_g <- ggplot(prevalence_gambiae, aes(x = dpi, y = mean_prevalence, fill = treatment, color = treatment)) +
    ylim(0, 1) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.5, color = NA) +
    geom_errorbar(aes(ymin = mean_prevalence - se_prevalence, ymax = mean_prevalence + se_prevalence, color = treatment),
        position = "dodge", alpha = 0.4
    ) +
    # geom_smooth(
    #     aes(group = treatment),
    #     method = "loess",
    #     se = FALSE,
    #     size = 1.5
    # ) +
    scale_color_manual(values = selected_colors) +
    guides(color = "none") + # comment this if lines are shown
    scale_fill_manual(values = selected_colors) +
    scale_x_continuous(breaks = seq(min(prevalence_gambiae$dpi), max(prevalence_gambiae$dpi), by = 2)) +
    theme_minimal() +
    labs(
        title = "An. gambiae",
        x = "Days Post Infection (DPI)",
        y = "Prevalence (% positive)",
        fill = "Temperature",
        # color = "Temperature",
    ) +
    theme(
        panel.grid.minor = element_line(color = "gray"),
        axis.text = element_text(size = 26),
        axis.title = element_text(size = 30),
        plot.margin = margin(10, 10, 10, 10), # margins (t, r, b, l)
        plot.title = element_text(size = 35, hjust = 1),
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 24),
        legend.key.size = unit(1.5, "cm")
    )

prevalence_c <- ggplot(prevalence_coluzzii, aes(x = dpi, y = mean_prevalence, fill = treatment, color = treatment)) +
    ylim(0, 1) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.5, color = NA) +
    geom_errorbar(aes(ymin = mean_prevalence - se_prevalence + 0.001, ymax = mean_prevalence + se_prevalence, color = treatment),
        position = "dodge", alpha = 0.4
    ) +
    # geom_smooth(
    #     aes(group = treatment),
    #     method = "loess",
    #     se = FALSE,
    #     size = 1.5
    # ) +
    scale_color_manual(values = selected_colors) +
    guides(color = "none") + # comment this if lines are shown
    scale_fill_manual(values = selected_colors) +
    scale_x_continuous(breaks = seq(min(prevalence_coluzzii$dpi), max(prevalence_coluzzii$dpi), by = 2)) +
    theme_minimal() +
    labs(
        title = "An. coluzzii",
        x = "Days Post Infection (DPI)",
        y = "Prevalence (% positive)",
        fill = "Temperature",
        color = "Temperature",
    ) +
    theme(
        panel.grid.minor = element_line(color = "gray"),
        axis.text = element_text(size = 26),
        axis.title = element_text(size = 30),
        plot.margin = margin(10, 10, 10, 10), # margins (t, r, b, l)
        plot.title = element_text(size = 35, hjust = 1),
        legend.position = "none",
    )

legend <- get_legend(prevalence_g) # get legend from subplot
prevalence_g <- prevalence_g + theme(legend.position = "none") # remove it from there

# combine
subplots <- plot_grid(prevalence_c, prevalence_g, ncol = 1)
plots <- plot_grid(subplots, legend, ncol = 2, rel_widths = c(1, 0.2)) # add legend
# plots <- subplots
plots
ggsave(plot = plots, filename = "/Users/ivancasas/GitHub/Thesis/Chapters/03_SURV/pics/oocyst_prevalences_timeline.png", width = 15, height = 9)
ggsave(plot = plots, filename = "Figures/oocyst_prevalences_timeline.png", width = 15, height = 9)

##################################################################################################################################
# Oocyst prevalence by treatment
##################################################################################################################################

# summary dataset:
prev_summary <- replicates |> # replicates was created in the data compiler
    group_by(pot, species, mean_temp, temp_range) |>
    summarise(
        n_reps = n(),
        n_total = sum(n_midguts, na.rm = TRUE),
        mean_prev = weighted.mean(mean_prevalence, w = n_midguts, na.rm = TRUE),
        se_prev = { # weighted SD of replicate means
            w <- n_midguts / sum(n_midguts)
            sqrt(sum(w * (mean_prevalence - mean_prev)^2) / (n() - 1))
        },
        mean_cnt = weighted.mean(mean_count, w = n_midguts_positive, na.rm = TRUE),
        se_cnt = {
            w <- n_midguts_positive / sum(n_midguts_positive)
            sqrt(sum(w * (mean_count - mean_cnt)^2) / (n() - 1))
        },
        .groups = "drop"
    ) |>
    mutate(treatment = paste0(substr(pot, 3, 4), "±", substr(pot, 5, 5)))

# plot
ooPrev <- ggplot(prev_summary, aes(x = interaction(temp_range, mean_temp, species, sep = "!"), y = mean_prev, color = treatment)) +
    geom_point(size = 4) +
    geom_errorbar(
        aes(ymin = pmax(mean_prev - se_prev, 0), ymax = mean_prev + se_prev),
        width = 0.6
    ) +
    scale_color_manual(values = selected_colors) +
    scale_x_discrete( # fancy x axis
        guide = guide_axis_nested(key = "!"),
        name = NULL # hide default axis title
    ) +
    scale_y_continuous(limits = c(0, 0.8), breaks = seq(0, 0.8, 0.1)) +
    theme_minimal(base_size = 15) +
    labs(
        # title = "Oocyst prevalences across pots",
        x = "", # see scale_x_discrete
        y = "Prevalence (% positive)",
        color = "Temperature"
    ) +
    coord_cartesian(clip = "off") + # to allow marginal annotations for fancy axis
    annotation_custom(
        grob = textGrob("Range °C", gp = gpar(fontsize = 24), hjust = 1),
        xmin = 0.3, xmax = 0.3, ymin = -0.06, ymax = -0.06
    ) + # add marginal trange label
    annotation_custom(
        grob = textGrob("Mean °C", gp = gpar(fontsize = 24), hjust = 1),
        xmin = 0.3, xmax = 0.3, ymin = -0.1025, ymax = -0.1025
    ) + # marginal tmean label
    annotation_custom(
        grob = textGrob("Species", gp = gpar(fontsize = 24), hjust = 1),
        xmin = 0.3, xmax = 0.3, ymin = -0.145, ymax = -0.145
    ) + # species label
    theme(
        # axis.text.x = element_text(hjust = 1),
        plot.margin = margin(10, 10, 10, 60), # margins (t, r, b, l)
        panel.grid.minor = element_line(color = "gray"),
        axis.text = element_text(size = 26),
        axis.title = element_text(size = 30),
        plot.title = element_text(size = 35, hjust = 0.5),
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 24),
        legend.key.size = unit(1.5, "cm"),
        legend.position = "none"
    )

ooPrev
ggsave(filename = "/Users/ivancasas/GitHub/Thesis/Chapters/03_SURV/pics/oocyst_prevalences.png", plot = ooPrev, width = 12, height = 9)
ggsave(filename = "Figures/oocyst_prevalences.png", plot = ooPrev, width = 12, height = 9)

# stats
##################################################################################################################################
# prepare data
# view(dissection) # raw dissection data (1 row per msoquito)
ooprev_data <- filter(dissection, str_sub(pot, 2, 2) == "i" & midgut_dissection == TRUE) # only exposed + dissected mosquitoes
ooprev_data <- ooprev_data %>%
    group_by(replicate, pot) %>% # first add t_o and t_s for each replicte:pot
    mutate(
        t_o = dpi[which(oocysts != 0)][1], # day of the first oocyst seen
        t_s = dpi[which(sporozoites == TRUE)][1] # day of the first sporozoite seen
    ) %>%
    filter(dpi >= t_o & dpi < coalesce(t_s, Inf)) %>% # drop rows where dpi < t_o | dpi >= t_s  (if t_s is NA, then no upper bound to dpi)
    mutate(
        oo_positive = as.factor(ifelse(oocysts > 0, 1, 0)), # binary variable for oocyst presence (will fit a binomial model)
        species = case_when(
            str_sub(pot, 1, 1) == "k" ~ "An. gambiae",
            str_sub(pot, 1, 1) == "c" ~ "An. coluzzii"
        ),
        temp_mean = as.numeric(str_sub(pot, 3, 4)),
        temp_range = as.numeric(str_sub(pot, 5, 5)),
        treatment = case_when(
            str_sub(pot, 3, 5) == "210" ~ "21 ± 0",
            str_sub(pot, 3, 5) == "216" ~ "21 ± 6",
            str_sub(pot, 3, 5) == "270" ~ "27 ± 0",
            str_sub(pot, 3, 5) == "276" ~ "27 ± 6"
        )
    )

# and format them well
ooprev_data$pot <- as.factor(ooprev_data$pot)
ooprev_data$species <- as.factor(ooprev_data$species)
ooprev_data$temp_mean <- as.factor(ooprev_data$temp_mean)
ooprev_data$temp_range <- as.factor(ooprev_data$temp_range)

# gamabiae control 27±0 as reference for interpretability
ooprev_data$temp_mean <- relevel(ooprev_data$temp_mean, ref = "27")
ooprev_data$temp_range <- relevel(ooprev_data$temp_range, ref = "0")
ooprev_data$species <- relevel(ooprev_data$species, ref = "An. gambiae")
ooprev_data$pot <- relevel(ooprev_data$pot, ref = "ki270")


# binomial model for oocyst presence as a function of experimental group
prev_model <- glmer(oo_positive ~ species * temp_mean * temp_range + (1 | replicate) + (1 | pot_unique), data = ooprev_data, family = binomial)
summary(prev_model) # only for reporting coefs and SD of random effects - NOT pvalues


# p computation - type iii analysis framework. we do it manually for consistency across models:
###############################################################################################

# reencode vars into sum-to-0 contrasts  (needed for type iii analysis)
local_contrasts <- list(
    species = contr.sum(2),
    temp_mean = contr.sum(2),
    temp_range = contr.sum(2)
)

# create model matrix manually with newly encoded vars
prev_matrix <- model.matrix(~ species * temp_mean * temp_range, data = ooprev_data, contrasts.arg = local_contrasts)[, -1]
colnames(prev_matrix) <- gsub(":", "_", colnames(prev_matrix)) # remove colons from colnames
prev_data_sumto0 <- cbind(ooprev_data, as.data.frame(prev_matrix)) # bind to oo_prev data

# fit full model
prev_model_sumto0 <- glmer(oo_positive ~ species * temp_mean * temp_range + (1 | replicate) + (1 | pot_unique), data = prev_data_sumto0, family = binomial)
# this model is equivalent to the previous one, the only difference is the contrasts

# sanity checks
logLik(prev_model) # -290.8983 (df=10)
logLik(prev_model_sumto0) # -290.8983 (df=10)
AIC(prev_model) # 601.7966
AIC(prev_model_sumto0) # 601.7966

# fit nested models for each term
no_sp <- glmer(oo_positive ~ temp_mean1 + temp_range1 + species1_temp_mean1 + species1_temp_range1 + temp_mean1_temp_range1 + species1_temp_mean1_temp_range1 + (1 | replicate) + (1 | pot_unique), data = prev_data_sumto0, family = binomial)
no_tm <- glmer(oo_positive ~ species1 + temp_range1 + species1_temp_mean1 + species1_temp_range1 + temp_mean1_temp_range1 + species1_temp_mean1_temp_range1 + (1 | replicate) + (1 | pot_unique), data = prev_data_sumto0, family = binomial)
no_tr <- glmer(oo_positive ~ species1 + temp_mean1 + species1_temp_mean1 + species1_temp_range1 + temp_mean1_temp_range1 + species1_temp_mean1_temp_range1 + (1 | replicate) + (1 | pot_unique), data = prev_data_sumto0, family = binomial)
no_sp_tm <- glmer(oo_positive ~ species1 + temp_mean1 + temp_range1 + species1_temp_range1 + temp_mean1_temp_range1 + species1_temp_mean1_temp_range1 + (1 | replicate) + (1 | pot_unique), data = prev_data_sumto0, family = binomial)
no_sp_tr <- glmer(oo_positive ~ species1 + temp_mean1 + temp_range1 + species1_temp_mean1 + temp_mean1_temp_range1 + species1_temp_mean1_temp_range1 + (1 | replicate) + (1 | pot_unique), data = prev_data_sumto0, family = binomial)
no_tm_tr <- glmer(oo_positive ~ species1 + temp_mean1 + temp_range1 + species1_temp_mean1 + species1_temp_range1 + species1_temp_mean1_temp_range1 + (1 | replicate) + (1 | pot_unique), data = prev_data_sumto0, family = binomial)
no_sp_tm_tr <- glmer(oo_positive ~ species1 + temp_mean1 + temp_range1 + species1_temp_mean1 + species1_temp_range1 + temp_mean1_temp_range1 + (1 | replicate) + (1 | pot_unique), data = prev_data_sumto0, family = binomial)

p_table_prev <- rbind(
    calc_lrt(prev_model_sumto0, no_sp, "species"),
    calc_lrt(prev_model_sumto0, no_tm, "temp_mean"),
    calc_lrt(prev_model_sumto0, no_tr, "temp_range"),
    calc_lrt(prev_model_sumto0, no_sp_tm, "species:temp_mean"),
    calc_lrt(prev_model_sumto0, no_sp_tr, "species:temp_range"),
    calc_lrt(prev_model_sumto0, no_tm_tr, "temp_mean:temp_range"),
    calc_lrt(prev_model_sumto0, no_sp_tm_tr, "species:temp_mean:temp_range")
)

print(p_table_prev)

##################################################################################################################################
# Oocyst count by treatment
##################################################################################################################################

ooCount <- ggplot(prev_summary, aes(x = interaction(temp_range, mean_temp, species, sep = "!"), y = mean_cnt, color = treatment)) +
    geom_point(size = 4) +
    geom_errorbar(
        aes(ymin = mean_cnt - se_cnt, ymax = mean_cnt + se_cnt),
        width = 0.6
    ) +
    scale_color_manual(values = selected_colors) +
    scale_y_continuous(limits = c(0, NA)) +
    scale_x_discrete(
        guide = guide_axis_nested(key = "!"),
        name = NULL # hide default axis title
    ) +
    theme_minimal(base_size = 15) +
    labs(
        # title = "Oocyst counts across pots",
        x = "",
        y = "Count",
        color = "Temperature"
    ) +
    coord_cartesian(clip = "off") + # to allow marginal annotations
    annotation_custom(
        grob = textGrob("Range °C",
            x = unit(-0.02, "npc"), y = unit(-0.02, "npc"),
            gp = gpar(fontsize = 24), hjust = 1
        )
    ) + # add marginal trange label
    annotation_custom(
        grob = textGrob("Mean °C",
            x = unit(-0.02, "npc"), y = unit(-0.07, "npc"),
            gp = gpar(fontsize = 24), hjust = 1
        )
    ) + # marginal tmean label
    annotation_custom(
        grob = textGrob("Species",
            x = unit(-0.02, "npc"), y = unit(-0.12, "npc"),
            gp = gpar(fontsize = 24), hjust = 1
        )
    ) + # species label
    theme(
        plot.margin = margin(10, 10, 10, 60), # margins (t, r, b, l)
        panel.grid.minor = element_line(color = "gray"),
        axis.text = element_text(size = 26),
        axis.title = element_text(size = 30),
        plot.title = element_text(size = 35, hjust = 0.5),
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 24),
        legend.key.size = unit(1.5, "cm"),
        legend.position = "none"
    )

ooCount
ggsave(filename = "/Users/ivancasas/GitHub/Thesis/Chapters/03_SURV/pics/oocyst_counts.png", plot = ooCount, width = 12, height = 9)
ggsave(filename = "Figures/oocyst_counts.png", plot = ooCount, width = 12, height = 9)

# extract legend for plotting
ooCount_legend <- ooCount + theme(legend.position = "bottom")
legend_only <- cowplot::get_legend(ooCount_legend)
cowplot::save_plot(
    "/Users/ivancasas/GitHub/Thesis/Chapters/03_SURV/pics/colormap.png",
    legend_only,
    base_width = 10, # wide enough for horizontal legend
    base_height = 0.8 # short since it's just one row
)
cowplot::save_plot(
    "Figures/colormap.png",
    legend_only,
    base_width = 10, # wide enough for horizontal legend
    base_height = 0.8 # short since it's just one row
)

# stats
##################################################################################################################################

# view(ooprev_data) # from earlier
oocount_data <- ooprev_data %>% filter(oocysts > 0) # only positives midguts here

# format well
oocount_data$species <- as.factor(oocount_data$species)
oocount_data$temp_mean <- as.factor(oocount_data$temp_mean)
oocount_data$temp_range <- as.factor(oocount_data$temp_range)

oocount_data$temp_mean <- relevel(oocount_data$temp_mean, ref = "27")
oocount_data$temp_range <- relevel(oocount_data$temp_range, ref = "0")
oocount_data$species <- relevel(oocount_data$species, ref = "An. gambiae")
oocount_data$pot <- relevel(oocount_data$pot, ref = "ki276")

# all interactions
count_model <- glmmTMB(oocysts ~ species * temp_mean * temp_range + (1 | replicate) + (1 | pot_unique),
    data = oocount_data, family = truncated_nbinom2
) # poisson with overdispersion -> nbinom. No zeros -> truncated
summary(count_model) # we use this only for coefs. NOT pvalues

# there are two nbinoms in glmmTMB, this code below shows that we chose the right one above...
# count_model_nb1 <- glmmTMB(
#     oocysts ~ species * temp_mean * temp_range + (1 | replicate) + (1 | pot_unique),
#     data = oocount_data,
#     family = nbinom1
# )
# AIC(count_model, count_model_nb1) # low is good

# p computation - type iii analysis framework. we do it manually for consistency across models:
##############################################################################################

# reencode vars into sum-to-0 contrasts  (needed for type iii analysis)
local_contrasts <- list(
    species = contr.sum(2),
    temp_mean = contr.sum(2),
    temp_range = contr.sum(2)
)

# create model matrix manually with newly encoded vars
count_matrix <- model.matrix(~ species * temp_mean * temp_range, data = oocount_data, contrasts.arg = local_contrasts)[, -1]
colnames(count_matrix) <- gsub(":", "_", colnames(count_matrix)) # remove colons from colnames
count_data_sumto0 <- cbind(oocount_data, as.data.frame(count_matrix)) # bind to oo_count dats

# fit full model
count_model_sumto0 <- glmmTMB(oocysts ~ species * temp_mean * temp_range + (1 | replicate) + (1 | pot_unique), data = count_data_sumto0, family = truncated_nbinom2)

# sanity check:
logLik(count_model) # -407.9872 (df=11)
logLik(count_model_sumto0) # -407.9872 (df=11)
AIC(count_model) # 837.9745
AIC(count_model_sumto0) # 837.9745

# fit nested models for each term
no_sp <- glmmTMB(oocysts ~ temp_mean1 + temp_range1 + species1_temp_mean1 + species1_temp_range1 + temp_mean1_temp_range1 + species1_temp_mean1_temp_range1 + (1 | replicate) + (1 | pot_unique), data = count_data_sumto0, family = truncated_nbinom2)
no_tm <- glmmTMB(oocysts ~ species1 + temp_range1 + species1_temp_mean1 + species1_temp_range1 + temp_mean1_temp_range1 + species1_temp_mean1_temp_range1 + (1 | replicate) + (1 | pot_unique), data = count_data_sumto0, family = truncated_nbinom2)
no_tr <- glmmTMB(oocysts ~ species1 + temp_mean1 + species1_temp_mean1 + species1_temp_range1 + temp_mean1_temp_range1 + species1_temp_mean1_temp_range1 + (1 | replicate) + (1 | pot_unique), data = count_data_sumto0, family = truncated_nbinom2)
no_sp_tm <- glmmTMB(oocysts ~ species1 + temp_mean1 + temp_range1 + species1_temp_range1 + temp_mean1_temp_range1 + species1_temp_mean1_temp_range1 + (1 | replicate) + (1 | pot_unique), data = count_data_sumto0, family = truncated_nbinom2)
no_sp_tr <- glmmTMB(oocysts ~ species1 + temp_mean1 + temp_range1 + species1_temp_mean1 + temp_mean1_temp_range1 + species1_temp_mean1_temp_range1 + (1 | replicate) + (1 | pot_unique), data = count_data_sumto0, family = truncated_nbinom2)
no_tm_tr <- glmmTMB(oocysts ~ species1 + temp_mean1 + temp_range1 + species1_temp_mean1 + species1_temp_range1 + species1_temp_mean1_temp_range1 + (1 | replicate) + (1 | pot_unique), data = count_data_sumto0, family = truncated_nbinom2)
no_sp_tm_tr <- glmmTMB(oocysts ~ species1 + temp_mean1 + temp_range1 + species1_temp_mean1 + species1_temp_range1 + temp_mean1_temp_range1 + (1 | replicate) + (1 | pot_unique), data = count_data_sumto0, family = truncated_nbinom2)

p_table_count <- rbind(
    calc_lrt(count_model_sumto0, no_sp, "species"),
    calc_lrt(count_model_sumto0, no_tm, "temp_mean"),
    calc_lrt(count_model_sumto0, no_tr, "temp_range"),
    calc_lrt(count_model_sumto0, no_sp_tm, "species:temp_mean"),
    calc_lrt(count_model_sumto0, no_sp_tr, "species:temp_range"),
    calc_lrt(count_model_sumto0, no_tm_tr, "temp_mean:temp_range"),
    calc_lrt(count_model_sumto0, no_sp_tm_tr, "species:temp_mean:temp_range")
)

print(p_table_count)
