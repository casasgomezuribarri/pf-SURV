# this script applies TAUS to mosquito data
# also plots the result overlapped with ocyst times
# TAUS: Time-Agnostic Unified Survival


##################################################################################################################################
# Environment (loading packages and data)
##################################################################################################################################
source("Code/compile_data_survival.r")

# load packages
packages <- c(
    "RColorBrewer",
    "betareg",
    "effectsize",
    "survminer",
    "nortest",
    "ggsurvfit",
    "cramer",
    "ggtext",
    "ggnewscale"
)

for (i in packages) {
    if (!require(i, character.only = TRUE)) install.packages(i)
    library(i, character.only = TRUE)
}

#################################################################################
# load the data, initialise useful variables
##################################################################################################################################

# view(surv)

# define a few useful things
var <- "pot" # variable of interest - try "treatment", some nice plots later
cats <- c("species", "exposed", "temp_mean", "temp_range", "treatment") # other categorical variables in the dataset (used for aggregating the data)
table(surv[[var]], useNA = "always") # check counts per level

time_var <- "dpi" # variable with time data
event_var <- "dead" # event variable
surv$event_factor <- factor(surv[[event_var]], levels = c(TRUE, FALSE), labels = c("Dead", "Censored")) # factor event variable for plotting
event_factor <- "event_factor"
max_time <- max(surv[[time_var]])
plot_names <- "phd_TAUS"

# quick descriptive plot of death and censoring times:
ggplot(surv, aes(x = !!sym(time_var), fill = event_factor)) +
    geom_histogram(position = "dodge", bins = 20, alpha = 1) +
    scale_fill_manual(values = c("Dead" = "red", "Censored" = "green")) +
    labs(
        title = "Histogram of survival times",
        x = "Time (days)",
        y = "Frequency",
        fill = "Exit reason"
    ) +
    facet_wrap(as.formula(paste("~", var)), ncol = 4) + # parse var
    # theme_minimal() +
    theme(
        plot.title = element_text(size = 30, face = "bold"),
        axis.title = element_text(size = 26),
        axis.text = element_text(size = 19),
        strip.text = element_text(size = 24),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 22)
    )
# ggsave(paste0("Figures/", plot_names, "_histograms.png"), width = 20, height = 20)

##################################################################################################################################
# conditional survival: calculate and visualise
##################################################################################################################################

# install.packages("remotes")
remotes::install_github("casasgomezuribarri/TAUS")
install.packages(".", repos = NULL, type = "source")
library(TAUS)

cond_surv <- cond_surv_mat(
    data = surv, # the dataset - a typical life table
    var = var, # variable of interest (string, must be categorical)
    cats = cats, # other categorical variables that matter (vector of strings, all categorical variables)
    time_var = time_var, # the time variable (string)
    event_var = event_var, # the event variable (string)
    res = 1, # resolution of the conditional survival grids for t and tau
    conf_int_level = 0.95, # confidence interval level (default 0.95)
    aggregate_by_cats = TRUE # if TRUE, will separate by each unique combination of values in var and cats (default FALSE)
)

##################################################################################################################################
# shaded plots of O_tau curves with significance shading
##################################################################################################################################

# make a dataframe with all valid pairwise comaprisons (when everything equal except exposure)
taus_analysis <- cond_surv$cond_surv %>%
    mutate(
        var_tau = as.factor(paste0(!!sym(var), "_", tau)) # unique id for later
    )

groups <- levels(taus_analysis$var_tau)

results <- data.frame(
    Group1 = character(),
    Group2 = character(),
    diff = numeric(),
    KS_p = numeric()
)

# this takes a minute (compares all relevant groups):
for (i in 1:length(groups)) {
    for (j in i:length(groups)) {
        # extract group names without exposure
        rest_1 <- paste0(substr(groups[i], 1, 1), substr(groups[i], 3, nchar(groups[i])))
        rest_2 <- paste0(substr(groups[j], 1, 1), substr(groups[j], 3, nchar(groups[j])))
        if (rest_1 != rest_2) {
            next # if different, skip
        }

        # extract exposure
        exp_1 <- substr(groups[i], 2, 2)
        exp_2 <- substr(groups[j], 2, 2)
        if (exp_1 == exp_2) {
            next # if they are the same, skip
        }

        # extrcat group and tau of interest
        unique_1 <- taus_analysis$unique_label[taus_analysis$var_tau == groups[i]][1]
        unique_2 <- taus_analysis$unique_label[taus_analysis$var_tau == groups[j]][1]

        tau_1 <- taus_analysis$tau[taus_analysis$var_tau == groups[i]][1]
        tau_2 <- taus_analysis$tau[taus_analysis$var_tau == groups[j]][1]


        stats <- pairwise_test(cond_surv,
            var_values = c(unique_1, unique_2) # the levels of the variable of interest to compare
            , tau_values = c(tau_1, tau_2) # the ages beyond which survival is interesting
        )
        results <- rbind(results, data.frame(
            Group1 = groups[i],
            Group2 = groups[j],
            diff = stats$effect_size,
            p_value = stats$p_value
        ))
    }
}

results <- results %>% # add species info columns)
    mutate(
        species = case_when(
            str_sub(Group1, 1, 1) == "k" & str_sub(Group2, 1, 1) == "k" ~ "An. gambiae",
            str_sub(Group1, 1, 1) == "c" & str_sub(Group2, 1, 1) == "c" ~ "An. coluzzii",
            TRUE ~ NA_character_
        ),
        tau = case_when(
            str_sub(Group1, 7, nchar(Group1)) == str_sub(Group2, 7, nchar(Group2)) ~ str_sub(Group2, 7, nchar(Group2)),
            TRUE ~ NA_character_
        ),
        treatment = case_when(
            str_sub(Group1, 3, 5) == "210" & str_sub(Group2, 3, 5) == "210" ~ "21±0°C",
            str_sub(Group1, 3, 5) == "216" & str_sub(Group2, 3, 5) == "216" ~ "21±6°C",
            str_sub(Group1, 3, 5) == "270" & str_sub(Group2, 3, 5) == "270" ~ "27±0°C",
            str_sub(Group1, 3, 5) == "276" & str_sub(Group2, 3, 5) == "276" ~ "27±6°C",
            TRUE ~ NA_character_
        )
    )

results$tau <- as.numeric(results$tau) # convert tau to numeric
# view(results)

# add to cond_surv for plotting
taus_analysis <- taus_analysis %>%
    left_join(results, by = c("species", "tau", "treatment"))

taus_analysis <- taus_analysis %>%
    mutate(
        signif = ifelse(p_value < 0.05, "Significant", "Not significant"),
        diff_thresh = ifelse(abs(diff) > 0.01, "Very different", "Not very different") # playing around to see if this matters much. Not really
    ) %>%
    mutate(relevant = ifelse(signif == "Significant", "p<0.05", "p>0.05"))

# compute oocyst phase for the horizontal bars
# view(replicates) # replicates is generated in the data compiler
oocyst_times <- replicates[, 1:4] %>% # only the first 4 cols
    mutate( # add species and treatment for facets
        species = case_when(
            str_sub(pot, 1, 1) == "k" ~ "An. gambiae",
            str_sub(pot, 1, 1) == "c" ~ "An. coluzzii",
            TRUE ~ NA_character_
        ),
        treatment = case_when(
            str_sub(pot, 3, 5) == "210" ~ "21±0°C",
            str_sub(pot, 3, 5) == "216" ~ "21±6°C",
            str_sub(pot, 3, 5) == "270" ~ "27±0°C",
            str_sub(pot, 3, 5) == "276" ~ "27±6°C",
            TRUE ~ NA_character_
        )
    ) %>%
    group_by(pot, species, treatment) %>%
    summarise(
        min_t_o = min(t_o),
        max_t_s = max(t_s, na.rm = TRUE)
    )

# view(oocyst_times)

# plot things
custom_labeller <- c("0" = "Constant temp", "6" = "Oscillating temp")
my_palette <- brewer.pal(6, "Dark2")
selected_colors <- my_palette[c(1, 2)]

combined_plot <- ggplot(taus_analysis, aes(x = tau, y = O_tau, color = exposed)) +
    # geom_tile for colouring the background according to pvalues
    geom_tile(
        aes(x = tau, y = O_tau, fill = factor(relevant)),
        width = 1, height = Inf, alpha = 0.002, inherit.aes = FALSE
    ) + # horizontal bars to mark time between first oocyst and last sporozoite observation
    geom_segment(
        data = oocyst_times,
        aes(x = min_t_o, xend = max_t_s, y = 0, yend = 0),
        linewidth = 2,
        inherit.aes = FALSE
    ) + # aesthetics
    scale_fill_manual(
        values = c("p<0.05" = "red", "p>0.05" = "transparent"),
        name = "Significance",
        na.translate = FALSE, # NA values (greay) where only one curve exists
        guide = guide_legend(
            order = 2, # we weant this legend to be the second to appear
            nrow = 2,
            override.aes = list(alpha = 0.4, color = "grey50", linewidth = 0.5)
        )
    ) +
    # the curves will use a different legend so initiate it
    ggnewscale::new_scale_fill() +
    # the curves. line + ribbon
    geom_line(size = 1) + # inherits ggplot aes no need to specify
    geom_ribbon(aes(ymin = O_tau_lo, ymax = O_tau_up, fill = exposed), alpha = 0.2) +
    xlim(0, max_time) +
    ylim(0, 1) +
    # colours for the curves and ribbons
    scale_color_manual(
        values = selected_colors,
        name = "P. falciparum status",
        labels = c("Control", "Exposed"),
        guide = guide_legend(
            order = 1,
            nrow = 2,
            override.aes = list(color = selected_colors)
        )
    ) +
    scale_fill_manual(
        values = selected_colors,
        name = "P. falciparum status",
        labels = c("Control", "Exposed"),
        guide = guide_legend(
            order = 1, # first legend
            nrow = 2,
            override.aes = list(alpha = 0.3, color = "grey50", linewidth = 0.5) # doesn't wanna colour the outline grey, ok...
        )
    ) +
    facet_grid(
        species ~ treatment,
        labeller = labeller(treatment = custom_labeller)
    ) +
    labs(
        x = "Days post exposure, τ",
        y = "Probability of outliving τ",
        title = "TAUS survival analysis"
    ) +
    theme_minimal() +
    theme(
        panel.grid.minor = element_line(color = "white"), # Customize minor grid lines
        panel.grid = element_line(color = "gray"),
        axis.text = element_text(size = 26), # Font size for axis ticks
        strip.text = element_text(size = 30),
        plot.caption = element_text(size = 27, hjust = 0.5),
        strip.background = element_rect(fill = "#bdbdef", color = "white"),
        axis.title = element_text(size = 30), # Adjust font of labels
        plot.margin = margin(10, 10, 10, 10), # Plot margins (t, r, b, l)
        plot.title = element_text(size = 45, hjust = 0.5), # Title settings
        legend.title = element_text(size = 30), # Font size for legend title
        legend.text = element_text(size = 24), # Font size for legend text
        legend.key.size = unit(1.5, "cm"), # Size of legend keys
        legend.position = "bottom"
    )
combined_plot

ggsave(plot = combined_plot, paste0("/Users/ivancasas/GitHub/Thesis/Chapters/03_SURV/pics/", plot_names, "_P(T>tau)_panel_shaded.png"), width = 20, height = 12, units = "in", dpi = 150)

##################################################################################################################################
# An example of TAUS for statistical analysis of survival
##################################################################################################################################
# # same tau for all levels of variable of interest

# stat_results <- pairwise_test(cond_surv,
#   tau_values = tau_value # the age beyond which survival is interesting
# )

# stat_results

# using TAUS to compare group survival to different ages
groups <- unique(cond_surv$cond_surv$unique_label) # these are all the groups than can be compared

# choose groups to compare
var_values <- c(groups[1], groups[3]) # let's compare unexposed gambaie and coluzzii at standard insectary conditions
tau_values <- c(10, 15) # let's compare their survival to these ages (assigned in the same order as var_values)

stat_results <- pairwise_test(cond_surv, # see Ln78
    var_values = var_values # the levels of the variable of interest to compare
    , tau_values = tau_values # the ages beyond which survival is interesting
)

stat_results


# how to interpret the results:

# group_n and tau_n: group and group-specific tau being compared
# O_tau_n = probability of outliving tau_n for a randomly selected individual from group_n
# effect_size = O_tau_1 - O_tau_2
#   an individual from group_1 has effect_size more chance of outliving tau_1 than an individual from group_2 of outliving tau_2
# effect_ratio = O_tau_1 / O_tau_2
#   an individual from group_1 has effect_ratio times the chance of outliving tau_1 than an individual from group_2 of outliving tau_2


# ks statistic - for comparing the O_tau_1 and O_tau_2 values:
# - 0: both are identical
# - 1: both are completely different

# ks statistic is analytically calculated as follows:
# 1. Betas are fitted to both O_tau (±95%CI) values, estimating alpha and beta parameters
# 2. ks = max(F1(x) - F2(x)) for all x in [0, 1], where F1 and F2 are the CDFs of the two Betas

# p-value is likelihood of H0 being true (both distributions are similar)
# calculated as follows
# 1. A null Beta is estimated by averaging the alpha and beta of the two betas fitted in the calculation of the ks statistic
# 2. n_draw (=10) pairs of values are sampled from the null Beta
# 3. the ks distance from those two sets of values is estimated and stored
# 4. Steps 3-4 are repeated B (=5000) times
# 5. p-value is the proportion of ks distances that are greater than the analytic ks.

# p-value answers these questions:
# how often is the KS statistic from the null pairs greater than or equal to the deterministic one?
# or, how likely is it to see a deterministic distance this big if they actually come from the same generative process?
# or, how likely is it that these two O(tau) values are fundamentally similar?
