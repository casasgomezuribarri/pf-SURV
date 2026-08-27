# analysis and plots for paraemtric survival
# author: ivan casas

##################################################################################################################################
# Environment (loading packages and data)
##################################################################################################################################

source("Code/compile_data_survival.r") # sets the right wd and calls view() on each dataset
source("Code/functions.r")

# load packages
packages <- c(
    "knitr",
    "dplyr",
    "survival",
    "emmeans",
    "ggplot2",
    "tibble",
    "devtools",
    "readr",
    "lubridate",
    "DT",
    "ggsurvfit",
    "gtsummary",
    "tidycmprsk",
    "RColorBrewer",
    "survminer",
    "coxme",
    "rms",
    "gridExtra",
    "flexsurv",
    "muhaz",
    "data.table",
    "readxl",
    "gridExtra"
)
for (i in packages) {
    if (!require(i, character.only = TRUE)) install.packages(i)
    library(i, character.only = TRUE)
}

surv$exposed <- as.factor(surv$exposed)
surv$species <- as.factor(surv$species)
surv$temp_mean <- as.factor(surv$temp_mean)
surv$temp_range <- as.factor(surv$temp_range)
surv$treatment <- as.factor(surv$treatment)
surv$age <- as.integer(surv$age)
surv$dpi <- as.integer(surv$dpi)
str(surv)

###################################################################################################################################
# parametric survival: all, effect of exposure
##################################################################################################################################

# flexsurvreg does not accept random effects

# choose a distribution
par_fit <- compare_parametric_fits(
    data = surv,
    time_var = "age",
    event_var = "dead",
    plot_title = paste0("Parametric fits")
)
# save plot
par_fit$plot
png(
    file = "/Users/ivancasas/GitHub/Thesis/Chapters/03_SURV/pics/parafits.png",
    width = 800, height = 800
)
par_fit$plot
dev.off()

par_fit$plot
png(
    file = "Figures/parafits.png",
    width = 800, height = 800
)
par_fit$plot
dev.off()

par_fit$comparison
# avoid GenF & Genγ (they're not very parsimonious options...)
best_dist <- "llogis"

# for interpreatbility, let's make ki270 the reference level...
surv$exposed <- relevel(surv$exposed, ref = "Control")
surv$temp_mean <- relevel(surv$temp_mean, ref = "27")
surv$temp_range <- relevel(surv$temp_range, ref = "0")
surv$species <- relevel(surv$species, ref = "An. gambiae")

m_exp_4 <- flexsurvreg(
    event ~ exposed * species * temp_mean * temp_range,
    data = surv, dist = best_dist
)

report_4 <- tidy(m_exp_4)
report_4 # use only for coefficient reporting
# view(report_4)

best_model <- m_exp_4

# p computation - type iii analysis framework. This is not supported in car pkg, so we do it manually:
##################################################################################################################################

# we'll need to reencode variables to sum-to-0
local_contrasts <- list(
    exposed = contr.sum(2),
    species = contr.sum(2),
    temp_mean = contr.sum(2),
    temp_range = contr.sum(2)
)

# using those contrasts, build the model matrix manually
X_matrix <- model.matrix(~ exposed * species * temp_mean * temp_range,
    data = surv,
    contrasts.arg = local_contrasts
)[, -1]

# gotta remove colons from col names
colnames(X_matrix) <- gsub(":", "_", colnames(X_matrix))

# bind to surv data
surv_data_sumto0 <- cbind(surv, as.data.frame(X_matrix)) # used to be surv_data_matrix


# full model formula (model.matrix warped names a bit but thats okay)
full_formula <- event ~ exposed1 + species1 + temp_mean1 + temp_range1 + # mains
    exposed1_species1 + exposed1_temp_mean1 + exposed1_temp_range1 + species1_temp_mean1 + species1_temp_range1 + temp_mean1_temp_range1 + # 2ways
    exposed1_species1_temp_mean1 + exposed1_species1_temp_range1 + exposed1_temp_mean1_temp_range1 + species1_temp_mean1_temp_range1 + # 3ways
    exposed1_species1_temp_mean1_temp_range1 # 4way

# this should be equivalent to the og model
best_model_sumto0 <- flexsurvreg(full_formula, data = surv_data_sumto0, dist = best_dist) # used to be full_surv

# sanity check that they are indeed the same model
logLik(best_model) # -16085.19 (df=17)
logLik(best_model_sumto0) # -16085.19 (df=17)
AIC(best_model) # 32204.39
AIC(best_model_sumto0) # 32204.39

# right. now fit every single relevant nested model (full model without 1 term)
no_ex <- flexsurvreg(event ~ species1 + temp_mean1 + temp_range1 + exposed1_species1 + exposed1_temp_mean1 + exposed1_temp_range1 + species1_temp_mean1 + species1_temp_range1 + temp_mean1_temp_range1 + exposed1_species1_temp_mean1 + exposed1_species1_temp_range1 + exposed1_temp_mean1_temp_range1 + species1_temp_mean1_temp_range1 + exposed1_species1_temp_mean1_temp_range1, data = surv_data_sumto0, dist = best_dist)
no_sp <- flexsurvreg(event ~ exposed1 + temp_mean1 + temp_range1 + exposed1_species1 + exposed1_temp_mean1 + exposed1_temp_range1 + species1_temp_mean1 + species1_temp_range1 + temp_mean1_temp_range1 + exposed1_species1_temp_mean1 + exposed1_species1_temp_range1 + exposed1_temp_mean1_temp_range1 + species1_temp_mean1_temp_range1 + exposed1_species1_temp_mean1_temp_range1, data = surv_data_sumto0, dist = best_dist)
no_tm <- flexsurvreg(event ~ exposed1 + species1 + temp_range1 + exposed1_species1 + exposed1_temp_mean1 + exposed1_temp_range1 + species1_temp_mean1 + species1_temp_range1 + temp_mean1_temp_range1 + exposed1_species1_temp_mean1 + exposed1_species1_temp_range1 + exposed1_temp_mean1_temp_range1 + species1_temp_mean1_temp_range1 + exposed1_species1_temp_mean1_temp_range1, data = surv_data_sumto0, dist = best_dist)
no_tr <- flexsurvreg(event ~ exposed1 + species1 + temp_mean1 + exposed1_species1 + exposed1_temp_mean1 + exposed1_temp_range1 + species1_temp_mean1 + species1_temp_range1 + temp_mean1_temp_range1 + exposed1_species1_temp_mean1 + exposed1_species1_temp_range1 + exposed1_temp_mean1_temp_range1 + species1_temp_mean1_temp_range1 + exposed1_species1_temp_mean1_temp_range1, data = surv_data_sumto0, dist = best_dist)
no_ex_sp <- flexsurvreg(event ~ exposed1 + species1 + temp_mean1 + temp_range1 + exposed1_temp_mean1 + exposed1_temp_range1 + species1_temp_mean1 + species1_temp_range1 + temp_mean1_temp_range1 + exposed1_species1_temp_mean1 + exposed1_species1_temp_range1 + exposed1_temp_mean1_temp_range1 + species1_temp_mean1_temp_range1 + exposed1_species1_temp_mean1_temp_range1, data = surv_data_sumto0, dist = best_dist)
no_ex_tm <- flexsurvreg(event ~ exposed1 + species1 + temp_mean1 + temp_range1 + exposed1_species1 + exposed1_temp_range1 + species1_temp_mean1 + species1_temp_range1 + temp_mean1_temp_range1 + exposed1_species1_temp_mean1 + exposed1_species1_temp_range1 + exposed1_temp_mean1_temp_range1 + species1_temp_mean1_temp_range1 + exposed1_species1_temp_mean1_temp_range1, data = surv_data_sumto0, dist = best_dist)
no_ex_tr <- flexsurvreg(event ~ exposed1 + species1 + temp_mean1 + temp_range1 + exposed1_species1 + exposed1_temp_mean1 + species1_temp_mean1 + species1_temp_range1 + temp_mean1_temp_range1 + exposed1_species1_temp_mean1 + exposed1_species1_temp_range1 + exposed1_temp_mean1_temp_range1 + species1_temp_mean1_temp_range1 + exposed1_species1_temp_mean1_temp_range1, data = surv_data_sumto0, dist = best_dist)
no_sp_tm <- flexsurvreg(event ~ exposed1 + species1 + temp_mean1 + temp_range1 + exposed1_species1 + exposed1_temp_mean1 + exposed1_temp_range1 + species1_temp_range1 + temp_mean1_temp_range1 + exposed1_species1_temp_mean1 + exposed1_species1_temp_range1 + exposed1_temp_mean1_temp_range1 + species1_temp_mean1_temp_range1 + exposed1_species1_temp_mean1_temp_range1, data = surv_data_sumto0, dist = best_dist)
no_sp_tr <- flexsurvreg(event ~ exposed1 + species1 + temp_mean1 + temp_range1 + exposed1_species1 + exposed1_temp_mean1 + exposed1_temp_range1 + species1_temp_mean1 + temp_mean1_temp_range1 + exposed1_species1_temp_mean1 + exposed1_species1_temp_range1 + exposed1_temp_mean1_temp_range1 + species1_temp_mean1_temp_range1 + exposed1_species1_temp_mean1_temp_range1, data = surv_data_sumto0, dist = best_dist)
no_tm_tr <- flexsurvreg(event ~ exposed1 + species1 + temp_mean1 + temp_range1 + exposed1_species1 + exposed1_temp_mean1 + exposed1_temp_range1 + species1_temp_mean1 + species1_temp_range1 + exposed1_species1_temp_mean1 + exposed1_species1_temp_range1 + exposed1_temp_mean1_temp_range1 + species1_temp_mean1_temp_range1 + exposed1_species1_temp_mean1_temp_range1, data = surv_data_sumto0, dist = best_dist)
no_ex_sp_tm <- flexsurvreg(event ~ exposed1 + species1 + temp_mean1 + temp_range1 + exposed1_species1 + exposed1_temp_mean1 + exposed1_temp_range1 + species1_temp_mean1 + species1_temp_range1 + temp_mean1_temp_range1 + exposed1_species1_temp_range1 + exposed1_temp_mean1_temp_range1 + species1_temp_mean1_temp_range1 + exposed1_species1_temp_mean1_temp_range1, data = surv_data_sumto0, dist = best_dist)
no_ex_sp_tr <- flexsurvreg(event ~ exposed1 + species1 + temp_mean1 + temp_range1 + exposed1_species1 + exposed1_temp_mean1 + exposed1_temp_range1 + species1_temp_mean1 + species1_temp_range1 + temp_mean1_temp_range1 + exposed1_species1_temp_mean1 + exposed1_temp_mean1_temp_range1 + species1_temp_mean1_temp_range1 + exposed1_species1_temp_mean1_temp_range1, data = surv_data_sumto0, dist = best_dist)
no_ex_tm_tr <- flexsurvreg(event ~ exposed1 + species1 + temp_mean1 + temp_range1 + exposed1_species1 + exposed1_temp_mean1 + exposed1_temp_range1 + species1_temp_mean1 + species1_temp_range1 + temp_mean1_temp_range1 + exposed1_species1_temp_mean1 + exposed1_species1_temp_range1 + species1_temp_mean1_temp_range1 + exposed1_species1_temp_mean1_temp_range1, data = surv_data_sumto0, dist = best_dist)
no_sp_tm_tr <- flexsurvreg(event ~ exposed1 + species1 + temp_mean1 + temp_range1 + exposed1_species1 + exposed1_temp_mean1 + exposed1_temp_range1 + species1_temp_mean1 + species1_temp_range1 + temp_mean1_temp_range1 + exposed1_species1_temp_mean1 + exposed1_species1_temp_range1 + exposed1_temp_mean1_temp_range1 + exposed1_species1_temp_mean1_temp_range1, data = surv_data_sumto0, dist = best_dist)
no_ex_sp_tm_tr <- flexsurvreg(event ~ exposed1 + species1 + temp_mean1 + temp_range1 + exposed1_species1 + exposed1_temp_mean1 + exposed1_temp_range1 + species1_temp_mean1 + species1_temp_range1 + temp_mean1_temp_range1 + exposed1_species1_temp_mean1 + exposed1_species1_temp_range1 + exposed1_temp_mean1_temp_range1 + species1_temp_mean1_temp_range1, data = surv_data_sumto0, dist = best_dist)

# custom function to show them all together
p_table_surv <- rbind(
    calc_lrt(best_model_sumto0, no_ex, "exposure"),
    calc_lrt(best_model_sumto0, no_sp, "species"),
    calc_lrt(best_model_sumto0, no_tm, "temp_mean"),
    calc_lrt(best_model_sumto0, no_tr, "temp_range"),
    calc_lrt(best_model_sumto0, no_ex_sp, "exposure:species"),
    calc_lrt(best_model_sumto0, no_ex_tm, "exposure:temp_mean"),
    calc_lrt(best_model_sumto0, no_ex_tr, "exposure:temp_range"),
    calc_lrt(best_model_sumto0, no_sp_tm, "species:temp_mean"),
    calc_lrt(best_model_sumto0, no_sp_tr, "species:temp_range"),
    calc_lrt(best_model_sumto0, no_tm_tr, "temp_mean:temp_range"),
    calc_lrt(best_model_sumto0, no_ex_sp_tm, "exposure:species:temp_mean"),
    calc_lrt(best_model_sumto0, no_ex_sp_tr, "exposure:species:temp_range"),
    calc_lrt(best_model_sumto0, no_ex_tm_tr, "exposure:temp_mean:temp_range"),
    calc_lrt(best_model_sumto0, no_sp_tm_tr, "species:temp_mean:temp_range"),
    calc_lrt(best_model_sumto0, no_ex_sp_tm_tr, "exposure:species:temp_mean:temp_range")
)

print(p_table_surv)
view(report_4) # og for comparison


# visualise predictions
##################################################################################################################################

# synthetic datasets for predictions
nd1 <- bind_rows(
    expand.grid(
        species = c("An. gambiae", "An. coluzzii"),
        temp_range = as.factor(c(0, 6)),
        temp_mean = as.factor(c(21, 27)),
        exposed = c("Control", "Exposed")
    )
)

# add counts for plot annotations
nd1$count <- numeric(nrow(nd1))
for (i in 1:nrow(nd1)) {
    # pot info
    sp <- nd1$species[i]
    tm <- nd1$temp_mean[i]
    tr <- nd1$temp_range[i]
    ex <- nd1$exposed[i]

    # grab pot data from surv
    pot <- subset(surv, species == sp & temp_mean == tm & temp_range == tr & exposed == ex)
    # count rows
    nd1$count[i] <- nrow(pot)
}


# sanity check - these should be the same
sum(nd1$count)
nrow(surv)

# make predictions on the new synthetic datasets
pred1 <- summary(best_model, newdata = nd1, type = "survival", ci = TRUE, tidy = TRUE)

# format for plotting
pred1 <- pred1 %>%
    mutate(
        treatment = paste0(temp_mean, "±", temp_range, "°C")
    )
my_palette <- brewer.pal(6, "Dark2")
selected_colors <- my_palette[c(1, 2)]


# kms to overlap with predictions (custom function for computing kms)
km1 <- get_km_data(surv, nd1, "exposed")

# calculate median survival per facet and group
median_survival <- km1 %>%
    group_by(species, treatment, exposed) %>%
    summarise(
        median_t = min(time[est <= 0.5], na.rm = TRUE),
        .groups = "drop"
    )

# add y positions for plotting
median_annotations <- median_survival %>%
    group_by(species, treatment) %>% # in each facet there's only 2 curves
    arrange(exposed) %>% # so we order alphabetically
    mutate(
        y_pos = case_when(
            row_number() == 1 ~ 0.20, # first one gets y = 0.2
            TRUE ~ 0.12 # second one (TRUE is 'all the rest') gets y =1.2
        )
    ) %>%
    ungroup() %>%
    mutate(label = paste0(median_t, " days"))

# probs good to add title to annotation
median_title <- median_annotations %>%
    distinct(species, treatment) %>%
    mutate(
        y_pos = 0.28,
        label = "Median:"
    )
# plot predictions + km + median annotations
parapreds1_km_ann <- ggplot(pred1, aes(x = time, y = est, colour = factor(exposed))) +
    geom_line(linewidth = 1) + # parametric fitted lines
    geom_ribbon( # confidence intervals
        aes(ymin = lcl, ymax = ucl, fill = factor(exposed)),
        alpha = 0.15, colour = NA
    ) +
    geom_step( # km curves in the same colours
        data = km1,
        aes(x = time, y = est, colour = factor(exposed)),
        linewidth = 1, inherit.aes = FALSE
    ) +
    geom_segment( # vertical lines from y=0.5 to y=0 at x=median_t
        data = median_survival,
        aes(
            x = median_t, xend = median_t,
            y = 0.5, yend = 0,
            colour = factor(exposed)
        ),
        linewidth = 0.8, linetype = "solid", inherit.aes = FALSE
    ) +
    geom_text( # 'median' titles (black font bold)
        data = median_title,
        aes(x = 1, y = y_pos, label = label),
        hjust = 0, vjust = 0,
        size = 6, fontface = "bold",
        colour = "black",
        inherit.aes = FALSE
    ) +
    geom_text( # median values in respective colours
        data = median_annotations,
        aes(x = 1, y = y_pos, label = label, colour = factor(exposed)),
        hjust = 0, vjust = 0,
        size = 6, fontface = "bold",
        inherit.aes = FALSE,
        show.legend = FALSE
    ) +
    scale_color_manual(values = selected_colors) +
    scale_fill_manual(values = selected_colors) +
    facet_grid(fct_rev(species) ~ treatment) + # coluzzii first because it's easier to explain and that helps the narrative
    scale_y_continuous(limits = c(0, 1)) +
    labs(
        title = "Parametric survival analysis",
        x = "Days post exposure",
        y = "Survival probability",
        colour = "P. falciparum status",
        fill = "P. falciparum status",
        caption = "Smooth = parametric fit; Step = Kaplan-Meier"
    ) +
    theme_minimal() +
    theme(
        panel.grid.minor = element_line(color = "white"),
        panel.grid = element_line(color = "gray"),
        axis.text = element_text(size = 26),
        strip.text = element_text(size = 30),
        plot.caption = element_text(size = 27, hjust = 0.5),
        strip.background = element_rect(fill = "#bdbdef", color = "white"),
        axis.title = element_text(size = 30),
        plot.margin = margin(10, 10, 10, 10),
        plot.title = element_text(size = 45, hjust = 0.5),
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 24),
        legend.key.size = unit(1.5, "cm"),
        legend.position = "bottom"
    )

parapreds1_km_ann

ggsave(plot = parapreds1_km_ann, "Figures/parametric_predictions_km_exposure_ann.png", width = 20, height = 12, units = "in", dpi = 150)
ggsave(plot = parapreds1_km_ann, "/Users/ivancasas/GitHub/Thesis/Chapters/03_SURV/pics/parametric_exposed_km_ann.png", width = 20, height = 12, units = "in", dpi = 150)

# actual coefficients (output is in log scale)
coefs <- tidy(best_model, conf.int = TRUE) %>%
    filter(!term %in% c("shape", "scale")) %>% # remove dist params
    mutate(
        TR = exp(estimate), # time ratio = exp(coef)
        TR_lo = exp(conf.low),
        TR_hi = exp(conf.high),
        term = factor(term, levels = rev(unique(term))) # make it a factor with levels in (rev) row order so plotting is in order (as opposed to alphabetically)
    )
view(coefs)

# plot parameter estimates
coefs_plot <- coefs %>%
    ggplot(aes(x = TR, y = term)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_errorbarh(aes(xmin = TR_lo, xmax = TR_hi), height = 0.2) +
    geom_point(size = 3, colour = "#1a6fa8") +
    scale_x_log10() +
    labs(x = "Time Ratio (log scale)", y = NULL) +
    theme_bw(base_size = 20)
coefs_plot

ggsave(plot = coefs_plot, "Figures/parametric_coefs.png", width = 15, height = 17, units = "in", dpi = 150)
ggsave(plot = coefs_plot, "/Users/ivancasas/GitHub/Thesis/Chapters/03_SURV/pics/parametric_coefs.png", width = 10, height = 12, units = "in", dpi = 150)
