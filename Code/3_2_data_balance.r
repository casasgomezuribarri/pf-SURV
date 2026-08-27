# this script plots the balance of the data
# author: ivan casas


##################################################################################################################################
# Environment
##################################################################################################################################

source("Code/compile_data_survival.r") # sets the right wd and prepares the data nicely
source("Code/functions.r")

# load packages
packages <- c(
  "ggplot2",
  "gridExtra",
  "grid",
  "cowplot",
  "MASS",
  "Matching",
  "dagitty"
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

##################################################################################################################################
# plot
##################################################################################################################################

T <- c("exposed")
C <- c("temp_mean", "temp_range", "species")
C_labels <- c("Mean Temperature", "Temperature Range", "Species")

plots <- list()
for (i in seq_along(C)) {
  p <- ggplot(surv, aes_string(x = C[i], fill = T[1])) +
    geom_bar(position = "fill") +
    labs(title = "", x = C_labels[i], y = "Percentage") +
    theme_minimal(base_size = 17) +
    theme(
      legend.position = "none"
    )
  plots[[C[i]]] <- p
}

# get legend from one plot (before hiding it)
legend_plot <- ggplot(surv, aes_string(x = C[1], fill = T[1])) +
  geom_bar(position = "fill") +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 14)
  )
shared_legend <- get_legend(legend_plot)

# plots + shared legend
covariates <- grid.arrange(
  arrangeGrob(grobs = plots, ncol = 3),
  shared_legend,
  nrow = 2,
  heights = c(10, 1)
)
ggsave(filename = "/Users/ivancasas/GitHub/Thesis/Chapters/03_SURV/pics/distribution_confounding_factors.png", plot = covariates, width = 12, height = 5)
ggsave(filename = "Figures/distribution_confounding_factors.png", plot = covariates, width = 12, height = 5)
