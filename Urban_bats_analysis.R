#------------------------------------------------------------------------------
# Analysis of insectivorous bat species-specific activity and species diversity
# Please ensure this script and all associated CSV files are located in the same 
# working directory, or open the provided .Rproj file before running.
#-------------------------------------------------------------------------------

# For reproducibility
set.seed(123)

library(tidyverse)
library(EcoCountHelper)
library(ggeffects)
library(patchwork)
library (broom.mixed)
library(corrplot)
library(performance)
library(see)
library(vegan)

# 1. Data cleaning
# Importing the data containing bat passes extracted from Raven Pro annotations
# The bat passes are assigned as per site their respective variables
Bat_data <- read.csv("Raw data_bat passes.csv")

# Importing location data
Sampled_locations <- read.csv("Sampled_locations.csv")

# Importing summary data and dropping the diversity metrics we don't need at the moment
call_summary_locations <- read.csv("Season wise bat community data.csv") %>%
  dplyr::select(-Shannon_index, -Species_richness, -Simpson_index)

# Removing any white spaces that exist
Sampled_locations$Site      <- trimws(Sampled_locations$Site)
call_summary_locations$Site <- trimws(call_summary_locations$Site)

# Selecting the required habitat metrics from the sampled locations data
corrected_habitat_data <- Sampled_locations %>%
  dplyr::select(Site, 
                Vegetation_percentage, 
                Built_up_percentage, 
                Open_barren_percentage, 
                Tar_road_percentage, 
                Unclassified_habitat_percentage) %>%
  # Converting NA, if any exist, to zero
  dplyr::mutate(across(ends_with("_percentage"), ~tidyr::replace_na(., 0))) %>%
  dplyr::distinct(Site, .keep_all = TRUE)

# Bringing in proper habitat metrics from corrected habitat data
# Ensuring there is no duplication in the process
call_summary_locations <- call_summary_locations %>%
  dplyr::select(-Vegetation_percentage, 
                -Built_up_percentage, 
                -Open_barren_percentage, 
                -Tar_road_percentage) %>%
  dplyr::left_join(corrected_habitat_data, by = "Site")

# Creating factors for season and renaming proper four letter codes for the bat species
species_data <- Bat_data %>%
  dplyr::mutate(
    Season = factor(Season, levels = c("winter", "non-winter")),
    Species = dplyr::case_when(
      # Alionoctula ceylonica / Scotophilus heathii
      Species == "PICE" ~ "ACSH",
      
      # Pipistrellus tenuis to Alionoctula tenuis
      Species == "PITE" ~ "ALTE",
      
      # Nyctinomus aegyptiacus / Mops plicatus
      Species == "TAAE" ~ "NAMP",
      Species == "MOPL" ~ "NAMP",
      
      # Taphozous melanopogon/ Taphozous nudiventris
      Species == "TAME" ~ "TAPH",
      Species == "TANU" ~ "TAPH",
      
      TRUE ~ Species 
    )
  )

# Aggregating total calls at the site-season-species level
species_summary <- species_data %>%
  dplyr::group_by(Site, Season, Species) %>%
  dplyr::summarise(
    Total_Calls = dplyr::n(),
    .groups = "drop"
  )

cat("Rows in species_summary:", nrow(species_summary), "\n")

# Filtering for four modellable species
species_model_data <- species_summary |> 
  dplyr::filter(Species %in% c( "ALTE", "ACSH", "TAPH", "NAMP"))

# Creating a combination of all site, season and species
# Full grid to make sure "zero" counts are included
actual_combinations <- call_summary_locations %>%
  dplyr::select(Site, Season) |> 
  dplyr::distinct() |> 
  dplyr::distinct() %>%
  dplyr::cross_join(
    data.frame(
      Species = c("ACSH", "ALTE", "TAPH", "NAMP"), # Using updated clean codes
      stringsAsFactors = FALSE
    )
  )

cat("Total rows after crossing with species:", nrow(actual_combinations), "\n")

# Collinearity testing
# Subsetting environmental and climatic variables
collinearity_data <- call_summary_locations %>%
  dplyr::select(Vegetation_percentage, Tar_road_percentage, Built_up_percentage,
                Mean_lux, Moon_percentage, Temperature, Season) %>%
  dplyr::distinct() %>%
  tidyr::drop_na() 

# Creating pearson correlation matrix
# Removing Season since it's a categorical variable
continuous_vars <- collinearity_data %>%
  dplyr::select(-Season)

cor_matrix <- cor(continuous_vars, use = "complete.obs", method = "pearson")

# Proceeding to plot
corrplot::corrplot(cor_matrix, method = "color", type = "lower", 
                   addCoef.col = "black", tl.col = "black", tl.srt = 45, 
                   diag = FALSE, title = "Pearson Correlation Matrix", mar = c(0,0,1,0))


# Testing VIF (Variance inflation factor)
# Creating a model containing all variables including season
vif_model <- lm(Total_Calls ~ Vegetation_percentage + Tar_road_percentage + 
                  Built_up_percentage + Mean_lux + Moon_percentage + 
                  Temperature + Season, data = call_summary_locations)

# Checking VIF scores of different variables
vif_results <- performance::check_collinearity(vif_model)
print(vif_results)

# Plotting and saving the diagnostic results
vif_plot <- plot(vif_results) + 
  theme_bw() +
  labs(title = NULL, subtitle = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

ggsave("Predictor_VIF_Plot.png", vif_plot, width = 8, height = 6, dpi = 600)

# We omit temperature because of its high collinearity with season

# Constructing the final modelling dataset
Species_model_input <- call_summary_locations |> 
  dplyr::select(Site, Season, Vegetation_percentage, Tar_road_percentage,
                Mean_lux, Built_up_percentage, Moon_percentage) |> 
  # Ensuring no duplicate rows exist
  dplyr::distinct() |> 
  # Bringing in species per every site visit
  dplyr::cross_join(
    data.frame(Species = c("ALTE", "ACSH", "TAPH", "NAMP"))
    ) |> 
  # Merging actual bat counts into the grid
  dplyr::left_join(
    species_summary |> 
      dplyr::select(Site, Season, Species, Total_Calls),
      by = c("Site", "Season", "Species")
  ) |> 
  # Adding zeros for true non detections
  dplyr::mutate(Total_Calls = tidyr::replace_na(Total_Calls, 0))

# Creating factors
Species_model_input <- Species_model_input %>%
  dplyr::mutate(
    # Season as factors
    Season = factor(Season, levels = c("winter", "non-winter")),
    
    # Treating Built_up percentage as a factor
    Builtup_cat = factor(
      dplyr::case_when(
        Built_up_percentage == 0 ~ "zero built-up",
        Built_up_percentage > 0  ~ "built-up"
      ),
      levels = c("zero built-up", "built-up")
    ),
    
    # Treating light as factor at 1 lux
    Lux_cat_1 = factor(
      dplyr::case_when(
        Mean_lux < 1  ~ "Low",
        Mean_lux >= 1 ~ "High"
      ),
      levels = c("Low", "High")
    ),
    
    #Treating light as factor at 5 lux
    Lux_cat_5 = factor(
      dplyr::case_when(
        Mean_lux < 5  ~ "Low",
        Mean_lux >= 5 ~ "High"
      ),
      levels = c("Low", "High")
    )
  )

cat("\n=== Lux threshold < 5 ===\n")
table(Species_model_input$Lux_cat_5)
cat("\n=== Lux threshold < 1 ===\n")
table(Species_model_input$Lux_cat_1)
cat("\n=== Builtup threshold 0 ===\n")
table(Species_model_input$Builtup_cat)


# Checking zeros per species 
Species_model_input %>%
  dplyr::group_by(Species) %>%
  dplyr::summarise(
    N         = dplyr::n(),
    Zeros     = sum(Total_Calls == 0),
    NonZeros  = sum(Total_Calls > 0),
    Pct_zeros = round(sum(Total_Calls == 0) / dplyr::n() * 100, 1),
    Mean      = round(mean(Total_Calls), 2),
    Variance  = round(var(Total_Calls), 2),
    Ratio     = round(var(Total_Calls) / mean(Total_Calls), 2), 
    .groups   = "drop"
  )

# Proceeding with zero inflation testing
# Creating a function which outputs whether a zero inflation model is required
# W.r.t each sonotype
selected_species <- c( "ALTE", "ACSH", "TAPH", "NAMP")

for(sp in selected_species) {
  cat("\n----------------------------------------\n")
  cat(" Zero-Inflation Test for:", sp, "\n")
  sp_data <- Species_model_input |> 
    dplyr::filter(Species == sp)
  
  # Model A: Standard Negative Binomial
  m <- glmmTMB(Total_Calls ~ Vegetation_percentage + Tar_road_percentage +
                 Lux_cat_1 + Builtup_cat + Moon_percentage + Season,
               family = nbinom2(), data = sp_data)
  
  # Model B: Zero-Inflated Negative Binomial
  m_zi <- glmmTMB(Total_Calls ~ Vegetation_percentage + Tar_road_percentage +
                    Lux_cat_1 + Builtup_cat + Moon_percentage + Season,
                  ziformula = ~1, family = nbinom2(), data = sp_data)
  
  # Simulating residuals to run the test
  sim     <- DHARMa::simulateResiduals(m, n = 1000)
  zi_test <- DHARMa::testZeroInflation(sim)
  
  cat("DHARMa ZI p-value:     ", round(zi_test$p.value, 3), "\n")
  cat("AIC Standard NB2:      ", round(AIC(m), 2), "\n")
  cat("AIC Zero-Inflated NB2: ", round(AIC(m_zi), 2), "\n")
  cat("ZI Model Required?     ", ifelse(AIC(m_zi) < AIC(m) - 2, "YES", "NO"), "\n")
}
# Therefore, in accordance with the results
# None of the sonotypes need a zero inflated model


# Next, Candidate model selection by testing models with families nbinom1 and nbinom2:
# All models are global models

# Pip. ceylonicus / Scot. heathii
# Modelling using site as a random effect
ACSH_data <- Species_model_input %>% dplyr::filter(Species == "ACSH")

ACSH_complete_Nb2 <- glmmTMB(Total_Calls ~ Vegetation_percentage + Tar_road_percentage + 
                               Lux_cat_1 + Builtup_cat + Moon_percentage + Season + (1|Site), 
                             family = nbinom2(), data = ACSH_data)

ACSH_complete_Nb1 <- glmmTMB(Total_Calls ~ Vegetation_percentage + Tar_road_percentage + 
                               Lux_cat_1 + Builtup_cat + Moon_percentage + Season + (1|Site), 
                             family = nbinom1(), data = ACSH_data)

# Alionoctula tenuis
ALTE_data <- Species_model_input %>% dplyr::filter(Species == "ALTE")

# Using fixed effects because st.dev when using Random effect was 7.0559e-05
# Therefore using random effects makes no difference
ALTE_complete_Nb2_fixed <- glmmTMB(Total_Calls ~ Vegetation_percentage + Tar_road_percentage + 
                                     Lux_cat_1 + Builtup_cat + Moon_percentage + Season, 
                                   family = nbinom2(), data = ALTE_data)

ALTE_complete_Nb1_fixed <- glmmTMB(Total_Calls ~ Vegetation_percentage + Tar_road_percentage + 
                                     Lux_cat_1 + Builtup_cat + Moon_percentage + Season, 
                                   family = nbinom1(), data = ALTE_data)

# Taphozous spp.
TAPH_data <- Species_model_input %>% dplyr::filter(Species == "TAPH")

TAPH_complete_Nb2 <- glmmTMB(Total_Calls ~ Vegetation_percentage + Tar_road_percentage + 
                               Lux_cat_1 + Builtup_cat + Moon_percentage + Season + (1|Site), 
                             family = nbinom2(), data = TAPH_data)

TAPH_complete_Nb1 <- glmmTMB(Total_Calls ~ Vegetation_percentage + Tar_road_percentage + 
                               Lux_cat_1 + Builtup_cat + Moon_percentage + Season + (1|Site), 
                             family = nbinom1(), data = TAPH_data)

# Nyctinomops aegyptiacus/Mops plicatus
NAMP_data <- Species_model_input %>% dplyr::filter(Species == "NAMP")

# Using fixed effects because st.dev when using Random effect was 0.00010865
# Therefore using random effect makes no difference
NAMP_complete_Nb2_fixed <- glmmTMB(Total_Calls ~ Vegetation_percentage + Tar_road_percentage + 
                                     Lux_cat_1 + Builtup_cat + Moon_percentage + Season, 
                                   family = nbinom2(), data = NAMP_data)

NAMP_complete_Nb1_fixed <- glmmTMB(Total_Calls ~ Vegetation_percentage + Tar_road_percentage + 
                                     Lux_cat_1 + Builtup_cat + Moon_percentage + Season, 
                                   family = nbinom1(), data = NAMP_data)

# Proceeding to model selection as per EcoCountHelper package - Cole et al. 2022:
Groups <- c("ACSH_complete", "ALTE_complete", "TAPH_complete", "NAMP_complete")
ModelCompare(Groups, TopModOutName = "CompleteBestMods")

# AIC breakdown per sonotype 
ACSH_completeAIC
ALTE_completeAIC
TAPH_completeAIC
NAMP_completeAIC
# AIC wise, nbinom2 is better than nbinom1 for all four models

# Combined list of models
CompleteBestMods <- list(
  ACSH = ACSH_complete_Nb2,
  ALTE = ALTE_complete_Nb2_fixed,
  TAPH = TAPH_complete_Nb2,
  NAMP = NAMP_complete_Nb2_fixed
)

# Combined list
print(CompleteBestMods)
summary(ACSH_complete_Nb2)
summary(ALTE_complete_Nb2_fixed)
summary(TAPH_complete_Nb2)
summary(NAMP_complete_Nb2_fixed)


# Testing residuals of the best models using DHarma package
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

sim_ACSH   <- DHARMa::simulateResiduals(ACSH_complete_Nb2,        n = 1000, plot = TRUE)
sim_ALTE   <- DHARMa::simulateResiduals(ALTE_complete_Nb2_fixed,  n = 1000, plot = TRUE)
sim_TAPH   <- DHARMa::simulateResiduals(TAPH_complete_Nb2,        n = 1000, plot = TRUE)
sim_NAMP   <- DHARMa::simulateResiduals(NAMP_complete_Nb2_fixed,  n = 1000, plot = TRUE)

# Extracting Incidence Ratio Rates (IRR)
# IRRs determine the proportional multiplier effect of each predictor on bat activity
extract_irr <- function(model, species_name) {
  broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE, exponentiate = TRUE) %>%
    dplyr::filter(term != "(Intercept)") %>%
    dplyr::mutate(
      Species = species_name,
      term = dplyr::case_when(
        term == "Vegetation_percentage" ~ "Vegetation %",
        term == "Tar_road_percentage"   ~ "Tar Road %",
        term == "Lux_cat_1High"          ~ "Light Level (High)",
        term == "Builtup_catbuilt-up"   ~ "Built-up (Present)",
        term == "Moon_percentage"       ~ "Moon %",
        term == "Seasonnon-winter"      ~ "Season (Non-Winter)",
        TRUE ~ term
      ),
      IRR     = round(estimate, 3),
      Lower   = round(conf.low, 3),
      Upper   = round(conf.high, 3),
      p_value = round(p.value, 3),
      Sig = dplyr::case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01  ~ "**",
        p.value < 0.05  ~ "*",
        p.value < 0.1   ~ ".",
        TRUE            ~ ""
      )
    ) %>%
    dplyr::select(Species, Predictor = term, IRR, Lower, Upper, p_value, Sig)
}

# Combining the best models into a single dataframe
irr_complete <- dplyr::bind_rows(
  extract_irr(ACSH_complete_Nb2, "ACSH"),
  extract_irr(ALTE_complete_Nb2_fixed, "ALTE"),
  extract_irr(TAPH_complete_Nb2, "TAPH"),
  extract_irr(NAMP_complete_Nb2_fixed, "NAMP")
)

# Viewing and saving the IRRs
print(irr_complete)
write.csv(irr_complete, "IRR_complete_models.csv", row.names = FALSE)

# Proceeding to calulate R2 values
# 1. Nakagawa R2 for Mixed Models
r2_ACSH <- performance::r2_nakagawa(ACSH_complete_Nb2)
r2_TAPH <- performance::r2_nakagawa(TAPH_complete_Nb2)

# 2. Fit Intercept-Only Null Models for Fixed-Effect Species
null_ALTE <- glmmTMB(Total_Calls ~ 1, family = nbinom2(), data = ALTE_data)
null_NAMP <- glmmTMB(Total_Calls ~ 1, family = nbinom2(), data = NAMP_data)

# 3. Compile Hybrid Dataframe safely
r2_summary <- data.frame(
  Species = c("ACSH", "ALTE", "TAPH", "NAMP"),
  
  R2_Marginal = c(
    round(r2_ACSH$R2_marginal, 3),
    round(1 - as.numeric(logLik(ALTE_complete_Nb2_fixed)) / as.numeric(logLik(null_ALTE)), 3),
    round(r2_TAPH$R2_marginal, 3),
    round(1 - as.numeric(logLik(NAMP_complete_Nb2_fixed)) / as.numeric(logLik(null_NAMP)), 3)
  ),
  
  R2_Conditional = c(
    round(r2_ACSH$R2_conditional, 3),
    NA, # Fixed GLM - No random effects partition
    round(r2_TAPH$R2_conditional, 3),
    NA  # Fixed GLM - No random effects partition
  )
)

# Looking at the summaries
print(r2_summary)
# Saving the summaries
write.csv(r2_summary, "R2_complete_models.csv", row.names = FALSE)


# Visualizing
# Creating a function for plotting
pred_plot_all <- function(model, predictor, x_label, species_name, 
                          sig = FALSE, bias_correction = FALSE) {
  
  is_cat <- predictor %in% c("Lux_cat_1", "Builtup_cat", "Season")
  
  # Pass the bias_correction parameter down to handle mixed model population-level predictions
  pred <- ggeffects::ggpredict(model, terms = predictor, bias_correction = bias_correction)
  
  # Adding asterisk to significant predictors in the plot title
  display_title <- if(sig) paste0(species_name, "*") else species_name
  
  if(is_cat) {
    p <- ggplot(pred, aes(x = x, y = predicted)) +
      geom_errorbar(aes(ymin = conf.low, ymax = conf.high), 
                    width = 0.15, colour = "grey40") +
      geom_point(size = 2.5)
  } else {
    p <- ggplot(pred, aes(x = x, y = predicted)) +
      geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
                  fill = "grey80", alpha = 0.6) +
      geom_line(linewidth = 0.8)
  }
  
  p +
    labs(title = display_title, x = x_label, y = "Predicted Total Calls") +
    theme_bw() +
    theme(
      plot.title       = element_text(face = "plain", hjust = 0.5, size = 12),
      panel.grid.minor = element_blank(),
      axis.title.x     = element_text(size = 14),
      axis.title.y     = element_text(size = 14),
      axis.text        = element_text(size = 12)
    )
}

# 1. Plots predicting the effects of vegetation percentage on each sonotype
p_veg_ACSH   <- pred_plot_all(ACSH_complete_Nb2, "Vegetation_percentage", 
                              "Vegetation %", "ACSH", sig = FALSE, bias_correction = TRUE) + 
  coord_cartesian(ylim = c(0, 125)) + scale_y_continuous(breaks = seq(0, 125, 25))
p_veg_ALTE   <- pred_plot_all(ALTE_complete_Nb2_fixed, "Vegetation_percentage", 
                              "Vegetation %", "ALTE", sig = FALSE) + 
  coord_cartesian(ylim = c(0, 125)) + scale_y_continuous(breaks = seq(0, 125, 25))
p_veg_TAPH   <- pred_plot_all(TAPH_complete_Nb2, "Vegetation_percentage", 
                              "Vegetation %", "TAPH", sig = FALSE, bias_correction = TRUE) + 
  coord_cartesian(ylim = c(0, 125)) + scale_y_continuous(breaks = seq(0, 125, 25))
p_veg_NAMP   <- pred_plot_all(NAMP_complete_Nb2_fixed, "Vegetation_percentage", 
                              "Vegetation %", "NAMP", sig = FALSE) + 
  coord_cartesian(ylim = c(0, 125)) + scale_y_continuous(breaks = seq(0, 125, 25))

# Creating a row of 4 plots representing each sonotype
p_veg_row <- wrap_plots(p_veg_ACSH, p_veg_ALTE, p_veg_TAPH, p_veg_NAMP, nrow = 1) +
  plot_layout(axes = "collect") +
  plot_annotation(tag_levels = list(c("A", "B", "C", "D"))) &
  theme(plot.tag = element_text(size = 10, face = "plain", hjust = 0.5, vjust = 0), 
        plot.tag.position = c(0.05, 1))

# Viewing the plot
p_veg_row

# Saving the plot
ggsave("pred_veg.png", p_veg_row, width = 12, height = 6, dpi = 600)


# 2. Plots predicting the effects of Tar percentage on each sonotype
p_tar_ACSH   <- pred_plot_all(ACSH_complete_Nb2, "Tar_road_percentage", 
                              "Tar Road %", "ACSH", sig = FALSE, bias_correction = TRUE) + 
  coord_cartesian(ylim = c(0, 100)) + scale_y_continuous(breaks = seq(0, 100, 25))

p_tar_ALTE   <- pred_plot_all(ALTE_complete_Nb2_fixed, "Tar_road_percentage", 
                              "Tar Road %", "ALTE", sig = TRUE) + 
  coord_cartesian(ylim = c(0, 100)) + scale_y_continuous(breaks = seq(0, 100, 25))

p_tar_TAPH   <- pred_plot_all(TAPH_complete_Nb2, "Tar_road_percentage", 
                              "Tar Road %", "TAPH", sig = FALSE, bias_correction = TRUE) +
  coord_cartesian(ylim = c(0, 100)) + scale_y_continuous(breaks = seq(0, 100, 25))

p_tar_NAMP   <- pred_plot_all(NAMP_complete_Nb2_fixed, "Tar_road_percentage", 
                              "Tar Road %", "NAMP", sig = FALSE) + 
  coord_cartesian(ylim = c(0, 100)) + scale_y_continuous(breaks = seq(0, 100, 25))

# Creating a row of 4 plots representing each sonotype
p_tar_row <- wrap_plots(p_tar_ACSH, p_tar_ALTE, p_tar_TAPH, p_tar_NAMP, nrow = 1) +
  plot_layout(axes = "collect") +
  plot_annotation(tag_levels = list(c("E", "F", "G", "H"))) &
  theme(plot.tag = element_text(size = 10, face = "plain", hjust = 0.5, vjust = 0), 
        plot.tag.position = c(0.05, 1))

# Viewing the plot
print(p_tar_row)

# Saving the plot
ggsave("pred_tar.png", p_tar_row, width = 12, height = 6, dpi = 600)

# 3. Plots predicting the effects of light levels on each sonotype
p_lux_ACSH   <- pred_plot_all(ACSH_complete_Nb2, "Lux_cat_1", 
                              "Light Level", "ACSH", sig = TRUE, bias_correction = TRUE) + 
  coord_cartesian(ylim = c(0, 80)) + scale_y_continuous(breaks = seq(0, 80, 20))

p_lux_ALTE   <- pred_plot_all(ALTE_complete_Nb2_fixed, "Lux_cat_1", 
                              "Light Level", "ALTE", sig = FALSE) + 
  coord_cartesian(ylim = c(0, 80)) + scale_y_continuous(breaks = seq(0, 80, 20))

p_lux_TAPH   <- pred_plot_all(TAPH_complete_Nb2, "Lux_cat_1", 
                              "Light Level", "TAPH", sig = FALSE, bias_correction = TRUE) + 
  coord_cartesian(ylim = c(0, 80)) + scale_y_continuous(breaks = seq(0, 80, 20))

p_lux_NAMP   <- pred_plot_all(NAMP_complete_Nb2_fixed, "Lux_cat_1", "Light Level", "NAMP", sig = FALSE) + 
  coord_cartesian(ylim = c(0, 80)) + scale_y_continuous(breaks = seq(0, 80, 20))

# Creating a row of 4 plots representing each sonotype
p_lux_row <- wrap_plots(p_lux_ACSH, p_lux_ALTE, p_lux_TAPH, p_lux_NAMP, nrow = 1) +
  plot_layout(axes = "collect") +
  plot_annotation(tag_levels = list(c("I", "J", "K", "L"))) &
  theme(plot.tag = element_text(size = 10, face = "plain", hjust = 0.5, vjust = 0), 
        plot.tag.position = c(0.05, 1))

# Viewing the plot
print(p_lux_row)

# Saving the plot
ggsave("pred_lux.png", p_lux_row, width = 12, height = 6, dpi = 600)

# 4. Plots predicting the effects of presence of built-up on each sonotype
p_built_ACSH <- pred_plot_all(ACSH_complete_Nb2, "Builtup_cat", 
                              "Built-up", "ACSH", sig = TRUE, bias_correction = TRUE) + 
  coord_cartesian(ylim = c(0, 80)) + scale_y_continuous(breaks = seq(0, 80, 20))

p_built_ALTE <- pred_plot_all(ALTE_complete_Nb2_fixed, "Builtup_cat",
                              "Built-up", "ALTE", sig = FALSE) + coord_cartesian(ylim = c(0, 80)) + 
  scale_y_continuous(breaks = seq(0, 80, 20))

p_built_TAPH <- pred_plot_all(TAPH_complete_Nb2, "Builtup_cat", 
                              "Built-up", "TAPH", sig = FALSE, bias_correction = TRUE) + 
  coord_cartesian(ylim = c(0, 80)) + scale_y_continuous(breaks = seq(0, 80, 20))

p_built_NAMP <- pred_plot_all(NAMP_complete_Nb2_fixed, "Builtup_cat", 
                              "Built-up", "NAMP", sig = FALSE) + 
  coord_cartesian(ylim = c(0, 80)) + scale_y_continuous(breaks = seq(0, 80, 20))

# Creating a row of 4 plots representing each sonotype
p_built_row <- wrap_plots(p_built_ACSH, p_built_ALTE, p_built_TAPH, p_built_NAMP, nrow = 1) +
  plot_layout(axes = "collect") +
  plot_annotation(tag_levels = list(c("M", "N", "O", "P"))) &
  theme(plot.tag = element_text(size = 10, face = "plain", hjust = 0.5, vjust = 0), 
        plot.tag.position = c(0.05, 1))

# Viewing the plot
print (p_built_row)

# Saving the plot
ggsave("pred_builtup.png", p_built_row, width = 12, height = 6, dpi = 600)


# 5. Plots predicting the effects of moon light percentage on each sonotype
p_moon_ACSH  <- pred_plot_all(ACSH_complete_Nb2, "Moon_percentage", 
                              "Moon %", "ACSH", sig = FALSE, bias_correction = TRUE) + 
  coord_cartesian(ylim = c(0, 80)) + scale_y_continuous(breaks = seq(0, 80, 20))

p_moon_ALTE  <- pred_plot_all(ALTE_complete_Nb2_fixed, "Moon_percentage", 
                              "Moon %", "ALTE", sig = FALSE) + coord_cartesian(ylim = c(0, 80)) + 
  scale_y_continuous(breaks = seq(0, 80, 20))

p_moon_TAPH  <- pred_plot_all(TAPH_complete_Nb2, "Moon_percentage", 
                              "Moon %", "TAPH", sig = FALSE, bias_correction = TRUE) + 
  coord_cartesian(ylim = c(0, 80)) + scale_y_continuous(breaks = seq(0, 80, 20))

p_moon_NAMP  <- pred_plot_all(NAMP_complete_Nb2_fixed, "Moon_percentage", 
                              "Moon %", "NAMP", sig = FALSE) + 
  coord_cartesian(ylim = c(0, 80)) + scale_y_continuous(breaks = seq(0, 80, 20))

# Creating a row of 4 plots representing each sonotype
p_moon_row <- wrap_plots(p_moon_ACSH, p_moon_ALTE, p_moon_TAPH, p_moon_NAMP, nrow = 1) +
  plot_layout(axes = "collect") +
  plot_annotation(tag_levels = list(c("Q", "R", "S", "T"))) &
  theme(plot.tag = element_text(size = 10, face = "plain", hjust = 0.5, vjust = 0), 
        plot.tag.position = c(0.05, 1))

# Viewing the plot
print(p_moon_row)

# Saving the plot
ggsave("pred_moon.png", p_moon_row, width = 12, height = 6, dpi = 600)


# 6. Plots predicting the effects of season (winter and non-winter) on each sonotype

p_season_ACSH <- pred_plot_all(ACSH_complete_Nb2, "Season", "Season", 
                               "ACSH", sig = FALSE, bias_correction = TRUE) + 
  coord_cartesian(ylim = c(0, 80)) + scale_y_continuous(breaks = seq(0, 80, 20))

p_season_ALTE <- pred_plot_all(ALTE_complete_Nb2_fixed, "Season", 
                               "Season", "ALTE", sig = FALSE) + 
  coord_cartesian(ylim = c(0, 80)) + scale_y_continuous(breaks = seq(0, 80, 20))

p_season_TAPH <- pred_plot_all(TAPH_complete_Nb2, "Season", 
                               "Season", "TAPH", sig = TRUE, bias_correction = TRUE) + 
  coord_cartesian(ylim = c(0, 80)) + scale_y_continuous(breaks = seq(0, 80, 20))

p_season_NAMP <- pred_plot_all(NAMP_complete_Nb2_fixed, "Season", 
                               "Season", "NAMP", sig = TRUE) + 
  coord_cartesian(ylim = c(0, 80)) + scale_y_continuous(breaks = seq(0, 80, 20))

# Creating a row of 4 plots representing each sonotype
p_season_row <- wrap_plots(p_season_ACSH, p_season_ALTE, p_season_TAPH, p_season_NAMP, nrow = 1) +
  plot_layout(axes = "collect") +
  plot_annotation(tag_levels = list(c("U", "V", "W", "X"))) &
  theme(plot.tag = element_text(size = 10, face = "plain", hjust = 0.5, vjust = 0), 
        plot.tag.position = c(0.05, 1))

# Viewing the plot
print(p_season_row)

# Saving the plot
ggsave("pred_season.png", p_season_row, width = 12, height = 6, dpi = 600)


# Extracting coefficients and plotting them
# Creating a function for this purpose
extract_coefs <- function(model, species_name) {
  broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE) %>%
    dplyr::filter(term != "(Intercept)") %>%
    dplyr::mutate(
      Species = species_name,
      term = dplyr::case_when(
        term == "Vegetation_percentage" ~ "Vegetation %",
        term == "Tar_road_percentage"   ~ "Tar Road %",
        term == "Lux_cat_1High"          ~ "Light Level",
        term == "Builtup_catbuilt-up"    ~ "Built-up",
        term == "Moon_percentage"        ~ "Moon %",
        term == "Seasonnon-winter"       ~ "Season",
        TRUE ~ term
      ),
      Significant = ifelse(p.value < 0.05, "Significant", "Non-significant"),
      term = factor(term, levels = rev(c("Vegetation %", "Tar Road %", "Built-up",
                                         "Light Level", "Moon %", "Season")))
    )
}

# Combining and organizing fixed effects across species
coef_all <- dplyr::bind_rows(
  extract_coefs(ACSH_complete_Nb2,        "ACSH"),
  extract_coefs(ALTE_complete_Nb2_fixed,  "ALTE"),
  extract_coefs(TAPH_complete_Nb2,        "TAPH"),
  extract_coefs(NAMP_complete_Nb2_fixed,  "NAMP")
)

coef_all$Species <- factor(coef_all$Species, levels = c("ACSH", "ALTE", "TAPH", "NAMP"))

# Generating a Coefficient Plot ---
coef_plot <- ggplot(coef_all, aes(x = estimate, y = term)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "red", linewidth = 0.5) +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high), width = 0.2, colour = "grey40") +
  geom_point(aes(fill = Significant), shape = 21, size = 2.5) +
  scale_fill_manual(values = c("Significant" = "black", "Non-significant" = "white")) +
  facet_wrap(~ Species, ncol = 2) +
  labs(x = "Log-linear Estimate", y = "Predictor Term", fill = "") +
  theme_bw(base_size = 10) +
  theme(
    strip.text       = element_text(face = "bold.italic"),
    panel.grid.minor = element_blank(),
    axis.text.y      = element_text(size = 9),
    legend.position  = "bottom"
  )

# Viewing the plot
print(coef_plot)

# Saving the plot
ggsave("coefficient_plot_complete.png", coef_plot, width = 10, height = 8, dpi = 600)


# Moving on to community level analysis from species specific analysis
# Here we are considering all six sonotypes
# Rebuilding complete wide-format species summary (incorporates rare species)
species_summary_allbats <- species_data %>%
  dplyr::group_by(Site, Season, Species) %>%
  dplyr::summarise(Total_Calls = dplyr::n(), .groups = "drop")

species_wide_allbats <- species_summary_allbats %>%
  tidyr::pivot_wider(names_from = Species, values_from = Total_Calls, values_fill = 0)

# Merging back with environmental landscape vectors
diversity_data <- species_wide_allbats %>%
  dplyr::left_join(
    call_summary_locations %>%
      dplyr::select(Site, Season, Vegetation_percentage, Tar_road_percentage,
                    Mean_lux, Built_up_percentage, Moon_percentage),
    by = c("Site", "Season")
  ) %>%
  dplyr::mutate(
    Lux_cat_1   = factor(ifelse(Mean_lux < 1, "Low", "High"), levels = c("Low", "High")),
    Builtup_cat = factor(ifelse(Built_up_percentage == 0, "zero built-up", "built-up"), 
                         levels = c("zero built-up", "built-up")),
    Season      = factor(Season, levels = c("winter", "non-winter"))
  )

# Separating absolute species data frame for vegan format
sp_matrix <- diversity_data %>%
  dplyr::select(HETI, ACSH, ALTE, RHINO, NAMP, TAPH) %>%
  as.data.frame()

rownames(sp_matrix) <- paste(diversity_data$Site, diversity_data$Season, sep = "_")

# Performing beta diversity calculations
diversity_data$Shannon  <- vegan::diversity(sp_matrix, index = "shannon")
diversity_data$Richness <- vegan::specnumber(sp_matrix)


# Fitting three candidate models
Shannon_fixed <- glmmTMB(Shannon ~ Vegetation_percentage + Tar_road_percentage +
                           Lux_cat_1 + Builtup_cat + Moon_percentage + Season,
                         family = gaussian(), data = diversity_data)

Shannon_re    <- glmmTMB(Shannon ~ Vegetation_percentage + Tar_road_percentage +
                           Lux_cat_1 + Builtup_cat + Moon_percentage + Season + (1|Site),
                         family = gaussian(), data = diversity_data)

null_Shannon  <- glmmTMB(Shannon ~ 1, family = gaussian(), data = diversity_data)

# Testing the Gaussian assumption of the Shannon diversity model
shannon_sim <- DHARMa::simulateResiduals(Shannon_re, n = 1000, plot = TRUE)

# Extracting comparative diagnostics
print(AIC(null_Shannon, Shannon_fixed, Shannon_re))

# Proceeding with model containing site as a random effect
summary(Shannon_re)

# Performance descriptive validation metrics
performance::r2(Shannon_re)

# Visualizing bat call diversity predictions
# Creating a function for plotting
plot_shannon_pred <- function(model, predictor, x_label) {
  is_cat <- predictor %in% c("Lux_cat_1", "Builtup_cat", "Season")
  pred   <- ggeffects::ggpredict(model, terms = predictor, bias_correction = TRUE)
  
  if(is_cat) {
    p <- ggplot(pred, aes(x = x, y = predicted)) +
      geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.15, colour = "grey40") +
      geom_point(size = 2.5)
  } else {
    p <- ggplot(pred, aes(x = x, y = predicted)) +
      geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "grey80", alpha = 0.6) +
      geom_line(linewidth = 0.8)
  }
  
  p +
    labs(x = x_label, y = "Bat call diversity index") +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      axis.title.x     = element_text(size = 14),
      axis.title.y     = element_text(size = 14),
      axis.text        = element_text(size = 12)
    )
}

# Building each component (predictor) of the plot
p1 <- plot_shannon_pred(Shannon_re, "Vegetation_percentage", "Vegetation %") + coord_cartesian(ylim = c(0, 1)) + scale_y_continuous(breaks = seq(0, 1, 0.2))
p2 <- plot_shannon_pred(Shannon_re, "Tar_road_percentage",   "Tar Road %")   + coord_cartesian(ylim = c(0, 1)) + scale_y_continuous(breaks = seq(0, 1, 0.2))
p3 <- plot_shannon_pred(Shannon_re, "Lux_cat_1",             "Light level")  + coord_cartesian(ylim = c(0, 1)) + scale_y_continuous(breaks = seq(0, 1, 0.2))
p4 <- plot_shannon_pred(Shannon_re, "Builtup_cat",           "Built-up")     + coord_cartesian(ylim = c(0, 1)) + scale_y_continuous(breaks = seq(0, 1, 0.2))
p5 <- plot_shannon_pred(Shannon_re, "Moon_percentage",       "Moon %")       + coord_cartesian(ylim = c(0, 1)) + scale_y_continuous(breaks = seq(0, 1, 0.2))
p6 <- plot_shannon_pred(Shannon_re, "Season",                "Season")       + coord_cartesian(ylim = c(0, 1)) + scale_y_continuous(breaks = seq(0, 1, 0.2))

# Assembling grids
# Creating lables as per Journal's guidelines
combined_diversity_plot <- wrap_plots(p1, p2, p3, p4, p5, p6, ncol = 3) + 
  plot_layout(axes = "collect_y") + 
  plot_annotation(tag_levels = list(c("A", "B", "C*", "D", "E", "F"))) & 
  theme(
    plot.tag          = element_text(size = 10, face = "plain"), 
    plot.tag.position = "top"
  )

# Viewing and saving the plot
print(combined_diversity_plot)
ggsave("combined_shannon_diversity.png", combined_diversity_plot, width = 12, height = 8, dpi = 600)



