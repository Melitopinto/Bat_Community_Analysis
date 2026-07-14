library(tidyverse)
library(EcoCountHelper)
library(ggeffects)
library(patchwork)
library(broom.mixed)
library(glmmTMB)
library(performance)
library(vegan)

# ==============================================================================
# SECTION 0: OUTPUT DIRECTORY SETUP
# ==============================================================================
output_dir <- "reworked_analysis"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}


# SECTION 1: DATA CLEANING & MASTER HARMONIZATION
# Clean up trailing spaces in Site names
Sampled_locations$Site      <- trimws(Sampled_locations$Site)
Bat_data$Site               <- trimws(Bat_data$Site)
call_summary_locations$Site <- trimws(call_summary_locations$Site)

# 1.1 Isolate unique static landscape metrics per unique geographic Site ID 
site_habitat_matrix <- Sampled_locations %>%
  dplyr::group_by(Site) %>%
  dplyr::summarise(
    Vegetation_percentage   = dplyr::first(Vegetation_percentage),
    Built_up_percentage     = dplyr::first(Built_up_percentage),
    Open_barren_percentage  = dplyr::first(Open_barren_percentage),
    Tar_road_percentage     = dplyr::first(Tar_road_percentage),
    .groups = "drop"
  )

# 1.2 Replace landscape metrics in call summary with verified metrics via Site ID
call_summary_locations <- call_summary_locations %>%
  dplyr::select(-any_of(c("Vegetation_percentage", "Built_up_percentage", 
                          "Open_barren_percentage", "Tar_road_percentage"))) %>% 
  dplyr::left_join(site_habitat_matrix, by = "Site")

# 1.3 Fix factor columns and resolve taxonomic splits in raw pass data
species_data <- Bat_data %>%
  dplyr::mutate(
    Season = factor(Season, levels = c("winter", "non-winter")),
    Species = dplyr::case_when(
      Species == "PICE" ~ "PCSH", # Pipistrellus ceylonicus / Scotophilus heathii
      Species == "TAAE" ~ "NAMP", # Nyctinomus aegyptiacus / Mops plicatus
      Species == "MOPL" ~ "NAMP",
      Species == "TAME" ~ "TAPH", # Taphozous melanopogon / Taphozous nudiventris
      Species == "TANU" ~ "TAPH",
      TRUE ~ Species              # Retains PITE as-is
    )
  )

# 1.4 Aggregate total counts at the site-season-species level
species_summary <- species_data %>%
  dplyr::group_by(Site, Season, Species) %>%
  dplyr::summarise(Total_Calls = dplyr::n(), .groups = "drop")

# 1.5 Assemble full modeling grid using the complete 87-row visit matrix
Species_model_input <- Sampled_locations %>% 
  dplyr::select(Site, Season, Vegetation_percentage, Tar_road_percentage,
                Mean_lux, Built_up_percentage, Moon_percentage) %>% 
  dplyr::distinct() %>% 
  dplyr::cross_join(data.frame(Species = c("PITE", "PCSH", "TAPH", "NAMP"))) %>% 
  dplyr::left_join(
    species_summary %>% dplyr::select(Site, Season, Species, Total_Calls),
    by = c("Site", "Season", "Species")
  ) %>% 
  dplyr::mutate(
    Total_Calls = tidyr::replace_na(Total_Calls, 0),
    Season = factor(Season, levels = c("winter", "non-winter")),
    
    # Categorize urban metrics explicitly
    Builtup_cat = factor(ifelse(Built_up_percentage == 0, "zero built-up", "built-up"), 
                         levels = c("zero built-up", "built-up")),
    
    # Categorize light thresholds safely
    Lux_cat_1 = factor(ifelse(Mean_lux < 1, "Low", "High"), levels = c("Low", "High")),
    Lux_cat_5 = factor(ifelse(Mean_lux < 5, "Low", "High"), levels = c("Low", "High"))
  )

# Output summary tables to verify thresholds and data structure
cat("\n=== Lux threshold < 5 ===\n")
print(table(Species_model_input$Lux_cat_5))
cat("\n=== Lux threshold < 1 ===\n")
print(table(Species_model_input$Lux_cat_1))
cat("\n=== Builtup threshold 0 ===\n")
print(table(Species_model_input$Builtup_cat))

# Show zero metrics breakdown per species
print(
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
)

# ==============================================================================
# SECTION 2: ZERO-INFLATION DIAGNOSTIC LOOPS
# ==============================================================================

selected_species <- c("PITE", "PCSH", "TAPH", "NAMP")

for(sp in selected_species) {
  cat("\n========================================\n")
  cat(" Zero-Inflation Test for:", sp, "\n")
  sp_data <- Species_model_input %>% dplyr::filter(Species == sp)
  
  m <- glmmTMB(Total_Calls ~ Vegetation_percentage + Tar_road_percentage +
                 Lux_cat_1 + Builtup_cat + Moon_percentage + Season,
               family = nbinom2(), data = sp_data)
  
  m_zi <- glmmTMB(Total_Calls ~ Vegetation_percentage + Tar_road_percentage +
                    Lux_cat_1 + Builtup_cat + Moon_percentage + Season,
                  ziformula = ~1, family = nbinom2(), data = sp_data)
  
  sim     <- DHARMa::simulateResiduals(m, n = 1000)
  zi_test <- DHARMa::testZeroInflation(sim)
  
  cat("DHARMa ZI p-value:      ", round(zi_test$p.value, 3), "\n")
  cat("AIC Standard NB2:       ", round(AIC(m), 2), "\n")
  cat("AIC Zero-Inflated NB2: ", round(AIC(m_zi), 2), "\n")
  cat("ZI Model Required?     ", ifelse(AIC(m_zi) < AIC(m) - 2, "YES", "NO"), "\n")
}

# ==============================================================================
# SECTION 3: CANDIDATE MODEL ESTIMATION & SELECTION
# ==============================================================================

# --- 3.1 PCSH Models ---
PCSH_data <- Species_model_input %>% dplyr::filter(Species == "PCSH")
PCSH_complete_Nb2 <- glmmTMB(Total_Calls ~ Vegetation_percentage + Tar_road_percentage + 
                               Lux_cat_1 + Builtup_cat + Moon_percentage + Season + (1|Site), 
                             family = nbinom2(), data = PCSH_data)
PCSH_complete_Nb1 <- glmmTMB(Total_Calls ~ Vegetation_percentage + Tar_road_percentage + 
                               Lux_cat_1 + Builtup_cat + Moon_percentage + Season + (1|Site), 
                             family = nbinom1(), data = PCSH_data)

# --- 3.2 PITE Models (Fixed-Effects only due to random boundary constraints) ---
PITE_data <- Species_model_input %>% dplyr::filter(Species == "PITE")
PITE_complete_Nb2_fixed <- glmmTMB(Total_Calls ~ Vegetation_percentage + Tar_road_percentage + 
                                     Lux_cat_1 + Builtup_cat + Moon_percentage + Season, 
                                   family = nbinom2(), data = PITE_data)
PITE_complete_Nb1_fixed <- glmmTMB(Total_Calls ~ Vegetation_percentage + Tar_road_percentage + 
                                     Lux_cat_1 + Builtup_cat + Moon_percentage + Season, 
                                   family = nbinom1(), data = PITE_data)

# --- 3.3 TAPH Models ---
TAPH_data <- Species_model_input %>% dplyr::filter(Species == "TAPH")
TAPH_complete_Nb2 <- glmmTMB(Total_Calls ~ Vegetation_percentage + Tar_road_percentage + 
                               Lux_cat_1 + Builtup_cat + Moon_percentage + Season + (1|Site), 
                             family = nbinom2(), data = TAPH_data)
TAPH_complete_Nb1 <- glmmTMB(Total_Calls ~ Vegetation_percentage + Tar_road_percentage + 
                               Lux_cat_1 + Builtup_cat + Moon_percentage + Season + (1|Site), 
                             family = nbinom1(), data = TAPH_data)

# --- 3.4 NAMP Models (Fixed-Effects only due to random boundary constraints) ---
NAMP_data <- Species_model_input %>% dplyr::filter(Species == "NAMP")
NAMP_complete_Nb2_fixed <- glmmTMB(Total_Calls ~ Vegetation_percentage + Tar_road_percentage + 
                                     Lux_cat_1 + Builtup_cat + Moon_percentage + Season, 
                                   family = nbinom2(), data = NAMP_data)
NAMP_complete_Nb1_fixed <- glmmTMB(Total_Calls ~ Vegetation_percentage + Tar_road_percentage + 
                                     Lux_cat_1 + Builtup_cat + Moon_percentage + Season, 
                                   family = nbinom1(), data = NAMP_data)

# Run Framework Comparisons
Groups <- c("PCSH_complete", "PITE_complete", "TAPH_complete", "NAMP_complete")
ModelCompare(Groups, TopModOutName = "CompleteBestMods")

# Save and summarize final model definitions
CompleteBestMods <- list(
  PCSH = PCSH_complete_Nb2,
  PITE = PITE_complete_Nb2_fixed,
  TAPH = TAPH_complete_Nb2,
  NAMP = NAMP_complete_Nb2_fixed
)

summary(PCSH_complete_Nb2)
summary(PITE_complete_Nb2_fixed)
summary(TAPH_complete_Nb2)
summary(NAMP_complete_Nb2_fixed)

# ==============================================================================
# SECTION 4: DHARMA RESIDUAL EVALUATION
# ==============================================================================

try(dev.off(), silent = TRUE) 
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

sim_PCSH <- DHARMa::simulateResiduals(PCSH_complete_Nb2,        n = 1000, plot = TRUE)
sim_PITE <- DHARMa::simulateResiduals(PITE_complete_Nb2_fixed,  n = 1000, plot = TRUE)
sim_TAPH <- DHARMa::simulateResiduals(TAPH_complete_Nb2,        n = 1000, plot = TRUE)
sim_NAMP <- DHARMa::simulateResiduals(NAMP_complete_Nb2_fixed,  n = 1000, plot = TRUE)

# ==============================================================================
# SECTION 5: STATISTICAL METRIC EXTRACTION (IRR & R2 TABLES)
# ==============================================================================

# Extract incidence ratio rates (IRR) helper function
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

irr_complete <- dplyr::bind_rows(
  extract_irr(PCSH_complete_Nb2,        "PCSH"),
  extract_irr(PITE_complete_Nb2_fixed,  "PITE"),
  extract_irr(TAPH_complete_Nb2,        "TAPH"),
  extract_irr(NAMP_complete_Nb2_fixed,  "NAMP")
)
print(irr_complete)
write.csv(irr_complete, file.path(output_dir, "IRR_complete_models.csv"), row.names = FALSE)

# Generate pseudo-R2 values
r2_PCSH <- performance::r2_nakagawa(PCSH_complete_Nb2)
r2_TAPH <- performance::r2_nakagawa(TAPH_complete_Nb2)

null_PITE <- glmmTMB(Total_Calls ~ 1, family = nbinom2(), data = PITE_data)
null_NAMP <- glmmTMB(Total_Calls ~ 1, family = nbinom2(), data = NAMP_data)

r2_summary <- data.frame(
  Species = c("PCSH", "PITE", "TAPH", "NAMP"),
  R2_Marginal = c(
    round(r2_PCSH$R2_marginal, 3),
    round(1 - as.numeric(logLik(PITE_complete_Nb2_fixed)) / as.numeric(logLik(null_PITE)), 3),
    round(r2_TAPH$R2_marginal, 3),
    round(1 - as.numeric(logLik(NAMP_complete_Nb2_fixed)) / as.numeric(logLik(null_NAMP)), 3)
  ),
  R2_Conditional = c(
    round(r2_PCSH$R2_conditional, 3),
    NA, 
    round(r2_TAPH$R2_conditional, 3),
    NA  
  )
)
print(r2_summary)
write.csv(r2_summary, file.path(output_dir, "R2_complete_models.csv"), row.names = FALSE)

# ==============================================================================
# SECTION 6: SPECIES MARGINAL RESPONSE PLOTS
# ==============================================================================

pred_plot_all <- function(model, predictor, x_label, species_name, 
                          sig = FALSE, bias_correction = FALSE) {
  is_cat <- predictor %in% c("Lux_cat_1", "Builtup_cat", "Season")
  pred   <- ggeffects::ggpredict(model, terms = predictor, bias_correction = bias_correction)
  display_title <- if(sig) paste0(species_name, "*") else species_name
  
  if(is_cat) {
    p <- ggplot(pred, aes(x = x, y = predicted)) +
      geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.15, colour = "grey40") +
      geom_point(size = 2.5)
  } else {
    p <- ggplot(pred, aes(x = x, y = predicted)) +
      geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "grey80", alpha = 0.6) +
      geom_line(linewidth = 0.8)
  }
  p + labs(title = display_title, x = x_label, y = "Predicted Total Calls") +
    theme_bw() +
    theme(
      plot.title       = element_text(face = "plain", hjust = 0.5, size = 12),
      panel.grid.minor = element_blank(),
      axis.title.x     = element_text(size = 14),
      axis.title.y     = element_text(size = 14),
      axis.text        = element_text(size = 12)
    )
}

# Row 1: Vegetation
p_veg_PCSH   <- pred_plot_all(PCSH_complete_Nb2, "Vegetation_percentage", "Vegetation %", "PCSH", sig = FALSE, bias_correction = TRUE) + coord_cartesian(ylim = c(0, 125)) + scale_y_continuous(breaks = seq(0, 125, 25))
p_veg_PITE   <- pred_plot_all(PITE_complete_Nb2_fixed, "Vegetation_percentage", "Vegetation %", "PITE", sig = FALSE) + coord_cartesian(ylim = c(0, 125)) + scale_y_continuous(breaks = seq(0, 125, 25))
p_veg_TAPH   <- pred_plot_all(TAPH_complete_Nb2, "Vegetation_percentage", "Vegetation %", "TAPH", sig = FALSE, bias_correction = TRUE) + coord_cartesian(ylim = c(0, 125)) + scale_y_continuous(breaks = seq(0, 125, 25))
p_veg_NAMP   <- pred_plot_all(NAMP_complete_Nb2_fixed, "Vegetation_percentage", "Vegetation %", "NAMP", sig = FALSE) + coord_cartesian(ylim = c(0, 125)) + scale_y_continuous(breaks = seq(0, 125, 25))

p_veg_row <- wrap_plots(p_veg_PCSH, p_veg_PITE, p_veg_TAPH, p_veg_NAMP, nrow = 1) +
  plot_layout(axes = "collect") + plot_annotation(tag_levels = list(c("A", "B", "C", "D"))) &
  theme(plot.tag = element_text(size = 10, face = "plain", hjust = 0.5), plot.tag.position = c(0.05, 1))
ggsave(file.path(output_dir, "pred_veg.png"), p_veg_row, width = 12, height = 6, dpi = 600)

# Row 2: Tar Road
p_tar_PCSH   <- pred_plot_all(PCSH_complete_Nb2, "Tar_road_percentage", "Tar Road %", "PCSH", sig = FALSE, bias_correction = TRUE) + coord_cartesian(ylim = c(0, 100)) + scale_y_continuous(breaks = seq(0, 100, 25))
p_tar_PITE   <- pred_plot_all(PITE_complete_Nb2_fixed, "Tar_road_percentage", "Tar Road %", "PITE", sig = TRUE) + coord_cartesian(ylim = c(0, 100)) + scale_y_continuous(breaks = seq(0, 100, 25))
p_tar_TAPH   <- pred_plot_all(TAPH_complete_Nb2, "Tar_road_percentage", "Tar Road %", "TAPH", sig = FALSE, bias_correction = TRUE) + coord_cartesian(ylim = c(0, 100)) + scale_y_continuous(breaks = seq(0, 100, 25))
p_tar_NAMP   <- pred_plot_all(NAMP_complete_Nb2_fixed, "Tar_road_percentage", "Tar Road %", "NAMP", sig = FALSE) + coord_cartesian(ylim = c(0, 100)) + scale_y_continuous(breaks = seq(0, 100, 25))

p_tar_row <- wrap_plots(p_tar_PCSH, p_tar_PITE, p_tar_TAPH, p_tar_NAMP, nrow = 1) +
  plot_layout(axes = "collect") + plot_annotation(tag_levels = list(c("E", "F", "G", "H"))) &
  theme(plot.tag = element_text(size = 10, face = "plain", hjust = 0.5), plot.tag.position = c(0.05, 1))
ggsave(file.path(output_dir, "pred_tar.png"), p_tar_row, width = 12, height = 6, dpi = 600)

# Row 3: Light Levels
p_lux_PCSH   <- pred_plot_all(PCSH_complete_Nb2, "Lux_cat_1", "Light Level", "PCSH", sig = TRUE, bias_correction = TRUE) + coord_cartesian(ylim = c(0, 80)) + scale_y_continuous(breaks = seq(0, 80, 20))
p_lux_PITE   <- pred_plot_all(PITE_complete_Nb2_fixed, "Lux_cat_1", "Light Level", "PITE", sig = FALSE) + coord_cartesian(ylim = c(0, 80)) + scale_y_continuous(breaks = seq(0, 80, 20))
p_lux_TAPH   <- pred_plot_all(TAPH_complete_Nb2, "Lux_cat_1", "Light Level", "TAPH", sig = FALSE, bias_correction = TRUE) + coord_cartesian(ylim = c(0, 80)) + scale_y_continuous(breaks = seq(0, 80, 20))
p_lux_NAMP   <- pred_plot_all(NAMP_complete_Nb2_fixed, "Lux_cat_1", "Light Level", "NAMP", sig = FALSE) + coord_cartesian(ylim = c(0, 80)) + scale_y_continuous(breaks = seq(0, 80, 20))

p_lux_row <- wrap_plots(p_lux_PCSH, p_lux_PITE, p_lux_TAPH, p_lux_NAMP, nrow = 1) +
  plot_layout(axes = "collect") + plot_annotation(tag_levels = list(c("I", "J", "K", "L"))) &
  theme(plot.tag = element_text(size = 10, face = "plain", hjust = 0.5), plot.tag.position = c(0.05, 1))
ggsave(file.path(output_dir, "pred_lux.png"), p_lux_row, width = 12, height = 6, dpi = 600)

# Row 4: Built-up
p_built_PCSH <- pred_plot_all(PCSH_complete_Nb2, "Builtup_cat", "Built-up", "PCSH", sig = TRUE, bias_correction = TRUE) + coord_cartesian(ylim = c(0, 80)) + scale_y_continuous(breaks = seq(0, 80, 20))
p_built_PITE <- pred_plot_all(PITE_complete_Nb2_fixed, "Builtup_cat", "Built-up", "PITE", sig = FALSE) + coord_cartesian(ylim = c(0, 80)) + scale_y_continuous(breaks = seq(0, 80, 20))
p_built_TAPH <- pred_plot_all(TAPH_complete_Nb2, "Builtup_cat", "Built-up", "TAPH", sig = FALSE, bias_correction = TRUE) + coord_cartesian(ylim = c(0, 80)) + scale_y_continuous(breaks = seq(0, 80, 20))
p_built_NAMP <- pred_plot_all(NAMP_complete_Nb2_fixed, "Builtup_cat", "Built-up", "NAMP", sig = FALSE) + coord_cartesian(ylim = c(0, 80)) + scale_y_continuous(breaks = seq(0, 80, 20))

p_built_row <- wrap_plots(p_built_PCSH, p_built_PITE, p_built_TAPH, p_built_NAMP, nrow = 1) +
  plot_layout(axes = "collect") + plot_annotation(tag_levels = list(c("M", "N", "O", "P"))) &
  theme(plot.tag = element_text(size = 10, face = "plain", hjust = 0.5), plot.tag.position = c(0.05, 1))
ggsave(file.path(output_dir, "pred_builtup.png"), p_built_row, width = 12, height = 6, dpi = 600)

# Row 5: Moon %
p_moon_PCSH  <- pred_plot_all(PCSH_complete_Nb2, "Moon_percentage", "Moon %", "PCSH", sig = FALSE, bias_correction = TRUE) + coord_cartesian(ylim = c(0, 80)) + scale_y_continuous(breaks = seq(0, 80, 20))
p_moon_PITE  <- pred_plot_all(PITE_complete_Nb2_fixed, "Moon_percentage", "Moon %", "PITE", sig = FALSE) + coord_cartesian(ylim = c(0, 80)) + scale_y_continuous(breaks = seq(0, 80, 20))
p_moon_TAPH  <- pred_plot_all(TAPH_complete_Nb2, "Moon_percentage", "Moon %", "TAPH", sig = FALSE, bias_correction = TRUE) + coord_cartesian(ylim = c(0, 80)) + scale_y_continuous(breaks = seq(0, 80, 20))
p_moon_NAMP  <- pred_plot_all(NAMP_complete_Nb2_fixed, "Moon_percentage", "Moon %", "NAMP", sig = FALSE) + coord_cartesian(ylim = c(0, 80)) + scale_y_continuous(breaks = seq(0, 80, 20))

p_moon_row <- wrap_plots(p_moon_PCSH, p_moon_PITE, p_moon_TAPH, p_moon_NAMP, nrow = 1) +
  plot_layout(axes = "collect") + plot_annotation(tag_levels = list(c("Q", "R", "S", "T"))) &
  theme(plot.tag = element_text(size = 10, face = "plain", hjust = 0.5), plot.tag.position = c(0.05, 1))
ggsave(file.path(output_dir, "pred_moon.png"), p_moon_row, width = 12, height = 6, dpi = 600)

# Row 6: Season
p_season_PCSH <- pred_plot_all(PCSH_complete_Nb2, "Season", "Season", "PCSH", sig = FALSE, bias_correction = TRUE) + coord_cartesian(ylim = c(0, 80)) + scale_y_continuous(breaks = seq(0, 80, 20))
p_season_PITE <- pred_plot_all(PITE_complete_Nb2_fixed, "Season", "Season", "PITE", sig = FALSE) + coord_cartesian(ylim = c(0, 80)) + scale_y_continuous(breaks = seq(0, 80, 20))
p_season_TAPH <- pred_plot_all(TAPH_complete_Nb2, "Season", "Season", "TAPH", sig = TRUE, bias_correction = TRUE) + coord_cartesian(ylim = c(0, 80)) + scale_y_continuous(breaks = seq(0, 80, 20))
p_season_NAMP <- pred_plot_all(NAMP_complete_Nb2_fixed, "Season", "Season", "NAMP", sig = TRUE) + coord_cartesian(ylim = c(0, 80)) + scale_y_continuous(breaks = seq(0, 80, 20))

p_season_row <- wrap_plots(p_season_PCSH, p_season_PITE, p_season_TAPH, p_season_NAMP, nrow = 1) +
  plot_layout(axes = "collect") + plot_annotation(tag_levels = list(c("U", "V", "W", "X"))) &
  theme(plot.tag = element_text(size = 10, face = "plain", hjust = 0.5), plot.tag.position = c(0.05, 1))
ggsave(file.path(output_dir, "pred_season.png"), p_season_row, width = 12, height = 6, dpi = 600)

# ==============================================================================
# SECTION 7: FIXED EFFECTS COEFFICIENT GRID PLOT
# ==============================================================================

extract_coefs <- function(model, species_name) {
  broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE) %>%
    dplyr::filter(term != "(Intercept)") %>%
    dplyr::mutate(
      Species = species_name,
      term = dplyr::case_when(
        term == "Vegetation_percentage" ~ "Vegetation %",
        term == "Tar_road_percentage"   ~ "Tar Road %",
        term == "Lux_cat_1High"          ~ "Light Level",
        term == "Builtup_catbuilt-up"   ~ "Built-up",
        term == "Moon_percentage"       ~ "Moon %",
        term == "Seasonnon-winter"      ~ "Season",
        TRUE ~ term
      ),
      Significant = ifelse(p.value < 0.05, "Significant", "Non-significant"),
      term = factor(term, levels = rev(c("Vegetation %", "Tar Road %", "Built-up", "Light Level", "Moon %", "Season")))
    )
}

coef_all <- dplyr::bind_rows(
  extract_coefs(PCSH_complete_Nb2,        "PCSH"),
  extract_coefs(PITE_complete_Nb2_fixed,  "PITE"),
  extract_coefs(TAPH_complete_Nb2,        "TAPH"),
  extract_coefs(NAMP_complete_Nb2_fixed,  "NAMP")
)
coef_all$Species <- factor(coef_all$Species, levels = c("PCSH", "PITE", "TAPH", "NAMP"))

coef_plot <- ggplot(coef_all, aes(x = estimate, y = term)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "red", linewidth = 0.5) +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high), width = 0.2, colour = "grey40") +
  geom_point(aes(fill = Significant), shape = 21, size = 2.5) +
  scale_fill_manual(values = c("Significant" = "black", "Non-significant" = "white")) +
  facet_wrap(~ Species, ncol = 2) +
  labs(x = "Log-linear Estimate", y = "Predictor Term", fill = "") +
  theme_bw(base_size = 10) +
  theme(strip.text = element_text(face = "bold.italic"), panel.grid.minor = element_blank(), legend.position = "bottom")

ggsave(file.path(output_dir, "coefficient_plot_complete.png"), coef_plot, width = 10, height = 8, dpi = 600)

# ==============================================================================
# SECTION 8: COMMUNITY DATA ASSEMBLY & BIODIVERSITY METRICS (ALL 6 SPECIES)
# ==============================================================================

species_summary_allbats <- species_data %>%
  dplyr::group_by(Site, Season, Species) %>%
  dplyr::summarise(Total_Calls = dplyr::n(), .groups = "drop")

species_wide_allbats <- species_summary_allbats %>%
  tidyr::pivot_wider(names_from = Species, values_from = Total_Calls, values_fill = 0)

# Merged safely against the master 87 survey visit log matrix
diversity_data <- species_wide_allbats %>%
  dplyr::left_join(
    Sampled_locations %>%
      dplyr::select(Site, Season, Vegetation_percentage, Tar_road_percentage,
                    Mean_lux, Built_up_percentage, Moon_percentage),
    by = c("Site", "Season")
  ) %>%
  dplyr::mutate(
    Lux_cat_1   = factor(ifelse(Mean_lux < 1, "Low", "High"), levels = c("Low", "High")),
    Builtup_cat = factor(ifelse(Built_up_percentage == 0, "zero built-up", "built-up"), levels = c("zero built-up", "built-up")),
    Season      = factor(Season, levels = c("winter", "non-winter"))
  )

# Separate species data matrix for strict community calculation
sp_matrix <- diversity_data %>%
  dplyr::select(any_of(c("HETI", "PCSH", "PITE", "RHINO", "NAMP", "TAPH"))) %>%
  dplyr::mutate(across(everything(), ~tidyr::replace_na(., 0))) %>% 
  as.data.frame()

rownames(sp_matrix) <- paste(diversity_data$Site, diversity_data$Season, sep = "_")

# Compute Ecological Community Indices
diversity_data$Shannon  <- vegan::diversity(sp_matrix, index = "shannon")
diversity_data$Richness <- vegan::specnumber(sp_matrix)

# ==============================================================================
# SECTION 9: ECO-COMMUNITY SHANNON DIVERSITY REGRESSION
# ==============================================================================

# Explicitly using Fixed effects setup for Shannon index
Shannon_fixed <- glmmTMB(Shannon ~ Vegetation_percentage + Tar_road_percentage +
                           Lux_cat_1 + Builtup_cat + Moon_percentage + Season,
                         family = gaussian(), data = diversity_data)

Shannon_re    <- glmmTMB(Shannon ~ Vegetation_percentage + Tar_road_percentage +
                           Lux_cat_1 + Builtup_cat + Moon_percentage + Season + (1|Site),
                         family = gaussian(), data = diversity_data)

null_Shannon  <- glmmTMB(Shannon ~ 1, family = gaussian(), data = diversity_data)

print(AIC(null_Shannon, Shannon_fixed, Shannon_re))
summary(Shannon_re)
print(performance::r2(Shannon_re))

# ==============================================================================
# SECTION 10: COMMUNITY VISUALIZATIONS
# ==============================================================================

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
  p + labs(x = x_label, y = "Bat call diversity index") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(), axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14), axis.text = element_text(size = 12))
}

p1 <- plot_shannon_pred(Shannon_re, "Vegetation_percentage", "Vegetation %") + coord_cartesian(ylim = c(0, 1)) + scale_y_continuous(breaks = seq(0, 1, 0.2))
p2 <- plot_shannon_pred(Shannon_re, "Tar_road_percentage",   "Tar Road %")   + coord_cartesian(ylim = c(0, 1)) + scale_y_continuous(breaks = seq(0, 1, 0.2))
p3 <- plot_shannon_pred(Shannon_re, "Lux_cat_1",             "Light level")  + coord_cartesian(ylim = c(0, 1)) + scale_y_continuous(breaks = seq(0, 1, 0.2))
p4 <- plot_shannon_pred(Shannon_re, "Builtup_cat",           "Built-up")     + coord_cartesian(ylim = c(0, 1)) + scale_y_continuous(breaks = seq(0, 1, 0.2))
p5 <- plot_shannon_pred(Shannon_re, "Moon_percentage",       "Moon %")       + coord_cartesian(ylim = c(0, 1)) + scale_y_continuous(breaks = seq(0, 1, 0.2))
p6 <- plot_shannon_pred(Shannon_re, "Season",                "Season")       + coord_cartesian(ylim = c(0, 1)) + scale_y_continuous(breaks = seq(0, 1, 0.2))

combined_diversity_plot <- wrap_plots(p1, p2, p3, p4, p5, p6, ncol = 3) + 
  plot_layout(axes = "collect_y") + 
  plot_annotation(tag_levels = list(c("A", "B", "C*", "D", "E", "F"))) & 
  theme(plot.tag = element_text(size = 10, face = "plain"), plot.tag.position = "top")

ggsave(file.path(output_dir, "combined_shannon_diversity.png"), combined_diversity_plot, width = 12, height = 8, dpi = 600)

# Save project objects explicitly inside the target directory
save(Shannon_re, CompleteBestMods, diversity_data, combined_diversity_plot,
     file = file.path(output_dir, "bat_diversity_results.RData"))


# 1. Create a structured dataframe containing your final model summary results
manuscript_table <- data.frame(
  Species = c(
    "PCSH (nbinom2)", rep("", 6),
    "PITE (nbinom2 - Fixed only)", rep("", 6),
    "TAPH (nbinom2)", rep("", 6),
    "NAMP (nbinom2 - Fixed only)", rep("", 6)
  ),
  Predictor_Variable = rep(c(
    "(Intercept)", "Vegetation %", "Tar Road %", 
    "Light (High vs Low)", "Urbanization (Built-up)", "Moon %", "Season (Non-winter)"
  ), 4),
  Estimate = c(
    3.4874, -0.0040, -0.0164, -0.9573, -1.0374, -0.0007,  0.4850,  # PCSH
    4.0725, -0.0052, -0.0196,  0.2465,  0.0998,  0.0024,  0.2657,  # PITE
    2.0683, -0.0070, -0.0130, -0.0361, -0.9612, -0.0049,  1.4668,  # TAPH
    1.0769, -0.0089, -0.0241, -1.2810, -1.0137, -0.0016,  2.4860   # NAMP
  ),
  Std_Error = c(
    0.5824, 0.0068, 0.0095, 0.3245, 0.4433, 0.0040, 0.2501,        # PCSH
    0.5212, 0.0067, 0.0091, 0.3263, 0.4005, 0.0042, 0.3486,        # PITE
    0.9696, 0.0115, 0.0160, 0.6154, 0.7595, 0.0079, 0.5043,        # TAPH
    1.2092, 0.0118, 0.0171, 0.7621, 0.7566, 0.0118, 0.7844         # NAMP
  ),
  z_value = c(
    5.988, -0.600, -1.735, -2.950, -2.340, -0.181,  1.939,         # PCSH
    7.814, -0.768, -2.155,  0.755,  0.249,  0.572,  0.762,         # PITE
    2.133, -0.611, -0.814, -0.059, -1.266, -0.621,  2.909,         # TAPH
    0.891, -0.765, -1.408, -1.681, -1.340, -0.140,  3.169          # NAMP
  ),
  P_Value = c(
    "< 0.001", "0.5488", "0.0827", "0.0032", "0.0193", "0.8566", "0.0525", # PCSH
    "< 0.001", "0.4426", "0.0311", "0.4500", "0.8032", "0.5673", "0.4459", # PITE
    "0.0329",  "0.5409", "0.4158", "0.9532", "0.2057", "0.5348", "0.0036", # TAPH
    "0.3731",  "0.4445", "0.1591", "0.0928", "0.1803", "0.8890", "0.0015"  # NAMP
  ),
  Significance = c(
    "***", "", ".", "**", "*", "", ".",                           # PCSH
    "***", "", "*", "", "", "", "",                               # PITE
    "*", "", "", "", "", "", "**",                                # TAPH
    "", "", "", ".", "", "", "**"                                 # NAMP
  )
)

# 2. Export to a clean CSV file in your working directory
write.csv(manuscript_table, "GLMM_Bat_Model_Results.csv", row.names = FALSE)

cat("Success! 'GLMM_Bat_Model_Results.csv' has been generated for your manuscript.\n")