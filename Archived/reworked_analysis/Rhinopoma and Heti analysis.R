# ==============================================================================
# DIAGNOSTIC: DETECTION PROFILE FOR RHINO & HETI (87 SURVEY VISITS)
# ==============================================================================

# 1. Ensure we have the full 87-visit background grid
visit_grid <- Sampled_locations %>% 
  dplyr::select(Site, Season) %>% 
  dplyr::distinct()

# 2. Extract raw presence/counts from your bat data for these two species
rare_species_counts <- Bat_data %>%
  dplyr::filter(Species %in% c("RHINO", "HETI")) %>%
  dplyr::group_by(Site, Season, Species) %>%
  dplyr::summarise(Total_Calls = dplyr::n(), .groups = "drop")

# 3. Merge onto the 87 visits so we can explicitly see the zeros (absences)
rare_species_matrix <- visit_grid %>%
  dplyr::cross_join(data.frame(Species = c("RHINO", "HETI"))) %>%
  dplyr::left_join(rare_species_counts, by = c("Site", "Season", "Species")) %>%
  dplyr::mutate(
    Total_Calls = tidyr::replace_na(Total_Calls, 0),
    Presence = ifelse(Total_Calls > 0, 1, 0)
  )

# --- REPORT 1: Overall Summary Stats ---
cat("\n==================================================\n")
cat(" OVERALL DETECTION SUMMARY (Out of 87 Total Visits)\n")
cat("==================================================\n")
summary_stats <- rare_species_matrix %>%
  dplyr::group_by(Species) %>%
  dplyr::summarise(
    Visits_Present = sum(Presence),
    Visits_Absent  = sum(Presence == 0),
    Percent_Present = round((sum(Presence) / 87) * 100, 1),
    Total_Passes_Recorded = sum(Total_Calls)
  )
print(summary_stats)

# --- REPORT 2: Seasonal Breakdown ---
cat("\n==================================================\n")
cat(" SEASONAL DETECTION BREAKDOWN\n")
cat("==================================================\n")
seasonal_stats <- rare_species_matrix %>%
  dplyr::group_by(Species, Season) %>%
  dplyr::summarise(
    N_Visits_In_Season = dplyr::n(),
    N_Detections       = sum(Presence),
    Detection_Rate_Pct = round((sum(Presence) / dplyr::n()) * 100, 1),
    Total_Calls        = sum(Total_Calls),
    .groups = "drop"
  )
print(seasonal_stats)

# --- REPORT 3: Specific Hotspots (Where they actually show up) ---
cat("\n==================================================\n")
cat(" DETECTION HOTSPOTS (Sites with $>0$ passes)\n")
cat("==================================================\n")
hotspots <- rare_species_matrix %>%
  dplyr::filter(Presence == 1) %>%
  dplyr::select(Species, Site, Season, Total_Calls) %>%
  dplyr::arrange(Species, desc(Total_Calls))

if(nrow(hotspots) > 0) {
  print(as.data.frame(hotspots))
} else {
  cat("No detections found for either species in the raw dataset.\n")
}

# ==============================================================================
# HABITAT PROFILE FOR RARE SPECIES DETECTIONS
# ==============================================================================

hotspot_sites <- c("POPA4", "PERU3", "BERO1", "RUSH1", "RUBE2", "ITHI1", "ENMA4", "LACO1", "RADI2", "ENMA2", "BERO2", "CHBA1")

habitat_profiles <- Sampled_locations %>%
  dplyr::filter(Site %in% hotspot_sites) %>%
  dplyr::select(Site, Vegetation_percentage, Tar_road_percentage, Built_up_percentage, Mean_lux) %>%
  dplyr::distinct() %>%
  dplyr::arrange(desc(Vegetation_percentage))

cat("\n====================================================================\n")
cat(" ENVIRONMENT PROFILES FOR THE DETECTION SITES\n")
cat("====================================================================\n")
print(as.data.frame(habitat_profiles))