library(ggplot2)
library(dplyr)
library(tidyr)

# Step 1: Load data
Sampled_locations <- read.csv("sampled_locations.csv")

# Step 2: Clean NAs, aggregate unique sites, and reshape
plot_data <- Sampled_locations %>%
  select(Site, 
         Vegetation_percentage,
         Built_up_percentage, 
         Open_barren_percentage, 
         Tar_road_percentage, 
         Unclassified_habitat_percentage) %>%
  
  # Replace any NA values in the percentage columns with 0
  mutate(across(ends_with("_percentage"), ~replace_na(., 0))) %>%
  
  # Group by Site and calculate the mean to get your unique sites
  group_by(Site) %>%
  summarise(across(everything(), mean)) %>%
  ungroup() %>%
  
  # Pivot to long format for ggplot
  pivot_longer(
    cols = -Site, 
    names_to = "Habitat", 
    values_to = "Percentage"
  )

habitat_comp_plot <- ggplot(plot_data, aes(x = Site, y = Percentage, fill = Habitat)) +
  geom_bar(stat = "identity", width = 0.85) +
  
  # Explicitly mapping the exact pivoted names to their colors and clean labels
  scale_fill_manual(
    values = c(
      "Built_up_percentage"             = "#F8766D", # Coral Red
      "Open_barren_percentage"          = "#7CAE00", # Olive Green
      "Tar_road_percentage"             = "#00BFC4", # Teal
      "Vegetation_percentage"           = "#C77CFF", # Purple
      "Unclassified_habitat_percentage" = "#999999"  # Grey
    ),
    breaks = c(
      "Built_up_percentage", 
      "Open_barren_percentage", 
      "Tar_road_percentage", 
      "Vegetation_percentage", 
      "Unclassified_habitat_percentage"
    ),
    labels = c("Built Up", "Open Barren", "Tar Road", "Vegetation", "Unclassified")
  ) +
  
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 25)) +
  theme_classic() +
  labs(x = "Sampled sites", y = "Habitat type percentage", fill = "Habitat") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "right"
  )

print(habitat_comp_plot)

# 2. Save the plot using ggeffects/ggplot's native saving engine
ggsave(
  filename = "habitat_composition_final.png", 
  plot     = habitat_comp_plot,
  width    = 16,       # Wide enough so 53 sites don't squish together
  height   = 8,        # Balanced height
  units    = "in",     # Dimensions in inches
  dpi      = 600       # High-resolution print/publication quality
)
