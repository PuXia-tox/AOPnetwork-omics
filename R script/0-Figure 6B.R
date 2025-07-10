# Required libraries
library(dplyr)
library(tidyr)
library(ggplot2)

# Updated dataset including species
data <- data.frame(
  Sample = c("TBT", "PFOS_2", "FBSA", "EE2", "PFOA", "PFOS_1", "ADD", "TCPP", 
             "DBP", "DiBP", "DGD", "DINP", "DINCH", "GTA"),
  Measured = c(-1.251811973, 2.300917094, 2.382926234, -2.241433674, 2.141308175,
               -0.108018299, -1.698970004, 2.086755871, 0.633468456, 0.740362689,
               1.089905111, 1.956648579, 1.912222057, 3.491361694),
  Estimated = c(-1.299988938, 1.424768735, 2.122972958, -2.285333072, -1.139952512,
                -0.699082906, -3.22184875, 1.279252934, 0.899820502, -1.698970004,
                -0.337242168, 0.930949031, 1.628491105, 1.334252642),
  Species = c("Zebrafish", "Zebrafish", "Zebrafish", "Rainbow Trout", "Zebrafish",
              "Zebrafish", "Zebrafish", "Yellowstripe goby", "Zebrafish", "Zebrafish",
              "Zebrafish", "Zebrafish", "Zebrafish", "Zebrafish")
)

# Create ±1 and ±2 log unit ranges
data <- data %>%
  arrange(Measured) %>%
  mutate(Sample = factor(Sample, levels = Sample),
         LowerBound1 = Measured - 1,
         UpperBound1 = Measured + 1,
         LowerBound2 = Measured - 2,
         UpperBound2 = Measured + 2)

# Long format for plotting points
long_data <- pivot_longer(data, cols = c("Measured", "Estimated"),
                          names_to = "Type", values_to = "PODa")

# Shape mapping for species
shape_map <- c("Zebrafish" = 16, "Rainbow Trout" = 15, "Yellowstripe goby" = 17)

# Journal-style plot with shape-coded species
ggplot() +
  # ±2 log unit zone (outer band)
  geom_rect(data = data, aes(xmin = as.numeric(Sample) - 0.3,
                             xmax = as.numeric(Sample) + 0.3,
                             ymin = LowerBound2, ymax = UpperBound2),
            fill = "#FDEACF", alpha = 0.4) +
  # ±1 log unit zone (inner band)
  geom_rect(data = data, aes(xmin = as.numeric(Sample) - 0.3,
                             xmax = as.numeric(Sample) + 0.3,
                             ymin = LowerBound1, ymax = UpperBound1),
            fill = "#DCE6F2", alpha = 0.6) +
  # Measured & Estimated points with shape based on species
  geom_point(data = long_data, aes(x = Sample, y = PODa, color = Type, shape = Species),
             size = 3, position = position_dodge(width = 0.4)) +
  scale_color_manual(values = c("Measured" = "#2166ac", "Estimated" = "#b2182b")) +
  scale_shape_manual(values = shape_map) +
  labs(x = "Compounds", y = "PODa Values", color = "Type", shape = "Species") +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.grid.major = element_line(color = "gray90", linetype = "dashed"),
    panel.grid.minor = element_line(color = "gray90", linetype = "dashed"),
    axis.line.x = element_line(color = "black", linewidth = 0.6),
    axis.line.y = element_line(color = "black", linewidth = 0.6),
    panel.border = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank()
  )