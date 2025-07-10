# run this after you run 0-Figure-updated

# ============================
# 📦 Load Required Libraries
# ============================
library(ggplot2)
library(pls)
library(dplyr)
library(readxl)
library(ggtext)

# ============================
# 🗂 Define Paths
# ============================
dir_results <- "/rds/projects/x/xiap-xia-transcriptomics/AOP network/Results/"
tAOP_dir    <- file.path(dir_results, "tAOP")

# ============================
# 🔍 Extract MIE Terms from tAOP
# ============================
extract_tAOP_terms <- function(dir) {
  list.files(dir, pattern = "MIE", full.names = TRUE) %>%
    lapply(read.csv) %>%
    lapply(function(df) df[df[, 4] == "mie", 2]) %>%
    unlist() %>%
    unique()
}
tAOP_terms <- extract_tAOP_terms(tAOP_dir)

# ============================
# 🧬 Prepare Expression Data
# ============================
x_base <- t(hdata) * -1
colnames(x_base) <- path.db[[2]][match(colnames(x_base), path.db[[2]][, 2]), 1]
rownames(x_base) <- colnames(hdata)

anno <- data.nano[, c(2:4)]
rownames(anno) <- anno[, 3]

# ============================
# 🧪 Define PLS Model Configurations
# ============================
model_setups <- list(
  list(
    name   = "All KEs",
    color  = "#E41A1C",
    columns = colnames(x_base)
  ),
  list(
    name   = "Neuro CRKEs",
    color  = "#377EB8",
    columns = tAOP_terms
  ),
  list(
    name   = "Neuro CR+TRKEs",
    color  = "#4DAF4A",
    columns = {
      neuro_cols <- grep("NEURO|LEARN|AHR", colnames(x_base), value = TRUE)
      union(neuro_cols, tAOP_terms)
    }
  )
)

# ============================
# ⚗️ Run PLS Models and Collect Results
# ============================
all_models <- lapply(model_setups, function(cfg) {
  selected <- na.omit(match(cfg$columns, colnames(x_base)))
  x <- x_base[, selected, drop = FALSE]
  y <- anno[rownames(x), ]
  y <- y[y[,2]!='No',]
  y1 <- y
  y <- y[,1]
  names(y) <- rownames(y1)
  x <- x[rownames(x)%in%y1[,3],]
  x <- x[order(rownames(x)),]
  y <- y[order(names(y))]

  pls_mod <- pls(x, y, center = TRUE, scale = TRUE, cv = 1)
  ncomp   <- pls_mod$ncomp.selected
  
  df <- data.frame(
    Measured  = pls_mod$cvres$y.ref[, 1],
    Predicted = pls_mod$cvres$y.pred[, ncomp, 1],
    Compound  = rownames(pls_mod$cvres$y.pred[, , 1]),
    Model     = cfg$name,
    Color     = cfg$color
  )
  
  fit <- lm(Predicted ~ Measured, data = df)
  df$Slope <- round(coef(fit)[2], 2)
  df$R2    <- round(summary(fit)$r.squared, 2)
  df$Pval  <- formatC(summary(fit)$coefficients[2, 4], format = "e", digits = 2)
  df$RMSE  <- round(sqrt(mean((df$Predicted - df$Measured)^2)), 2)
  
  return(df)
})

merged_df <- bind_rows(all_models)
merged_df$Model <- factor(merged_df$Model, levels = c("All KEs", "Neuro CRKEs", "Neuro CR+TRKEs"))

# ============================
# 🎨 Format Legend Labels with HTML Styling
# ============================
legend_labels <- merged_df %>%
  group_by(Model, Color) %>%
  summarise(
    Slope = first(Slope),
    R2    = first(R2),
    Pval  = first(Pval),
    RMSE  = first(RMSE),
    .groups = "drop"
  ) %>%
  mutate(
    html_label = paste0(
      "<span style='color:", Color, "'><b>", Model, "</b></span><br>",
      "Slope = ", Slope, ", R² = ", R2, "<br>",
      "P = ", Pval, ", RMSE = ", RMSE
    )
  )

# Named color and label mappings
legend_colors <- setNames(legend_labels$Color, legend_labels$Model)
legend_texts  <- setNames(legend_labels$html_label, legend_labels$Model)

# ============================
# 📈 ggplot2 Visualization with ggtext Legend
# ============================
gg <- ggplot(merged_df, aes(x = Measured, y = Predicted, color = Model)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray60") +
  scale_color_manual(
    values = legend_colors,
    labels = legend_texts,
    guide  = guide_legend(override.aes = list(size = 4))
  ) +
  labs(
    title = "",
    x = "Measured PODa",
    y = "Predicted PODa",
    color = NULL
  ) +
  theme_minimal(base_size = 20) +
  theme(
    legend.position = "bottom",
    legend.text     = element_markdown(size = 24),  # Increased from 20 to 24
    plot.title      = element_text(face = "bold", hjust = 0.5),
    axis.title.x    = element_text(size = 25, face = "bold"),
    axis.title.y    = element_text(size = 25, face = "bold"),
    axis.line       = element_line(color = "black", linewidth = 1),
    axis.ticks      = element_line(linewidth = 1),
    axis.text       = element_text(size = 20)
  )
# ============================
# 💾 Save Annotated Plot
# ============================
ggsave(
  filename = file.path(dir_results, "Figure5_Merged_ggplot_coloredlegend.pdf"),
  plot     = gg,
  width    = 16,
  height   = 16
)
