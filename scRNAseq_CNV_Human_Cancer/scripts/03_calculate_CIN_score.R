# =====================================================
# Quantification of chromosomal instability (CIN score)
# =====================================================

library(dplyr)
library(ggplot2)

# ---------- Load copyKAT results ----------
cnv_raw <- read.delim(
  "results/copykat/CNA_raw_gene_by_cell.txt",
  row.names = 1
)

pred <- read.delim(
  "results/copykat/copykat_prediction.txt",
  stringsAsFactors = FALSE
)

# ---------- Align cell names ----------
common_cells <- intersect(colnames(cnv_raw), pred$cell.names)

cnv_raw <- cnv_raw[, common_cells]
pred <- pred[pred$cell.names %in% common_cells, ]

# ---------- CIN score definition ----------
# Mean absolute CNV deviation per cell
cin_score <- colMeans(abs(cnv_raw), na.rm = TRUE)

cin_df <- data.frame(
  cell_id = names(cin_score),
  CIN_score = cin_score
) %>%
  left_join(pred, by = c("cell_id" = "cell.names"))

# ---------- Save results ----------
write.csv(
  cin_df,
  file = "results/CIN_score_per_cell.csv",
  row.names = FALSE
)

# ---------- Visualization ----------
p <- ggplot(cin_df, aes(x = copykat.pred, y = CIN_score, fill = copykat.pred)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  theme_classic(base_size = 14) +
  labs(
    title = "Chromosomal Instability (CIN) per Cell",
    x = "Cell classification (copyKAT)",
    y = "Mean absolute CNV deviation"
  )

ggsave(
  filename = "figures/CIN_score_boxplot.png",
  plot = p,
  width = 6,
  height = 4
)

# Session information for reproducibility
# ------------------------------
sessionInfo()