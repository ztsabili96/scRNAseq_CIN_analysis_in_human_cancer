# =====================================================
# Infer large-scale CNVs using copyKAT
# Windows-safe single-core execution
# =====================================================

library(Seurat)
library(copykat)

# ---------- Load Seurat object ----------
seu <- readRDS("results/seurat_object_qc.rds")

# ---------- Extract raw counts ----------
raw_counts <- GetAssayData(seu, slot = "counts")

# ---------- Run copyKAT ----------
ck_res <- copykat(
  rawmat = raw_counts,
  id.type = "S",          # Gene symbols
  ngene.chr = 5,
  win.size = 25,
  KS.cut = 0.1,
  sam.name = "Sample_01",
  distance = "euclidean",
  norm.cell.names = NULL,
  n.cores = 1             # set to 1 for stability and lower memory usage on Windows systems
)

# ---------- Save outputs ----------
dir.create("results/copykat", showWarnings = FALSE)

saveRDS(
  ck_res,
  file = "results/copykat/copykat_results.rds"
)

write.table(
  ck_res$prediction,
  file = "results/copykat/copykat_prediction.txt",
  sep = "\t",
  quote = FALSE
)

write.table(
  ck_res$CNAmat,
  file = "results/copykat/CNA_raw_gene_by_cell.txt",
  sep = "\t",
  quote = FALSE
)
