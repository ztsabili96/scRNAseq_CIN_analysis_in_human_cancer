# =====================================================
# Load raw scRNA-seq data and perform basic QC
# Dataset: Human breast cancer scRNA-seq
# =====================================================

library(Seurat)
library(Matrix)

# ---------- File paths (example, platform-independent) ----------
data_dir <- "data/raw"

counts <- readMM(file.path(data_dir, "counts.mtx"))
genes  <- read.delim(file.path(data_dir, "genes.tsv"), header = FALSE)
cells  <- read.delim(file.path(data_dir, "barcodes.tsv"), header = FALSE)

rownames(counts) <- genes$V1
colnames(counts) <- cells$V1

# ---------- Create Seurat object ----------
seu <- CreateSeuratObject(
  counts = counts,
  project = "Human_scRNA_CNV",
  min.cells = 3,
  min.features = 200
)

# ---------- Quality control metrics ----------
seu[["percent.mt"]] <- PercentageFeatureSet(
  seu,
  pattern = "^MT-"
)

# ---------- Filtering (CIN-friendly thresholds) ----------
seu <- subset(
  seu,
  subset =
    nFeature_RNA > 300 &
    nFeature_RNA < 7000 &
    percent.mt < 20
)

# ---------- Save QC-filtered object ----------
saveRDS(seu, file = "results/seurat_object_qc.rds")
