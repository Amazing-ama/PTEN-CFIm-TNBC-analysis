# PTEN and CFIm complex expression in basal-like breast cancer (TCGA-BRCA)
# Genes: PTEN, NUDT21 (CFIm25), CPSF6 (CFIm68), CPSF7 (CFIm59)

# ---- 1. Load packages ----
library(TCGAbiolinks)
library(SummarizedExperiment)
library(tidyr)
library(ggplot2)

# ---- 2. Download TCGA-BRCA RNA-seq data ----
query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
GDCdownload(query)
data <- GDCprepare(query)

# ---- 3. Select the four genes (TPM, log2-transformed) ----
genes <- c("PTEN", "NUDT21", "CPSF6", "CPSF7")

rows <- which(rowData(data)$gene_name %in% genes)
expr <- assay(data, "tpm_unstrand")[rows, ]
rownames(expr) <- rowData(data)$gene_name[rows]
expr <- log2(expr + 1)

# ---- 4. Label samples: Normal vs Basal-like tumor (PAM50) ----
subtypes <- TCGAquery_subtype("BRCA")

pam50 <- subtypes$BRCA_Subtype_PAM50[
  match(colData(data)$patient, subtypes$patient)
]

group <- ifelse(
  colData(data)$sample_type == "Solid Tissue Normal", "Normal",
  ifelse(pam50 == "Basal", "Basal-like tumor", NA)
)
print(table(group, useNA = "ifany"))

# ---- 5. Build one table, drop other subtypes ----
df <- as.data.frame(t(expr))
df$group <- group
df <- df[!is.na(df$group), ]

long <- pivot_longer(df, cols = all_of(genes),
                     names_to = "gene", values_to = "expression")

# ---- 6. Boxplot: tumor vs normal ----
p <- ggplot(long, aes(x = group, y = expression, fill = group)) +
  geom_boxplot() +
  facet_wrap(~ gene, scales = "free_y") +
  labs(y = "log2(TPM + 1)", x = NULL) +
  theme_minimal() +
  theme(legend.position = "none")

print(p)
ggsave("tailor_team_boxplot.png", plot = p, width = 8, height = 6, dpi = 300)

# ---- 7. Wilcoxon tests with multiple-testing correction ----
res <- do.call(rbind, lapply(genes, function(g) {
  data.frame(
    gene = g,
    median_tumor  = median(df[[g]][df$group == "Basal-like tumor"]),
    median_normal = median(df[[g]][df$group == "Normal"]),
    p = wilcox.test(df[[g]] ~ df$group)$p.value
  )
}))
res$p_adj <- p.adjust(res$p, method = "BH")
print(res)
write.csv(res, "wilcoxon_results.csv", row.names = FALSE)

# ---- 8. Correlation inside basal-like tumors ----
tumor <- df[df$group == "Basal-like tumor", ]
print(cor.test(tumor$CPSF6, tumor$PTEN, method = "spearman"))
