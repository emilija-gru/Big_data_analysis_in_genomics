library(dplyr)
library(ggplot2)
library(tidyr)
library(pheatmap)
library(DESeq2)
library(RColorBrewer)
library(clusterProfiler)
library(org.Hs.eg.db)
library(reshape2)
library(msigdbr)
library(enrichplot)
library(aPEAR)
library(ggrepel)
library(arules)

# nustatau, kad HW1 būtų darbinė direktorija
getwd()
setwd("~/HW1")

################ mapping rate for each sample plot ################

samples <- c("SRR11647692", "SRR11647693", "SRR11647694",
             "SRR11647702", "SRR11647703", "SRR11647704")

mapping_rate <- c(95.2, 94.8, 95.5, 94.9, 95.1, 94.6)

png("mapping_rate_plot.png", width = 1000, height = 700)

par(mar = c(10, 5, 4, 2))

barplot(mapping_rate,
        names.arg = samples,
        las = 2,
        col = "steelblue",
        ylab = "Mapping rate (%)",
        ylim = c(0, 100),
        main = "Read mapping rate per sample")

dev.off()

################ feature assignment rates (from featureCounts) ################

summary <- read.delim("counts.txt.summary")

summary

# pakeičių duomenis į log formatą

summary_long <- summary |>
  pivot_longer(
    cols = -Status,
    names_to = "sample",
    values_to = "reads"
  )

head(summary_long)

# apskaičiuoju procentus

summary_long <- summary_long |>
  group_by(sample) |>
  mutate(percent = reads / sum(reads) * 100)

# nubraižom bar plot

### rodo tik assigned ir unassigned
summary_simple <- summary_long |>
  mutate(Status = ifelse(Status == "Assigned", "Assigned", "Unassigned")) |>
  group_by(sample, Status) |>
  summarise(percent = sum(percent))

p <- ggplot(summary_simple, aes(sample, percent, fill = Status)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Assigned vs unassigned reads",
    y = "Percentage of reads"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("feature_assignment_rates.png", plot = p, width = 8, height = 6)

### rodo visas grupes
# p <- ggplot(summary_long, aes(sample, percent, fill = Status)) +
# geom_col() +
# labs(
#    title = "Feature assignment rates",
#    x = "Sample",
#    y = "Percentage of reads"
#  ) +
#  theme_minimal()

# ggsave("feature_assignment_rates.png", plot = p, width = 8, height = 6)

# ((((((((((((((((((()))))))))))))))))))
# {{{{{{{{{{{{{{{{{{{}}}}}}}}}}}}}}}}}}}
# DATA STATISTICAL ANALYSIS
# ((((((((((((((((((()))))))))))))))))))
# {{{{{{{{{{{{{{{{{{{}}}}}}}}}}}}}}}}}}}

# įsikeliu counts failą
counts_raw <- read.table("counts.txt",
                         header = TRUE,
                         sep = "\t",
                         comment.char = "#",
                         check.names = FALSE)

counts_raw$Geneid <- gsub("\\..*$", "", counts_raw$Geneid)

# Patikrinu, ar viskas tikrai gerai įsikėlė
head(counts_raw[, 1:7])
dim(counts_raw)
colnames(counts_raw)

# pasidarau count_matrix

# paimu tik read count stulpelius ir genų ID padarau eilutėmis
count_matrix <- counts_raw[, 7:ncol(counts_raw)]

rownames(count_matrix) <- counts_raw$Geneid

dim(count_matrix)
head(count_matrix[, 1:3])

# Įsikeliu sampleInfo.txt

meta <- read.table(
  "sampleInfo.txt",
  header = TRUE,
  sep = " "
)

head(meta)
colnames(count_matrix)

# Patikrinu, ar meta sutampa su count_matrix

all(colnames(count_matrix) == meta$sample)
colnames(count_matrix) <- gsub("_sorted.bam$", "", colnames(count_matrix))

# 11111111

# LIBRARY SIZE EXPLORATION
colSums(count_matrix) # čia reikės total counts grafiko braižymui


# susikuriu lentelę su sample ir total reads
lib_df <- data.frame(
  sample = colnames(count_matrix),
  total_reads = colSums(count_matrix)
)

lib_df <- merge(lib_df, meta, by = "sample")

# braižau grafiką
lib_size_plot <- ggplot(lib_df, aes(x = sample, y = total_reads / 1e6, fill = condition)) +
  geom_col(color = "black") +
  labs(
    title = "Library Sizes per Sample",
    subtitle = "Total mapped reads (raw counts)",
    x = "Sample",
    y = "Total counts (millions)",
    fill = "Treatment"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("normal" = "skyblue",
                               "tumor" = "tomato")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("library_sizes_plot.png", plot = lib_size_plot, width = 8, height = 6, dpi = 300)

# 2222222

# PEARSON CORRELATION

# log transformacija
log_counts <- log2(count_matrix + 1)

# skaičiuojam Pearson koreliaciją
cor_matrix <- cor(log_counts, method = "pearson")

# braižom heatmap
png("correlation/correlation_heatmap.png", width = 900, height = 800)
pheatmap(cor_matrix)
dev.off()

# 3333333

# PCA

# PCA skaičiavimas

log_counts_2 <- log_counts[apply(log_counts, 1, var) != 0, ]

pca <- prcomp(t(log_counts_2), scale. = TRUE)

# grafiko braižymas ir išsaugojimas

pca_df <- data.frame(
  sample = colnames(count_matrix),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2]
)

pca_df <- merge(pca_df, meta, by = "sample")


PCA_plot_R2 <- ggplot(pca_df, aes(PC1, PC2, color = condition, label = sample)) +
  geom_point(size = 3) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  geom_text_repel(show.legend = FALSE) +
  scale_color_manual(values = c("normal" = "skyblue",
                                "tumor" = "tomato")) +
  labs(title = "PCA plot")

ggsave("PCA_plot_R.png", plot = PCA_plot_R2, width = 8, height = 6, dpi = 300)

# 444444444

# Raw and normalized counts visually
# count_matrix → raw counts, norm_counts → normalized counts
# Normalization koreguoja library size / sequencing depth bias, kad mėginiai su daugiau reads neatrodytų dirbtinai labiau ekspresuoti.


dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = meta,
  design = ~ condition
)

dds <- estimateSizeFactors(dds)

norm_counts <- counts(dds, normalized = TRUE)

# RAW COUNTS
raw_df <- as.data.frame(log2(count_matrix + 1))
raw_df$gene <- rownames(raw_df)

raw_long <- melt(
  raw_df,
  id.vars = "gene",
  variable.name = "sample",
  value.name = "value"
)
raw_long$type <- "Raw counts"

# NORMALIZED COUNTS
norm_df <- as.data.frame(log2(norm_counts + 1))
norm_df$gene <- rownames(norm_df)

norm_long <- melt(
  norm_df,
  id.vars = "gene",
  variable.name = "sample",
  value.name = "value"
)
norm_long$type <- "Normalized counts"

# Sujungiam
all_data <- rbind(raw_long, norm_long)

# Prijungiam condition
all_data <- merge(all_data, meta, by = "sample")

all_data$type <- factor(
  all_data$type,
  levels = c("Raw counts", "Normalized counts")
)

# GRAFIKAS
raw_norm_plot <- ggplot(all_data, aes(x = sample, y = value, fill = condition)) +
  geom_boxplot(outlier.size = 0.3) +
  stat_summary(
    fun = median,
    geom = "point",
    shape = 95,
    size = 6,
    color = "black"
  ) +
  facet_wrap(~type) +
  theme_minimal() +
  scale_fill_manual(values = c(
    "normal" = "skyblue",
    "tumor" = "tomato"
  )) +
  labs(
    title = "Raw vs normalized counts",
    x = "Sample",
    y = "log2(counts + 1)",
    fill = "Condition"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

ggsave(
  "raw_vs_normalized_counts_colored.png",
  plot = raw_norm_plot,
  width = 10,
  height = 6,
  dpi = 300
)


# 5555555555

# Volcano plot


# Differential expression
dds <- DESeq(dds)
res <- results(dds)

res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

# pašalinam NA
res_df <- res_df[!is.na(res_df$padj) & !is.na(res_df$log2FoldChange), ]

# kategorijos
res_df$Direction <- "NS"
res_df$Direction[res_df$padj < 0.05 & res_df$log2FoldChange > 1] <- "Up"
res_df$Direction[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Down"

# parenku reikšmingiausius genus žymėjimui grafike
label_df <- subset(res_df, padj < 0.05 & abs(log2FoldChange) > 1)
label_df <- label_df[order(label_df$padj), ]
label_df <- head(label_df, 15)

# grafikas
volcano_plot2 <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Direction)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_text_repel(
    data = label_df,
    aes(label = gene),
    size = 3,
    show.legend = FALSE
  ) +
  scale_color_manual(values = c(
    "Up" = "#F8766D",
    "Down" = "#619CFF",
    "NS" = "grey70"
  )) +
  theme_minimal() +
  labs(
    title = "Volcano Plot — Tumor vs. Normal",
    subtitle = "padj < 0.05 and |log2FC| > 1 threshold shown as dashed lines",
    x = "log2 Fold Change (Tumor / Normal)",
    y = "-log10(adjusted p-value)"
  ) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

ggsave("volcano_plot.png", plot = volcano_plot2, width = 8, height = 6, dpi = 300)


# 66666666666

# MA plot

# shrink log2FC (stabilesni rezultatai)
resLFC <- lfcShrink(dds, coef = 2, type = "normal")

# į dataframe
ma_df <- as.data.frame(resLFC)
ma_df$gene <- rownames(ma_df)

# pašalinam NA
ma_df <- ma_df[!is.na(ma_df$padj), ]

# reikšmingumo kategorija
ma_df$significant <- ifelse(ma_df$padj < 0.05, "DE", "NS")

# MA plot
MA_plot <- ggplot(ma_df, aes(x = baseMean, y = log2FoldChange, color = significant)) +
  geom_point(size = 1.2, alpha = 0.7) +
  scale_x_log10() +
  scale_color_manual(values = c("NS" = "grey70", "DE" = "blue")) +
  geom_hline(yintercept = 0, color = "black") +
  theme_minimal() +
  labs(
    title = "MA Plot",
    x = "Mean of normalized counts",
    y = "log2 fold change",
    color = "Significance"
  ) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

# išsaugoti
ggsave("MA_plot.png", plot = MA_plot, width = 8, height = 6, dpi = 300)


# 77777777

# DEG heatmap


# VST transformacija – stabilizuoja variaciją tarp genų su skirtingu ekspresijos lygiu
vsd <- vst(dds, blind = FALSE)
vst_mat <- assay(vsd) 

# į rezultatų lentelę įdedame genų pavadinimus
res_df$gene <- rownames(res_df)

# pašaliname genus su NA padj ir surikiuojame pagal mažiausią padj
top_genes_df <- res_df[!is.na(res_df$padj), ]
top_genes_df <- top_genes_df[order(top_genes_df$padj), ]

# paimame top 50 diferencijuotai ekspresuotų genų
top_genes <- head(top_genes_df$gene, 50)

# iš VST matricos ištraukiame tik šiuos genus
heat_mat <- vst_mat[top_genes, ]

# Z-score normalizacija per genus (eilutes), kad matytume santykinius ekspresijos pokyčius
heat_mat_z <- t(scale(t(heat_mat)))

# paruošiame stulpelių anotaciją (sample + condition)
anno_col <- meta[, c("sample", "condition")]
rownames(anno_col) <- anno_col$sample
anno_col <- anno_col[colnames(heat_mat_z), , drop = FALSE]
anno_col$sample <- NULL

# nustatome spalvas sąlygoms
anno_colors <- list(
  condition = c(
    "normal" = "skyblue",
    "tumor" = "tomato"
  )
)

# išsaugome heatmap kaip paveikslą
png("DEG_heatmap_top50.png", width = 1000, height = 1200, res = 150)

# braižome heatmap
pheatmap(
  heat_mat_z,
  annotation_col = anno_col,
  annotation_colors = anno_colors,
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(50),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 7,
  border_color = NA,
  main = "Top 50 DE genes"
)

dev.off()


# 8888888888

# raw p-value versus adjusted p-values

head(res_df[, c("log2FoldChange", "pvalue", "padj")])

sum(res_df$pvalue < 0.05, na.rm = TRUE)
sum(res_df$padj < 0.05, na.rm = TRUE)


# ((((((((((((((((((()))))))))))))))))))
# {{{{{{{{{{{{{{{{{{{}}}}}}}}}}}}}}}}}}}
# DATA BIOLOGICAL ANALYSIS
# ((((((((((((((((((()))))))))))))))))))
# {{{{{{{{{{{{{{{{{{{}}}}}}}}}}}}}}}}}}}

# 111111111111

# TOP 3 upregulated genes and TOP  3 downregulated genes

res_df1 <- as.data.frame(res)
res_df1 <- na.omit(res_df1)

# significant genes
sig <- res_df1 |>
  filter(padj < 0.05)

# TOP 3 upregulated
top_up <- sig |>
  arrange(desc(log2FoldChange)) |>
  head(3)

# TOP 3 downregulated
top_down <- sig |>
  arrange(log2FoldChange) |>
  head(3)


top_up$symbol <- mapIds(org.Hs.eg.db,
                        keys = rownames(top_up),
                        column = "SYMBOL",
                        keytype = "ENSEMBL",
                        multiVals = "first")

top_up

top_down$symbol <- mapIds(org.Hs.eg.db,
                          keys = rownames(top_down),
                          column = "SYMBOL",
                          keytype = "ENSEMBL",
                          multiVals = "first")

top_down

# 22222222222

# GO over-representation analysis

# Reikia vector su reikšmingų genų ID
sig_genes <- rownames(sig)
head(sig_genes)

# Kadangi GO analizė dažniausiai naudoja ENTREZ ID, konvertuojame
gene_entrez <- mapIds(org.Hs.eg.db,
                      keys = sig_genes,
                      column = "ENTREZID",
                      keytype = "ENSEMBL",
                      multiVals = "first")

gene_entrez <- na.omit(gene_entrez)

# Padarome GO enrichment

ego <- enrichGO(
  gene = gene_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "ALL",          # BP + MF + CC
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

png("GO_dotplot.png", width = 1200, height = 900)

dotplot(ego, showCategory = 15)

dev.off()

# 33333333333333

# GSEA analysis over GO and MSIGdb

# Paruošiam visų genų reitinguotą sąrašą pagal log2FoldChange
gsea_df <- as.data.frame(res)
gsea_df <- na.omit(gsea_df)

gsea_df$ENSEMBL <- rownames(gsea_df)

gsea_df$ENTREZID <- mapIds(org.Hs.eg.db,
                           keys = gsea_df$ENSEMBL,
                           column = "ENTREZID",
                           keytype = "ENSEMBL",
                           multiVals = "first")

gsea_df <- gsea_df[!is.na(gsea_df$ENTREZID), ]

# jei tas pats ENTREZ pasikartoja, paliekam eilutę su didžiausiu |log2FC|
gsea_df <- gsea_df |>
  group_by(ENTREZID) |>
  slice_max(order_by = abs(log2FoldChange), n = 1, with_ties = FALSE) |>
  ungroup()

gene_list <- gsea_df$log2FoldChange
names(gene_list) <- gsea_df$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)

# GSEA su GO
gsea_go <- gseGO(
  geneList = gene_list,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "ALL",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  verbose = FALSE
)

png("GSEA_GO_dotplot.png", width = 1200, height = 900)
dotplot(gsea_go, showCategory = 15) +
  ggtitle("GSEA — GO (Tumor vs. Normal)")
dev.off()

# GSEA su MSigDB Hallmark
msig_h <- msigdbr(species = "Homo sapiens", category = "H")
term2gene_h <- msig_h[, c("gs_name", "entrez_gene")]

gsea_msig <- GSEA(
  geneList = gene_list,
  TERM2GENE = term2gene_h,
  pvalueCutoff = 0.05,
  verbose = FALSE
)

png("GSEA_MSigDB_dotplot.png", width = 1200, height = 900)
dotplot(gsea_msig, showCategory = 15) +
  ggtitle("GSEA — MSigDB Hallmark (Tumor vs. Normal)")
dev.off()

# Peržiūra konsolėj
head(as.data.frame(gsea_go)[, c("ID", "Description", "NES", "p.adjust")])
head(as.data.frame(gsea_msig)[, c("ID", "Description", "NES", "p.adjust")])

# 4444444444

# map visualisation of GSEA results 

# aPEAR įsidiegimo kelias
# install.packages(
# "https://cran.r-project.org/src/contrib/Archive/arules/arules_1.7-9.tar.gz",
# repos = NULL,
# type = "source"

# install.packages("ggrepel")
# remotes::install_github("kerseviciute/aPEAR", dependencies = FALSE)

# konvertuoju GSEA rezultatą į dataframe
gsea_df <- as.data.frame(gsea_go)

# sukuriu pathway tinklą
pear <- enrichmentNetwork(gsea_df)

# išsaugau grafiką
png("GSEA_aPEAR_map.png", width = 1200, height = 900)

plot(pear)

dev.off()



# kad paleistume visą R kodo failą per remote - source("HW1_R_dalis.R")
