library(cowplot)
library(DESeq2)
library(ggdendro)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(RColorBrewer)

source("scripts/functions/dge_gsea.R")

sample_info <- readRDS("data/coldata.rds")
counts <- readRDS("data/counts_matrix.rds")
dds <- readRDS("data/dds.rds")
vst_dds <- readRDS("data/vst_dds.rds")
geneset_ls <- readRDS("data/genesets.rds")
dge_result <- readRDS("data/dge_result.rds")

#### General description ####
# N patients in each group?
pts_included <- substr(colnames(counts), 1, 3)
table(sample_info$temperature[sample_info$patient %in% pts_included])
# 37 41
# 20 17

#### DGE ####
res_proc <- as.data.frame(dge_result$temperature$`41°C vs 37°C`)
# res_proc<- as.data.frame(dge_result$dose$`120mg vs 75mg`)
res_proc$HGNC <- rownames(res_proc)
genes_DE <- get_DE(
    logFC = res_proc$log2FoldChange,
    P = res_proc$padj,
    genes = res_proc$HGNC,
    th_logFC = 1,
    th_logP = 2,
    curve = 0
)
genes_to_label <- intersect(genes_DE$DE, geneset_ls$TFT_ls$HSF2_TARGET_GENES)
p_volc <- plot_volcano(
    DE_results = res_proc,
    gene = genes_to_label,
    labelAll = TRUE,
    labelCol = "black",
    useHyperbolicTH = FALSE,
    logFC_cu = 1,
    plotTH = TRUE
)
p_volc <- p_volc +
    theme(
        plot.title = element_text(hjust = 0.5, size = 8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", size = 0.2),
        axis.ticks = element_line(colour = "black", size = 0.2),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 7)
    ) +
    scale_x_continuous(name = "log2(Fold Change)", limits = c(-4, 4)) +
    scale_y_continuous(name = "-log10(Padj)")

# n DE?
length(genes_DE$up) # 44
length(genes_DE$down) # 5

sort(genes_DE$up)
# [1] "ACTRT3"   "AHSA1"    "APOA1"    "BAG3"     "C11orf84" "C17orf67" "C19orf84" "CACYBP"
# [9] "CKS2"     "CLIC2"    "DEDD2"    "DLX2"     "DNAJA1"   "DNAJA4"   "DNAJB1"   "DNAJB4"
# [17] "FBXL14"   "FKBP4"    "GMPR"     "HSPA4L"   "HSPA6"    "HSPB1"    "HSPH1"    "IER5"
# [25] "JMJD6"    "LGR5"     "LMAN2L"   "LRIF1"    "MB21D1"   "MRPS6"    "NKRF"     "OSBP2"
# [33] "PNMAL1"   "PODXL2"   "RND1"     "SLC5A3"   "SPR"      "STIP1"    "TRAPPC5"  "TRIM16"
# [41] "TSHZ1"    "TSPYL4"   "USPL1"    "ZFAND2A"

sort(genes_DE$down)
# [1] "AIRE"     "ALDH1A3"  "C17orf58" "KLHL29"   "MAGEA9B"

#### GSEA ####
stat <- get_GSEA_stat(res_proc)
genesets <- c("Ha_ls", "TFT_ls", "GO_MF_ls", "Rea_ls")
pws <- c("HALLMARK_E2F_TARGETS", "HSF2_TARGET_GENES", "GOMF_HEAT_SHOCK_PROTEIN_BINDING", "REACTOME_HSF1_ACTIVATION")

# GSEA
GSEA_res <- sapply(genesets, \(gs) do_fGSEA(geneset_ls[[gs]], stat), simplify = FALSE)

# RS
p_RS_ls <- mapply(\(pw, gs) {
    p_pw <- format(signif(GSEA_res[[gs]][pathway == pw, "padj"], 3), scientific = TRUE)
    plotEnrichment(geneset_ls[[gs]][[pw]], stat, ticksSize = .1, ticksLength = 0.2) +
        # ggtitle(paste0(pw)) +
        ggtitle(paste0(pw, "\n(Padj=", p_pw, ")")) +
        geom_line(size = .5, col = "green") +
        theme(
            plot.title = element_text(hjust = 0.5, size = 7, face = "italic"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black", size = 0.2),
            axis.ticks = element_line(colour = "black", size = 0.2),
            axis.text = element_text(size = 6),
            axis.title = element_text(size = 7),
        ) +
        scale_x_continuous(name = "Rank") +
        scale_y_continuous(name = "Enrichment Score")
}, pws, genesets, SIMPLIFY = FALSE)

# Save GSEA + DGE
WriteXLS::WriteXLS(
    list(res_proc, GSEA_res$TFT_ls, GSEA_res$GO_MF_ls, GSEA_res$Rea_ls, GSEA_res$Ha_ls),
    "results/tables/manuscript_GSEA_table.xlsx",
    row.names = TRUE,
    SheetNames = c(
        "DGE analysis",
        "GSEA - Transcr. factor targets",
        "GSEA - GO Molecular Function",
        "GSEA - Reactome",
        "GSEA - Hallmark"
    )
)

#### Heatmap ####
# obtain relevant geneset genes
hsf2_target_genes <- geneset_ls$TFT_ls$HSF2_TARGET_GENES
reactome_hsf1 <- geneset_ls$Rea_ls$REACTOME_HSF1_ACTIVATION
gomf_hsp_binding <- geneset_ls$GO_MF_ls$GOMF_HEAT_SHOCK_PROTEIN_BINDING

# obtain and filter relevant temperature result
result <- dge_result[["temperature"]][["41°C vs 37°C"]]
result <- result[result$padj < 0.01 & abs(result$log2FoldChange) > 1, ]
result <- result[order(result$padj), ]

top_count <- min(nrow(result), 50)
top_genes <- result[1:top_count, "symbol"]

rownames(vst_dds) <- convert_ens(rownames(vst_dds))

matrix <- assay(vst_dds)[rownames(assay(vst_dds)) %in% top_genes, ]
# mean-center values
matrix <- matrix - rowMeans(matrix)

anno_data <- as.data.frame(colData(dds)[order(dds$temperature, dds$dose), c("dose", "temperature")])
samples_order <- rownames(anno_data)

# helper function
prepare_data <- function(data) {
    data <- as.data.frame(data)
    data$gene <- rownames(data)
    reshape2::melt(data, id.vars = "gene")
}

grid_color <- "grey40"

# cluster genes
clust_genes <- hclust(dist(matrix, method = "euclidean"), method = "complete")
genes_order <- clust_genes$labels[clust_genes$order]

# create gene dendrogram
hcdata <- dendro_data(clust_genes)
p_dendro_genes <- ggplot() +
    geom_segment(
        data = segment(hcdata),
        aes(x = x, y = y, xend = xend, yend = yend)
    ) +
    theme_void() +
    scale_y_reverse(expand = c(0, 0)) +
    scale_x_discrete(expand = c(0.0125, 0)) +
    coord_flip()

# create main heatmap plot
pd_vsd <- prepare_data(matrix)
pd_vsd$variable <- factor(pd_vsd$variable, levels = rownames(anno_data))
pd_vsd$gene <- factor(pd_vsd$gene, levels = genes_order)
p_heatmap <- ggplot(
    pd_vsd,
    aes(
        y = gene,
        x = variable
    )
) +
    geom_tile(
        color = NA,
        aes(fill = value)
    ) +
    scale_fill_gradient2(
        "VST counts",
        low = "gold1", mid = "white", high = "blueviolet",
        midpoint = 0,
    ) +
    theme_minimal() +
    theme(
        axis.text = element_text(size = 12, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6)
    ) +
    guides(fill = guide_colourbar(
        barheight = unit(5, "cm"),
        ticks.colour = "black",
        ticks.linewidth = 1
    ))

# create annotation plot
anno_data$id <- rownames(anno_data)
pd_anno <- reshape2::melt(anno_data, id.vars = "id")
pd_anno$id <- factor(pd_anno$id, levels = samples_order)

p_anno <- ggplot(pd_anno, aes(x = id, y = variable)) +
    geom_tile(
        data = pd_anno[pd_anno$variable == "dose", ],
        color = NA,
        aes(
            fill = factor(value, levels = c("75", "100", "120"))
        )
    ) +
    scale_fill_manual(
        name = "dose",
        values = c(brewer.pal(4, "Greens")[2:4]),
        labels = c("75 mg/m2", "100 mg/m2", "120 mg/m2")
    ) +
    ggnewscale::new_scale_fill() +
    geom_tile(
        data = pd_anno[pd_anno$variable == "temperature", ],
        color = NA,
        aes(
            fill = factor(value, levels = c("37", "41"))
        )
    ) +
    scale_fill_manual(
        name = "temperature",
        values = c(c("#377EB8", "#E41A1C")),
        labels = c("37 °C", "41 °C")
    ) +
    scale_y_discrete(position = "right") +
    theme_minimal() +
    theme(
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6)
    )

# create gene membership plot
gm_data <- data.frame(
    "gene" = top_genes,
    "HSF2 target genes" = top_genes %in% hsf2_target_genes,
    "HSF1 activation" = top_genes %in% reactome_hsf1,
    "HSP binding" = top_genes %in% gomf_hsp_binding
)
gm_data <- gm_data[match(genes_order, gm_data$gene), ]
gm_data <- reshape2::melt(gm_data, id.vars = "gene")

p_gm <- ggplot(
    gm_data,
    aes(
        y = factor(gene, levels = genes_order),
        x = variable
    )
) +
    geom_tile(
        color = grid_color,
        aes(fill = value)
    ) +
    scale_fill_manual(
        "geneset\nmembership",
        values = c("FALSE" = "white", "TRUE" = "grey40"),
        labels = c("FALSE" = "No", "TRUE" = "Yes")
    ) +
    scale_y_discrete(position = "right") +
    scale_x_discrete(
        labels = \(l) gsub("\\.", " ", l)
    ) +
    theme_minimal() +
    theme(
        legend.position = "none",
        axis.text.y = element_text(size = 6, color = "black", face = "italic"),
        axis.text.x = element_text(color = "black", angle = 45, hjust = 1, size = 6),
        axis.title = element_blank()
    )

layout <- "
    #AAAAAAAA##
    CBBBBBBBBD#
    CBBBBBBBBD#
    CBBBBBBBBD#
    CBBBBBBBBD#
    CBBBBBBBBD#
    CBBBBBBBBD#
    CBBBBBBBBD#
    CBBBBBBBBD#
    CBBBBBBBBD#
    CBBBBBBBBD#
    CBBBBBBBBD#
    CBBBBBBBBD#
    CBBBBBBBBD#
    CBBBBBBBBD#
    CBBBBBBBBD#
    CBBBBBBBBD#"

p_hm <- wrap_plots(
    A = p_anno,
    B = p_heatmap,
    C = p_dendro_genes,
    D = p_gm,
    design = layout
) +
    plot_layout(guides = "collect") &
    theme(
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()
    )

#### Create figure ####
# 3 RS plots
# Shorter titles
p_RS_ls$REACTOME_HSF1_ACTIVATION$labels$title <- gsub(
    "ACTIVATION", "ACT", p_RS_ls$REACTOME_HSF1_ACTIVATION$labels$title
)
p_RS_ls$GOMF_HEAT_SHOCK_PROTEIN_BINDING$labels$title <- gsub(
    "HEAT_SHOCK", "HS", p_RS_ls$GOMF_HEAT_SHOCK_PROTEIN_BINDING$labels$title
)

p_rs <- plot_grid(
    p_RS_ls$HSF2_TARGET_GENES,
    p_RS_ls$REACTOME_HSF1_ACTIVATION,
    p_RS_ls$GOMF_HEAT_SHOCK_PROTEIN_BINDING,
    ncol = 3
)

# Merge with volcano
p_volcrs <- plot_grid(
    p_volc, p_rs,
    ncol = 1,
    rel_heights = c(2, 1),
    scale = .95,
    labels = "AUTO"
)

# Merge with cibersort results
load("data/boxplot_cibersort_x_lm22_relative.RData")
p_cs <- plot_grid(icd_boxplot, scale = .95, labels = "D")

p <- plot_grid(
    p_volcrs, p_cs,
    ncol = 2,
    rel_widths = c(2, 1)
)

# Merge with heatmap
p_hm2 <- plot_grid(p_hm, scale = .95, labels = "C")
p <- plot_grid(
    p, p_hm2,
    ncol = 1,
    rel_heights = c(3, 4),
    labels = c(NA, "C")
)

ggsave("results/figs/fig_11.pdf", p, width = 178, height = 265, units = "mm")
ggsave("results/figs/fig_11.png", p, width = 178, height = 265, units = "mm", bg = "#ffffff")
