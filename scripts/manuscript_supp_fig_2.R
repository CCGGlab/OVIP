library(cowplot)
library(DESeq2)
library(ggplot2)

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
table(sample_info$dose[sample_info$patient %in% pts_included])
# 100 120  75
# 6  15  16

#### DGE ####
res_proc <- as.data.frame(dge_result$dose$`120mg vs 75mg`)
res_proc$HGNC <- rownames(res_proc)
genes_DE <- get_DE(
    logFC = res_proc$log2FoldChange,
    P = res_proc$padj,
    genes = res_proc$HGNC,
    th_logFC = 1,
    th_logP = 2,
    curve = 0
)
p_volc <- plot_volcano(
    DE_results = res_proc,
    gene = NA,
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

# Add outliers?
res_proc$isOut <- abs(res_proc$log2FoldChange) > 4 & res_proc$padj < 0.01
res_proc2 <- res_proc[res_proc$isOut, ]
res_proc2$log2FoldChange <- sign(res_proc2$log2FoldChange) * 4
p_volc <- p_volc + geom_point(aes(x = res_proc2$log2FoldChange[1], y = -log10(res_proc2$padj)[1]), shape = 2)
p_volc <- p_volc + geom_point(aes(x = res_proc2$log2FoldChange[2], y = -log10(res_proc2$padj)[2]), shape = 2)

# n DE?
length(genes_DE$up) # 5
length(genes_DE$down) # 2

sort(genes_DE$up)
# "CCL22"  "CD52"   "DUSP10" "LRRC15" "RNF5"

sort(genes_DE$down)
# "DLX2"    "MAGEA9B"

#### GSEA ####
stat <- get_GSEA_stat(res_proc)
pws <- c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_E2F_TARGETS", "HALLMARK_DNA_REPAIR")
# GSEA
GSEA_Ha <- do_fGSEA(geneset_ls$Ha_ls, stat)
# running score plots
p_RS_ls <- sapply(pws, \(pw) {
    p_pw <- format(signif(GSEA_Ha[GSEA_Ha$pathway == pw, "padj"], 3), scientific = TRUE)
    plotEnrichment(geneset_ls$Ha_ls[[pw]], stat, ticksSize = .1) +
        # ggtitle(paste0(pw)) +
        ggtitle(paste0(pw, "\n(Padj=", p_pw, ")")) +
        geom_line(linewidth = 0.5, col = "green") +
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
}, simplify = FALSE)

#### Save DGE + GSEA ####
WriteXLS::WriteXLS(
    c("res_proc", "GSEA_Ha"),
    "results/tables/manuscript_GSEA_table_dose.xlsx",
    row.names = TRUE,
    SheetNames = c("DGE analysis", "GSEA - Hallmark")
)

#### Create figure ####
p_RS_ls$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION$labels$title <-
    gsub(
        "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
        "HALLMARK_EPITH_MES_TRANSITION",
        p_RS_ls$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION$labels$title
    )
p_rs <- plot_grid(
    p_RS_ls$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION,
    p_RS_ls$HALLMARK_E2F_TARGETS,
    p_RS_ls$HALLMARK_DNA_REPAIR,
    ncol = 1
)

p <- plot_grid(
    p_volc, p_rs,
    rel_widths = c(2, 1),
    labels = "AUTO",
    scale = .9
)

p <- plot_grid(
    p, NA,
    ncol = 1,
    rel_heights = c(1, 1)
)

ggsave("results/figs/supp_fig_2.pdf", p, width = 178, height = 265, units = "mm")
