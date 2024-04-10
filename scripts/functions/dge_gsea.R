library(AnnotationDbi)
library(EnsDb.Hsapiens.v86)
library(fgsea)
library(ggplot2)
library(ggrepel)
library(grid)
library(gridExtra)

convert_ens <- function(ens, target = "symbol") {
    target <- tolower(target)
    stopifnot(
        "Provided vector contains invalid ensembl identifiers."
            = any(stringr::str_starts(ens, "^ENS.*?G")),
        "Target must be one of symbol, entrez" = target %in% c("symbol", "entrez")
    )
    if (target == "symbol") {
        mapIds(EnsDb.Hsapiens.v86, keys = ens, column = "SYMBOL", keytype = "GENEID")
    } else if (target == "entrez") {
        mapIds(EnsDb.Hsapiens.v86, keys = ens, column = "ENTREZID", keytype = "GENEID")
    }
}

get_DE <- function(logFC, P, genes, th_logFC = 1, th_logP = -log10(0.05), curve = 1) {
    idx_DE <- which(
        -log10(P) > sapply(logFC, \(x) get_hyperbolic_th_P(x, th_logFC = th_logFC, th_logP = th_logP, curve = curve))
    )
    idx_up <- intersect(idx_DE, which(logFC > 0))
    idx_down <- intersect(idx_DE, which(logFC < 0))
    genes_up <- genes[idx_up]
    genes_down <- genes[idx_down]
    genes_DE <- genes[idx_DE]
    res <- list(up = genes_up, down = genes_down, DE = genes_DE)
    return(res)
}

get_hyperbolic_th_P <- function(logFC, th_logFC = 1, th_logP = -log10(0.05), curve = 1) {
    if (abs(logFC) < th_logFC) {
        res <- Inf
    } else {
        res <- th_logP + (curve / sqrt(logFC^2 - th_logFC^2))
    }
    return(res)
}

get_GSEA_stat <- function(DE_results, isMice = FALSE, MGI_to_HGNC = NULL, isProt = FALSE) {
    res_proc <- DE_results
    if (isProt) {
        colnames(res_proc) <- c("log2FoldChange", "pvalue", "padj", "stat")
        res_proc$stat <- -1 * res_proc$stat
    }
    stat <- res_proc$stat + 10e-10 * res_proc$log2FoldChange # Addition ro prevent ties
    names(stat) <- rownames(res_proc)

    if (isMice) {
        MGI_to_HGNC <- MGI_to_HGNC[!duplicated(MGI_to_HGNC$MGI.symbol), ]
        rownames(MGI_to_HGNC) <- MGI_to_HGNC$MGI.symbol
        names(stat) <- MGI_to_HGNC[names(stat), "HGNC.symbol"]
    }

    stat <- stat[!duplicated(names(stat))]
    stat <- sort(stat, decreasing = TRUE)
    stat
}

do_fGSEA <- function(db, stat, minSize = 15, maxSize = 500) {
    Res <- fgsea::fgseaMultilevel(
        pathways = db,
        stats = stat,
        minSize = minSize,
        maxSize = maxSize,
        eps = 0
    )
    Res <- Res[order(Res$pval), ]
    Res
}

# Plot volcano
plot_volcano <- function(
    DE_results,
    gene,
    isProt = FALSE,
    useHyperbolicTH = FALSE,
    labelAll = FALSE,
    labelCol = "red",
    p_cu = 0.01,
    logFC_cu = 2,
    curve = 0.5,
    plotTH = TRUE,
    plot_nominal_p = FALSE,
    th_nominal_p = FALSE
) {
    # Get selected data
    res_proc <- DE_results
    if (isProt) colnames(res_proc) <- c("log2FoldChange", "pvalue", "padj", "stat")
    res_proc <- res_proc[!is.na(res_proc[, "padj"]), ] # remove NA (i.e. no expression before/after treatment)
    # if q=0, put at lowest possible value
    res_proc[res_proc[, "padj"] == 0, "padj"] <- min(res_proc[res_proc[, "padj"] != 0, "padj"]) / 100
    # if q=0, put at lowest possible value
    res_proc[res_proc[, "pvalue"] == 0, "pvalue"] <- min(res_proc[res_proc[, "pvalue"] != 0, "pvalue"]) / 100

    # DE thresholds
    if (th_nominal_p) {
        if (useHyperbolicTH) {
            res_proc$isDE <- rownames(res_proc) %in% get_DE(
                res_proc$log2FoldChange, res_proc$pvalue, rownames(res_proc), logFC_cu, -log10(p_cu), curve
            )$DE
        } else {
            res_proc$isDE <- abs(res_proc$log2FoldChange) >= logFC_cu & -log10(res_proc$pvalue) >= -log10(p_cu)
        }
    } else {
        if (useHyperbolicTH) {
            res_proc$isDE <- rownames(res_proc) %in% get_DE(
                res_proc$log2FoldChange, res_proc$padj, rownames(res_proc), logFC_cu, -log10(p_cu), curve
            )$DE
        } else {
            res_proc$isDE <- abs(res_proc$log2FoldChange) >= logFC_cu & -log10(res_proc$padj) >= -log10(p_cu)
        }
    }

    # Gene ids
    res_proc$gene_id <- rownames(res_proc)
    res_proc$isLabel <- res_proc$gene_id %in% gene

    # Plot
    res_proc <- res_proc[order(res_proc$isDE), ] # Blue before grey
    if (plot_nominal_p) {
        p <- ggplot(res_proc, aes(x = log2FoldChange, y = -log10(pvalue), color = isDE, key = gene_id))
    } else {
        p <- ggplot(res_proc, aes(x = log2FoldChange, y = -log10(padj), color = isDE, key = gene_id))
    }
    p <- p +
        geom_point() +
        scale_color_manual(values = c("#d0d3d4", "#3498db")) +
        theme(legend.position = "none")

    if (length(gene) >= 1) p <- p + geom_point(data = res_proc[gene, ], colour = labelCol) # this adds a red point

    if (labelAll == TRUE) {
        p <- p +
            geom_text_repel(
                data = subset(res_proc, res_proc$isLabel),
                aes(label = gene_id, fontface = 3),
                size = 2,
                colour = "black",
                box.padding = 0.5,
                point.padding = 0.5,
                max.overlaps = 25,
                segment.color = "grey50"
            )
        # this adds a label for the red point
        # p <- p + geom_text(data = res_proc[gene, ], label = gene, vjust = 0, hjust = 0, colour = labelCol)
    }

    # Hyperbolic treshold?
    if (plotTH == TRUE) {
        if (useHyperbolicTH) {
            x_lim <- c(round(min(res_proc$log2FoldChange) - 1), round(max(res_proc$log2FoldChange) + 1))
            if (plot_nominal_p) {
                y_lim <- c(0, round(max(-log10(res_proc$pvalue)) + 1))
            } else {
                y_lim <- c(0, round(max(-log10(res_proc$padj)) + 1))
            }
            th_data <- data.frame(
                x = seq(x_lim[1], x_lim[2], 0.01),
                y = get_hyperbolic_th_P(
                    seq(x_lim[1], x_lim[2], 0.01),
                    th_logFC = logFC_cu, th_logP = -log10(p_cu), curve = curve
                ),
                gene_id = "x",
                isDE = FALSE
            )
            th_data <- th_data[th_data$y < y_lim[2], ]
            p <- p +
                geom_line(data = th_data[th_data$x < 0, ], aes(x = x, y = y), linetype = "dashed", colour = "grey20") +
                geom_line(data = th_data[th_data$x > 0, ], aes(x = x, y = y), linetype = "dashed", colour = "grey20")
        } else {
            p <- p +
                geom_hline(yintercept = c(-log10(p_cu)), linetype = "dashed", colour = "grey20") +
                geom_vline(xintercept = c(-logFC_cu, logFC_cu), linetype = "dashed", colour = "grey20")
        }
    }

    # return
    return(p)
}

do_GSEA2 <- function(genes_retrieved, genes_all, GSEA_db, min_genes, isList = FALSE) {
    if (isList) {
        GSEA_names <- names(GSEA_db)
    } else {
        GSEA_names <- rownames(GSEA_db)
    }
    GSEA_table <- matrix(
        NA, length(GSEA_names), 9,
        dimnames = list(
            GSEA_names,
            c("n_genes_pw", "n_genes_pw_neg", "n_genes_pw_pos", "prop_neg", "prop_pos", "OR", "p", "q", "genes")
        )
    )
    for (i in 1:nrow(GSEA_table)) {
        # cat(i," ")
        if (isList) {
            genes_temp <- unlist(GSEA_db[[i]])
        } else {
            genes_temp <- unique(as.character(GSEA_db[i, -1]))
            genes_temp <- genes_temp[genes_temp != ""]
        }
        # genes_overlap<- genes_temp[genes_temp%in%genes_all]
        genes_overlap <- intersect(genes_temp, genes_all)
        GSEA_table[i, "n_genes_pw"] <- length(genes_overlap)
        if (length(genes_overlap) < min_genes) {
            # No use analysing
            next
        } else {
            isRetrieved <- genes_all %in% genes_retrieved
            inDB <- factor(genes_all %in% genes_overlap, levels = c(FALSE, TRUE))
            compare_pw_temp_t <- table(isRetrieved, inDB)
            GSEA_table[i, c("n_genes_pw_neg", "n_genes_pw_pos")] <- compare_pw_temp_t[, "TRUE"]
            GSEA_table[i, c("prop_neg", "prop_pos")] <- round(100 * prop.table(compare_pw_temp_t, 1)[, "TRUE"], 1)
            GSEA_table[i, "genes"] <- paste(genes_all[isRetrieved & inDB == TRUE], collapse = ",")
            GSEA_table[i, "OR"] <- fisher.test(compare_pw_temp_t, alternative = "greater")$estimate
            GSEA_table[i, "p"] <- fisher.test(compare_pw_temp_t, alternative = "greater")$p.value
        }
    }
    GSEA_table <- GSEA_table[!is.na(GSEA_table[, "p"]) & !duplicated(rownames(GSEA_table)), ]
    # GSEA_table[,"q"]<- qvalue(as.numeric(GSEA_table[,"p"]))$qvalues # To use Storey method (library qvalue)
    GSEA_table[, "q"] <- p.adjust(as.numeric(GSEA_table[, "p"]), "fdr")
    GSEA_table <- as.data.frame(GSEA_table[order(as.numeric(GSEA_table[, "p"])), ])
    return(GSEA_table)
}


# Rewrite fGSEA to make ticklengths adjustable
plotEnrichment <- function(
    pathway,
    stats,
    gseaParam = 1,
    ticksSize = 0.2,
    ticksLength = 0.1
) {
    rnk <- rank(-stats)
    ord <- order(rnk)

    statsAdj <- stats[ord]
    statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
    statsAdj <- statsAdj / max(abs(statsAdj))

    pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
    pathway <- sort(pathway)

    gseaRes <- calcGseaStat(statsAdj,
        selectedStats = pathway,
        returnAllExtremes = TRUE
    )

    bottoms <- gseaRes$bottoms
    tops <- gseaRes$tops

    n <- length(statsAdj)
    xs <- as.vector(rbind(pathway - 1, pathway))
    ys <- as.vector(rbind(bottoms, tops))
    toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))

    # Getting rid of NOTEs
    x <- y <- NULL
    g <- ggplot(toPlot, aes(x = x, y = y)) +
        geom_point(color = "green", size = 0.1) +
        geom_hline(yintercept = max(tops), colour = "red", linetype = "dashed") +
        geom_hline(yintercept = min(bottoms), colour = "red", linetype = "dashed") +
        geom_hline(yintercept = 0, colour = "black") +
        geom_line(color = "green") +
        theme_bw() +
        geom_segment(
            data = data.frame(x = pathway),
            mapping = aes(
                x = x, y = -ticksLength / 2,
                xend = x, yend = ticksLength / 2
            ),
            size = ticksSize
        ) +
        theme(
            panel.border = element_blank(),
            panel.grid.minor = element_blank()
        ) +
        labs(x = "rank", y = "enrichment score")
    g
}
