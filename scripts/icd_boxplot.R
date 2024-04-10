library(ggplot2)

source("scripts/functions/icd.R")

#### data preparation ####
# annotation data
coldata <- readRDS("data/coldata.rds")
# prepare and filter CIBERSORT data (CIBERSORTx, LM22, relative mode)
icd_data <- read.delim(
    "data/CIBERSORTx_Job4_Adjusted.txt",
    sep = "\t", header = TRUE
)
# define column subsets based on icd_data colnames
define_column_sets(icd_data, "x")
# simplify sample names (remove A-suffix)
icd_data[[sample_column]] <- gsub("A$", "", icd_data[[sample_column]])
icd_data[[sample_column]] <- as.factor(icd_data[[sample_column]])
# inner join annotation
icd_data <- base::merge(x = icd_data, y = coldata, by.x = sample_column, by.y = "patient")
# filter 
# a) only significant patients
# medium dose samples are kept for temperature tests and filtered out later for dose tests
icd_data <- icd_data[icd_data$P.value <= 0.05, ]

#### figure: grouped cell types ####
icd_data <- prepare_icd_data(icd_data, "x", melt = TRUE)
cell_groups <- list(
    "CD4 T cells" = c(
        cell_columns[grep("^T.cells.CD4", cell_columns)], 
        "T.cells.follicular.helper", "T.cells.regulatory..Tregs."
    ),
    "CD8 T cells" = c(
        cell_columns[grep("^T.cells.CD8", cell_columns)]
    ),
    "B cells" = c(
        cell_columns[grep("^B.cells", cell_columns)], 
        "Plasma.cells"
    ),
    "NK cells" = c(
        cell_columns[grep("^NK.cells", cell_columns)]
    ),
    "Monocytes/macrophages" = c(
        cell_columns[grep("^Macrophage|^Dendritic.cells", cell_columns)], 
        "Monocytes"
    ),
    "Granulocytes" = c(
        cell_columns[grep("^Mast.cells", cell_columns)], 
        "Neutrophils", "Eosinophils"
    )
)
# cell_columns[which(!cell_columns %in% cell_groups$value)]
# not classified: T.cells.gamma.delta

cg_m <- reshape2::melt(cell_groups)

icd_data$cell_group <- cg_m[match(icd_data$cell_type, cg_m$value, NA), "L1"]
icd_data$cell_group <- factor(icd_data$cell_group, levels = names(cell_groups))

icd_data <- na.omit(icd_data)

# check normality
setNames(lapply(unique(icd_data$cell_group), \(g) {
    shapiro.test(icd_data[icd_data$cell_group == g, "value"])
}), unique(icd_data$cell_group))
# --> neither of the cell groups has normally distributed values

group_sums <- aggregate(value ~ Mixture + cell_group, icd_data, sum)
icd_sum <- merge(unique(icd_data[, c("Mixture", "temperature")]), group_sums, by = "Mixture")
icd_sum <- icd_sum[order(icd_sum$Mixture), ]

# pairwise testing by group
# group values summed per patient
icd_pwt <- icd_wilcox_pairwise(
    icd_sum, 
    group_column = "cell_group", 
    condition = "temperature", 
    value_column = "value"
)
icd_pwt <- Reduce(rbind, icd_pwt)
# prepare for stat_pvalue_manual
colnames(icd_pwt)[which(colnames(icd_pwt) %in% c("37°C", "41°C"))] <- c("group1", "group2")
icd_pwt$group1 <- "37°C"
icd_pwt$group2 <- "41°C"

icd_boxplot <- ggplot(
    icd_sum,
    aes(
        x = temperature,
        y = value
    )
) +
    geom_violin(
        aes(fill = cell_group)
    ) +
    geom_dotplot(
        binaxis = "y",
        stackdir = "center",
        stackratio = 1.5,
        dotsize = 0.5
    ) +
    stat_summary(
        color = "red3",
        width = 0.33,
        size = 0.33,
        fun = "median",
        geom = "crossbar"
    ) +
    stat_pvalue_manual(
        icd_pwt,
        label = "p = {sprintf('%5.3f', P)}",
        remove.bracket = TRUE,
        y.position = 0.95,
        hjust = 1,
        # position = position_dodge(width = 0.5),
        size = 2
    ) +
    facet_wrap(
        . ~ factor(cell_group, levels = levels(icd_sum$cell_group)),
        nrow = 3,
        ncol = 2
    ) +
    labs(
        y = "score"
    ) +
    scale_fill_brewer(
        palette = "Set2"
    ) +
    scale_y_continuous(name = "Cibersort score", limits = c(0, 1)) +
    theme_minimal() +
    theme(
        axis.text = element_text(color = "black", size = 6),
        axis.title = element_text(color = "black", size = 7),
        strip.text = element_text(color = "black", size = 7),
        strip.clip = "off",
        legend.position = "none"
    ) +
    guides(
        fill = guide_legend(
            title = "cell group",
            direction = "horizontal",
            nrow = 1
        )
    )

save(icd_sum, icd_pwt, icd_boxplot, file = "data/boxplot_cibersort_x_lm22_relative.RData")
