prepare_icd_data <- function(data, mode, subset = "both", relevel.factors = TRUE, melt = FALSE) {
    stopifnot(
        "Argument subset must be one of cells, stats, both" =
            subset %in% c("cells", "stats", "both")
    )
    # resolve absolute score name ambiguity
    if ("Absolute.score" %in% colnames(data)) {
        colnames(data)[which(colnames(data) == "Absolute.score")] <- "Abs.score.(rel.mode)"
    }

    define_column_sets(data, mode)

    if (relevel.factors) {
        # relevel factors correctly
        data$dose <- factor(
            data$dose,
            levels = c("75", "100", "120"),
            labels = c("75mg", "100mg", "120mg"),
            ordered = TRUE
        )
        data$temperature <- factor(
            data$temperature,
            levels = c("37", "41"),
            labels = c("37°C", "41°C"),
            ordered = TRUE
        )
        data$treatment <- factor(
            data$treatment,
            levels = c("75_37", "75_41", "100_37", "100_41", "120_37", "120_41"),
            labels = c("75mg 37°C", "75mg 41°C", "100mg 37°C", "100mg 41°C", "120mg 37°C", "120mg 41°C"),
            ordered = TRUE
        )
    }

    if (subset == "cells") {
        data <- data[, c(id_columns, cell_columns)]
    } else if (subset == "stats") {
        idxs <- c(
            which(colnames(data) %in% id_columns),
            which(colnames(data) %in% stats_columns)
        )
        data <- data[, idxs]
    }

    if (melt) {
        data <- reshape2::melt(data, id.vars = id_columns)
        colnames(data)[which(colnames(data) == "variable")] <- "cell_type"
    }

    data
}

define_column_sets <- function(data, mode) {
    stopifnot("Mode must be one of x, legacy" = mode %in% c("x", "legacy"))
    stopifnot("data data appears to be molten" = !all(c("cell_type", "value") %in% colnames(data)))

    stats_columns <<- switch(mode,
        "x" = c(
            "P.value", "Correlation",
            "RMSE", "Absolute.score..sig.score."
        ),
        "legacy" = c(
            "Pearson.Correlation", "P.value",
            "Absolute.Score", "Absolute.score",
            "Abs.score.(rel.mode)", "RMSE"
        )
    )
    sample_column <<- switch(mode,
        "x" = "Mixture",
        "legacy" = "Input.Sample"
    )
    id_columns <<- c(
        sample_column, "dose", "temperature", "treatment"
    )
    cell_columns <<- base::setdiff(
        base::setdiff(colnames(data), id_columns),
        stats_columns
    )
}

icd_wilcox_pairwise <- function(data, group_column, condition_column, value_column = "value") {
    stopifnot(
        "Group column must be a column of data." = group_column %in% colnames(data),
        "Condition column must be a column of data." = condition_column %in% colnames(data),
        "Value column must be a column of data." = value_column %in% colnames(data)
    )

    # test fails if input is a tibble
    data <- as.data.frame(data)

    setNames(lapply(levels(data[[group_column]]), \(c) {
        groups <- data[data[[group_column]] == c, condition_column]
        # here we only test the group-wise comparisons of 'condition_column',
        # FDR correction should be done on all tests combined,
        # i.e. on the result of this function
        test <- pairwise.wilcox.test(
            data[data[[group_column]] == c, value_column],
            groups,
            p.adjust.method = "none"
        )
        # add median + IQR to table
        res <- reshape2::melt(test$p.value)
        colnames(res)[which(colnames(res) == "value")] <- "P"
        u_groups <- as.character(unique(groups))
        group_values <- sapply(u_groups, \(i) {
            data[data[[group_column]] == c & data[[condition_column]] == i, value_column]
        }, simplify = FALSE)
        group_quantiles <- sapply(u_groups, \(i) {
            quantile(group_values[[i]])
        }, simplify = FALSE)
        res <- cbind(
            res,
            setNames(lapply(u_groups, \(i) {
                paste0(
                    sprintf("%.3f", median(group_values[[i]])),
                    " (",
                    sprintf("%.3f", group_quantiles[[i]]["25%"]),
                    "-",
                    sprintf("%.3f", group_quantiles[[i]]["75%"]),
                    ")"
                )
            }), u_groups)
        )
        # add Ns to table
        freq <- table(groups)
        freq_names <- unlist(Map(paste, "n", u_groups, sep = "_"))
        res <- cbind(
            res,
            setNames(list(
                as.numeric(freq[as.character(res$Var1)]),
                as.numeric(freq[as.character(res$Var2)])
            ), freq_names)
        )
        res <- cbind(res, setNames(list(c), group_column))
        res[, c(group_column, u_groups, freq_names, "P")]
    }), levels(data[[group_column]]))
}
