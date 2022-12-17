#' Clustering by PCA
#'
#' Runs a set of features through PCA and saves results to reducedDims slot of input TapestriObject.
#'
#' @param TapestriExperiment TapestriExperiment object
#' @param alt.exp Chr string indicating altExp to use, NULL uses top-level experiment. Default "alleleFrequency".
#' @param assay Chr string indicating assay to use. NULL (default) selects first listed assay.
#' @param sd.min.threshold Numeric, minimum threshold for allelefreq.sd. Increase to run PCA on fewer, increasingly variable dimensions. Set to NULL if not using for alleleFrequency slot. Default NULL.
#' @param center Logical, whether the variables should be shifted to be zero centered. See [stats::prcomp()].
#' @param scale. Logical, whether the variables should be scaled to have unit variance before the analysis takes place. See [stats::prcomp()].
#'
#' @return TapestriExperiment with PCA results saved to reducedDims slot of altExp, and proportion of variance saved to metadata slot of altExp.
#' @export
#'
#' @examples
#' \dontrun{TapestriExperiment <- runPCA(TapestriExperiment,
#' alt.exp = "alleleFrequency", sd.min.threshold = 35)}
runPCA <- function(TapestriExperiment, alt.exp = "alleleFrequency", assay = NULL, sd.min.threshold = NULL, center = TRUE, scale. = TRUE){

    assay <- .SelectAssay(TapestriExperiment, alt.exp = alt.exp, assay = assay)

    message(paste("Running PCA on:", alt.exp, assay))

    # filter by sd.min.threshold
    if(!is.null(sd.min.threshold)){
        feature.set.filter <- SingleCellExperiment::rowData(SingleCellExperiment::altExp(TapestriExperiment, alt.exp))$allelefreq.sd >= sd.min.threshold
        if(all(!feature.set.filter)){
            stop("All features filtered out. Reduce threshold.")
        }
    }

    pca.assay <- getTidyData(TapestriExperiment, alt.exp = alt.exp, assay = assay, feature.id.as.factor = F)
    pca.assay <- pca.assay %>% tidyr::pivot_wider(id_cols = "feature.id", names_from = "cell.barcode", values_from = {{assay}}) %>% tibble::column_to_rownames("feature.id") %>%
        as.data.frame()

    if(!is.null(sd.min.threshold)){
        pca.assay <- pca.assay[feature.set.filter,]
    }

    #transpose matrix
    pca.assay <- t(pca.assay)

    # run PCA
    pca.result <- stats::prcomp(pca.assay, center = center, scale. = scale.)

    # save PCA to main experiment
    if(is.null(alt.exp)){
        SingleCellExperiment::reducedDim(TapestriExperiment, "PCA", withDimnames = T) <- pca.result$x
        S4Vectors::metadata(TapestriExperiment)$pca.proportion.of.variance <- summary(pca.result)$importance["Proportion of Variance", ]
        S4Vectors::metadata(TapestriExperiment)$pca.assay <- assay
    } else { # save PCA to alt experiment
        SingleCellExperiment::reducedDim(SingleCellExperiment::altExp(TapestriExperiment, alt.exp), "PCA", withDimnames = T) <- pca.result$x
        S4Vectors::metadata(SingleCellExperiment::altExp(TapestriExperiment, alt.exp))$pca.proportion.of.variance <- summary(pca.result)$importance["Proportion of Variance", ]
        S4Vectors::metadata(SingleCellExperiment::altExp(TapestriExperiment, alt.exp))$pca.assay <- assay
    }

    return(TapestriExperiment)
}

#' Plot PCA Variance Knee
#'
#' Draws knee plot of PCA variance explained to determine which PCs to include for downstream applications e.g. clustering.
#'
#' @param TapestriExperiment TapestriExperiment object
#' @param alt.exp Chr string indicating altExp to use, NULL uses top-level experiment. Default "alleleFrequency".
#' @param n.pcs Numeric vector, indicating number of PCs to plot, starting at 1. Default 10.
#'
#' @return ggplot knee plot
#' @export
#'
#' @examples
#' \dontrun{PCAKneePlot(TapestriExperiment, pcs = 1:5)}
PCAKneePlot <- function(TapestriExperiment, alt.exp = "alleleFrequency", n.pcs = 10){

    if(is.null(alt.exp)){
        knee.df <- data.frame("prop.variance" = S4Vectors::metadata(TapestriExperiment)$pca.proportion.of.variance)
        chart.title <- S4Vectors::metadata(TapestriExperiment)$pca.assay
    } else {
        knee.df <- data.frame("prop.variance" = S4Vectors::metadata(SingleCellExperiment::altExp(TapestriExperiment, alt.exp))$pca.proportion.of.variance)
        chart.title <- S4Vectors::metadata(SingleCellExperiment::altExp(TapestriExperiment, alt.exp))$pca.assay
    }

    knee.df$cumulative.prop <- cumsum(knee.df$prop.variance)
    knee.df$index <- seq_len(nrow(knee.df))

    if(n.pcs > nrow(knee.df)){
        warning("n.pcs exceeds number of principal components.")
        n.pcs <- nrow(knee.df)
    }

    knee.df <- knee.df[seq_len(n.pcs), ]
    knee.df$index <- as.factor(knee.df$index)

    g1 <- ggplot(knee.df) +
        geom_line(aes(x = .data$index, y = .data$prop.variance, group = 1)) +
        geom_point(aes(x = .data$index, y = .data$prop.variance)) +
        geom_col(aes(x = .data$index, y = .data$cumulative.prop), alpha = 0.3, width = 0.5) +
        theme_bw() +
        coord_cartesian(ylim = c(0,1)) +
        theme(panel.grid.minor.x = element_blank()) +
        labs(x = "Principal Component", y = "Percent/cumulative Variance Explained", title = "Variance Explained by Principal Components", subtitle = chart.title)

    return(g1)
}

#' Cluster Data by UMAP
#'
#' @param TapestriExperiment TapestriExperiment object
#' @param alt.exp Chr string indicating altExp to use, NULL uses top-level experiment. Default "alleleFrequency".
#' @param assay Chr string indicating assay to use. NULL (default) selects first listed assay. Npt used when `use.pca.dims = TRUE`.
#' @param use.pca.dims Logical, if TRUE, uses experiment PCA, otherwise uses assay data. Default TRUE.
#' @param pca.dims Numeric vector, indicating indices of PCs to use in UMAP. NULL (default) uses all dimensions.
#' @param ... Additional parameters to pass to umap, e.g. for configuration (see [`umap::umap.defaults`]).
#'
#' @return TapestriExperiment with UMAP embeddings saved to reducedDims slot of altExp.
#' @export
#'
#' @examples
#' \dontrun{TapestriExperiment <- runUMAP(TapestriExperiment, input.dims = 1:3)}
runUMAP <- function(TapestriExperiment, alt.exp = "alleleFrequency", assay = NULL, use.pca.dims = TRUE, pca.dims = NULL, ...){

    assay <- .SelectAssay(TapestriExperiment, alt.exp = alt.exp, assay = assay)

    message(paste("Running UMAP on:", alt.exp, assay))

    # PCA
    if(use.pca.dims){

        if(is.null(alt.exp)){
            umap.assay <- SingleCellExperiment::reducedDim(TapestriExperiment, type = "PCA")
        } else {
            umap.assay <- SingleCellExperiment::reducedDim(SingleCellExperiment::altExp(TapestriExperiment, alt.exp), type = "PCA")
        }

        if(is.null(pca.dims)){
            warning("pca.dims not set. Running on all PCs.")
        } else {
            umap.assay <- umap.assay[,pca.dims]
        }
    } else {

        if(is.null(alt.exp)){
            umap.assay <- SummarizedExperiment::assay(TapestriExperiment, assay)
        } else
            umap.assay <- SummarizedExperiment::assay(SingleCellExperiment::altExp(TapestriExperiment, alt.exp))
    }

    # run UMAP
    umap.result <- umap::umap(umap.assay, ... = ...)
    umap.embeddings <- as.data.frame(umap.result$layout)
    colnames(umap.embeddings) <- c("umap.1", "umap.2")

    # save UMAP to main experiment
    if(is.null(alt.exp)){
        SingleCellExperiment::reducedDim(TapestriExperiment, "UMAP", withDimnames = T) <- umap.embeddings
        S4Vectors::metadata(TapestriExperiment)$"umap.assay" <- ifelse(use.pca.dims, "PCA", assay)
    } else { # save UMAP to alt experiment
        SingleCellExperiment::reducedDim(altExp(TapestriExperiment, alt.exp), "UMAP", withDimnames = T) <- umap.embeddings
        S4Vectors::metadata(SingleCellExperiment::altExp(TapestriExperiment, alt.exp))$"umap.assay" <- ifelse(use.pca.dims, "PCA", assay)
    }

    return(TapestriExperiment)
}

#' Scatter plot for dimensional reduction results.
#'
#' @param TapestriExperiment TapestriExperiment object
#' @param dim.reduction Chr, which dimension reduction to plot, either "PCA" or "UMAP".
#' @param group.label Chr string indicating colData column for coloring samples. Default NULL.
#' @param alt.exp Chr string indicating altExp to use, NULL uses top-level experiment. Default "alleleFrequency".
#' @param dim.x Numeric, index of dimensional reduction data to plot on X axis. Default 1.
#' @param dim.y Numeric, index of dimensional reduction data to plot on Y axis. Default 2.
#'
#' @return ggplot scatter plot
#' @export
#'
#' @examples
#' \dontrun{reducedDimPlot(TapestriExperiment, dim.reduction = "pca")}
reducedDimPlot <- function(TapestriExperiment, alt.exp =  "alleleFrequency",
                           dim.reduction, dim.x = 1, dim.y = 2, group.label = NULL){

    dim.reduction <- toupper(dim.reduction)

    if(is.null(alt.exp)){
        to.plot <- reducedDim(TapestriExperiment, dim.reduction)
        to.plot <- to.plot[,c(dim.x, dim.y)]
    } else {
        to.plot <- reducedDim(altExp(TapestriExperiment, alt.exp), dim.reduction)
        to.plot <- to.plot[,c(dim.x, dim.y)]
    }

    if(is.null(group.label)){
        group.label.data <- NULL
    } else {
        #check for group label
        if(!group.label %in% colnames(colData(TapestriExperiment))){
            stop("group.label not found in colData.")
        }
        group.label.data <- colData(TapestriExperiment)[,group.label]
    }

    g1 <- simpleScatterPlot(x = to.plot[,1],
                            y = to.plot[,2],
                            labs.title = dim.reduction,
                            labs.x = colnames(to.plot)[1],
                            labs.y = colnames(to.plot[2]),
                            group.label = group.label.data,
                            group.label.legend = group.label)

    return(g1)

}

#' Cluster data
#'
#' Cluster data using dbscan method. Uses UMAP reduced dimensions to partition data into clusters and saves the clusters to colData.
#'
#' @param TapestriExperiment TapestriExperiment object
#' @param feature.set Chr string identifying the altExp feature set to use Default "alleleFrequency".
#' @param dim.reduction Chr string indicating the reduced dimension set to use. Default "UMAP".
#' @param eps Numeric, dbscan eps parameter. Change to adjust cluster granularity. See [`dbscan::dbscan()`]. Default 0.8.
#' @param ... Additional parameters to pass to `dbscan::dbscan()`.
#'
#' @return TapestriExperiment object with updated colData
#' @export
#'
#' @seealso [`dbscan::dbscan()`].
#'
#' @examples
#' \dontrun{TapestriExperiment <- getClusters(TapestriExperiment, dim.reduction = "UMAP", eps = 0.8)}
getClusters <- function(TapestriExperiment, feature.set = "alleleFrequency", dim.reduction = "UMAP", eps = 0.8, ...){

    dim.reduction <- toupper(dim.reduction)

    if(dim.reduction != "UMAP"){
        stop("dim.reduction currently only supports UMAP.")
    }

    if(!feature.set %in% altExpNames(TapestriExperiment)){
        stop("feature.set not found in alternative experiments.")
    }
    message(paste0("Finding clusters with ", dim.reduction, " using Alt Experiment: ", feature.set))

    dbscan.assay <- SingleCellExperiment::reducedDim(SingleCellExperiment::altExp(TapestriExperiment, feature.set), "UMAP")
    dbscan.result <- dbscan::dbscan(dbscan.assay, eps = eps, ...)
    dbscan.result.clusters <- data.frame(cell.barcode = rownames(dbscan.assay), cluster = as.factor(dbscan.result$cluster))

    dbscan.result.clusters$cluster <- .FactorNumericalOrder(dbscan.result.clusters$cluster)

    ##
    # get and merge colData in main
    existing.cell.data <- as.data.frame(SummarizedExperiment::colData(TapestriExperiment))

    # drop existing clusters if they exist to allow overwriting
    existing.cell.data <- existing.cell.data[,which(colnames(existing.cell.data) != "cluster"), drop = FALSE]

    # merge result and existing colData
    updated.cell.data <- merge(existing.cell.data, dbscan.result.clusters, by = "cell.barcode", all.x = TRUE, sort = FALSE)

    # reorder to match colData
    rownames(updated.cell.data) <- updated.cell.data$cell.barcode
    updated.cell.data <- updated.cell.data[rownames(SingleCellExperiment::colData(TapestriExperiment)),]

    # update TapestriExperiment
    SummarizedExperiment::colData(TapestriExperiment) <- S4Vectors::DataFrame(updated.cell.data)

    ##
    # get and merge colData in altExp
    existing.cell.data <- as.data.frame(SummarizedExperiment::colData(SingleCellExperiment::altExp(TapestriExperiment, feature.set)))

    # drop existing clusters if they exist to allow overwriting
    existing.cell.data <- existing.cell.data[,which(colnames(existing.cell.data) != "cluster"), drop = FALSE]

    # merge result and existing colData
    updated.cell.data <- merge(existing.cell.data, dbscan.result.clusters, by = "cell.barcode", all.x = TRUE, sort = FALSE)

    # reorder to match colData
    rownames(updated.cell.data) <- updated.cell.data$cell.barcode
    updated.cell.data <- updated.cell.data[rownames(SummarizedExperiment::colData(SingleCellExperiment::altExp(TapestriExperiment, feature.set))),]

    # update TapestriExperiment
    SummarizedExperiment::colData(SingleCellExperiment::altExp(TapestriExperiment, feature.set)) <- S4Vectors::DataFrame(updated.cell.data)

    return(TapestriExperiment)
}

.FactorNumericalOrder <- function(f){

    f <- stats::reorder(f,f,FUN=length)
    f <- factor(f, levels = rev(levels(f)))
    levels(f) <- seq_len(length(levels(f)))

    return(f)
}
