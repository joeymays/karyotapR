#' Clustering by PCA
#'
#' Runs a set of features through PCA and saves results to reducedDims slot of input TapestriObject.
#'
#' @param TapestriExperiment TapestriExperiment object
#' @param feature.set Chr string identifying the altExp feature set to perform PCA on. Default "alleleFrequency".
#' @param sd.min.threshold Numeric, minimum threshold for allelefreq.sd. Increase to run PCA on fewer, increasingly variable dimensions. Default 0.
#' @param center Logical, whether the variables should be shifted to be zero centered. See [stats::prcomp()].
#' @param scale. Logical, whether the variables should be scaled to have unit variance before the analysis takes place. See [stats::prcomp()].
#'
#' @return TapestriExperiment with PCA results saved to reducedDims slot of altExp, and proportion of variance saved to metadata slot of altExp.
#' @export
#'
#' @examples
#' \dontrun{TapestriExperiment <- runPCA(TapestriExperiment, sd.min.threshold = 35)}
runPCA <- function(TapestriExperiment, feature.set = "alleleFrequency", sd.min.threshold = 0, center = T, scale. = T){

    if(!feature.set %in% altExpNames(TapestriExperiment)){
        stop("feature.set not found in alternative experiments.")
    }
    message(paste0("Running PCA on Alt Experiment: ", feature.set))

    # filter by sd.min.threshold
    feature.set.filter <- SingleCellExperiment::rowData(SingleCellExperiment::altExp(TapestriExperiment, "alleleFrequency"))$allelefreq.sd >= sd.min.threshold

    if(all(!feature.set.filter)){
        stop("All features filtered out. Reduce threshold.")
    }

    pca.assay <- SummarizedExperiment::assay(SingleCellExperiment::altExp(TapestriExperiment, "alleleFrequency"), feature.set, withDimnames = T)
    pca.assay <- pca.assay[feature.set.filter,] #subset data
    pca.assay <- t(pca.assay) #transpose matrix

    # run PCA
    pca.result <- stats::prcomp(pca.assay, center = center, scale. = scale.)

    # save PCA to object
    SingleCellExperiment::reducedDim(SingleCellExperiment::altExp(TapestriExperiment, "alleleFrequency"), "PCA", withDimnames = T) <- pca.result$x
    S4Vectors::metadata(SingleCellExperiment::altExp(TapestriExperiment, "alleleFrequency"))$pca.proportion.of.variance <- summary(pca.result)$importance["Proportion of Variance", ]

    return(TapestriExperiment)
}

#' Plot PCA Variance Knee
#'
#' Draws knee plot of PCA variance explained to determine which PCs to include for downstream applications e.g. clustering.
#'
#' @param TapestriExperiment TapestriExperiment object
#' @param feature.set Chr string identifying the altExp feature set PCA was performed on. Default "alleleFrequency".
#' @param pcs Numeric vector, indicating indices of PCs to plot. Default 1:10.
#'
#' @return ggplot knee plot
#' @export
#'
#' @examples
#' \dontrun{PCAKneePlot(TapestriExperiment, pcs = 1:5)}
PCAKneePlot <- function(TapestriExperiment, feature.set = "alleleFrequency", pcs = 1:10){

    prop.vector <- S4Vectors::metadata(SingleCellExperiment::altExp(TapestriExperiment, feature.set))$pca.proportion.of.variance
    ylim <- c(0,1)

    simpleLinePlot(x = 1:length(prop.vector), y = prop.vector, labs.title = "Variance Explained by Principal Components",
                   labs.x = "Principal Component", labs.y = "Percent Variance Explained", xlim = c(min(pcs), max(pcs)), ylim = ylim)

}

#' Cluster Data by UMAP
#'
#' @param TapestriExperiment TapestriExperiment object
#' @param feature.set Chr string identifying the altExp feature set to perform UMAP on. Default "alleleFrequency".
#' @param input.dims Numeric vector, indicating indices of PCs to use in UMAP.
#' @param ... Additional parameters to pass to umap, e.g. for configuration (see [`umap::umap.defaults`]).
#' @param use.pca.dims Logical, only TRUE supported. Default TRUE.
#'
#' @return TapestriExperiment with UMAP embeddings saved to reducedDims slot of altExp.
#' @export
#'
#' @examples
#' \dontrun{TapestriExperiment <- runUMAP(TapestriExperiment, input.dims = 1:3)}
runUMAP <- function(TapestriExperiment, feature.set = "alleleFrequency", use.pca.dims = T, input.dims = NULL, ...){

    if(!feature.set %in% altExpNames(TapestriExperiment)){
        stop("feature.set not found in alternative experiments.")
    }
    message(paste0("Running UMAP on Alt Experiment: ", feature.set))

    # subset dims
    if(use.pca.dims){
    umap.assay <- SingleCellExperiment::reducedDim(SingleCellExperiment::altExp(TapestriExperiment, feature.set), type = "PCA")[,input.dims]
    } else {
        stop(paste0("Set use.pca.dims to TRUE, nothing else currently supported."))
    }

    # run UMAP
    umap.result <- umap::umap(umap.assay, ... = ...)
    umap.embeddings <- as.data.frame(umap.result$layout)
    colnames(umap.embeddings) <- c("umap.1", "umap.2")

    # save UMAP to object
    SingleCellExperiment::reducedDim(altExp(TapestriExperiment, feature.set), "UMAP", withDimnames = T) <- umap.embeddings

    return(TapestriExperiment)
}

#' Scatter plot for reduced dimensions
#'
#' @param TapestriExperiment TapestriExperiment object
#' @param dim.reduction Chr, which dimension reduction to plot, either "PCA" or "UMAP".
#' @param pc.x Numeric, if `dim.reduction == "PCA"`, index of PC to plot. Default 1 for PC1.
#' @param pc.y Numeric, if `dim.reduction == "PCA"`, index of PC to plot. Default 2 for PC2.
#' @param feature.set Chr string identifying the altExp feature set to plot. Default "alleleFrequency".
#' @param group.label Chr string indicating colData column for coloring samples. Default NULL.
#'
#' @return ggplot scatter plot
#' @export
#'
#' @examples
#' \dontrun{reducedDimPlot(TapestriExperiment, dim.reduction = "pca")}
reducedDimPlot <- function(TapestriExperiment, dim.reduction, pc.x = 1, pc.y = 2, feature.set = "alleleFrequency", group.label = NULL){

    dim.reduction <- toupper(dim.reduction)

    if(!is.null(group.label)){

        if(group.label %in% colnames(colData(TapestriExperiment))){
            group.label <- colData(TapestriExperiment)[,group.label]
        } else {
            stop(paste0("group.label '", group.label, "' not found in colData."))
        }
    }

    if(dim.reduction == "PCA"){

        to.plot <- reducedDim(altExp(TapestriExperiment, feature.set), "PCA")

        g1 <- simpleScatterPlot(x = to.plot[,pc.x],
                                y = to.plot[,pc.y],
                                labs.x = paste0("PC", pc.x),
                                labs.y = paste0("PC", pc.y),
                                labs.title = "PCA")

    } else if(dim.reduction == "UMAP"){

        to.plot <- reducedDim(altExp(TapestriExperiment, feature.set), "UMAP")

        g1 <- simpleScatterPlot(x = to.plot[,1],
                                y = to.plot[,2],
                                labs.x = "UMAP1",
                                labs.y = "UMAP2",
                                labs.title = "UMAP",
                                group.label = group.label)

    } else {
        stop(paste0("dim.reduction", dim.reduction, "not found in object"))
    }

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

    dbscan.assay <- SingleCellExperiment::reducedDim(SingleCellExperiment::altExp(TapestriExperiment, "alleleFrequency"), "UMAP")
    dbscan.result <- dbscan::dbscan(dbscan.assay, eps = eps, ...)
    dbscan.result.clusters <- data.frame(cell.barcode = rownames(dbscan.assay), cluster = as.factor(dbscan.result$cluster))

    # get and merge colData in main
    cell.data <- as.data.frame(SummarizedExperiment::colData(TapestriExperiment))
    updated.cell.data <- merge(cell.data, dbscan.result.clusters, by = "cell.barcode", all.x = T, sort = F)
    rownames(updated.cell.data) <- updated.cell.data$cell.barcode
    updated.cell.data <- updated.cell.data[rownames(SingleCellExperiment::colData(TapestriExperiment)),]
    SummarizedExperiment::colData(TapestriExperiment) <- S4Vectors::DataFrame(updated.cell.data)

    # get and merge colData in altExp
    cell.data <- as.data.frame(SummarizedExperiment::colData(SingleCellExperiment::altExp(TapestriExperiment, feature.set)))
    updated.cell.data <- merge(cell.data, dbscan.result.clusters, by = "cell.barcode", all.x = T, sort = F)
    rownames(updated.cell.data) <- updated.cell.data$cell.barcode
    updated.cell.data <- updated.cell.data[rownames(SummarizedExperiment::colData(SingleCellExperiment::altExp(TapestriExperiment, feature.set))),]
    SummarizedExperiment::colData(SingleCellExperiment::altExp(TapestriExperiment, feature.set)) <- S4Vectors::DataFrame(updated.cell.data)

    return(TapestriExperiment)
}

