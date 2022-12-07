#' Clustering by PCA
#'
#' Runs a set of features through PCA and saves results to reducedDims slot of input TapestriObject.
#'
#' @param TapestriExperiment TapestriExperiment object
#' @param feature.set Chr string identifying the altExp feature set to perform PCA on. Default "alleleFrequency".
#' @param sd.min.threshold Numeric, minimum threshold for allelefreq.sd. Increase to run PCA on more variable dimensions. Default 0.
#' @param center Logical, whether the variables should be shifted to be zero centered. See [stats::prcomp()].
#' @param scale. Logical, whether the variables should be scaled to have unit variance before the analysis takes place. See [stats::prcomp()].
#'
#' @return TapestriExperiment with PCA results saved to reducedDims slot of altExp, and proportion of variance saved to metadata slot of altExp.
#' @export
#'
#' @examples
#' \dontrun{TapestriExperiment <- runPCA(TapestriExperiment, sd.min.threshold = 35)}
runPCA <- function(TapestriExperiment, feature.set = "alleleFrequency", sd.min.threshold = 0, center = T, scale. = T){

    # if feature.set not main exp, search in alt exp and switch
    # if(feature.set == "main" || feature.set == mainExpName(TapestriExperiment)){
    #     message("Running PCA on Main Experiment")
    # } else {
    #
    #     if(!feature.set %in% altExpNames(TapestriExperiment)){
    #         stop("feature.set not found in main experiment or alternative experiments.")
    #     }
    #     message(paste0("Running PCA on Alt Experiment: ", feature.set))
    #     TapestriExperiment <- SingleCellExperiment::swapAltExp(TapestriExperiment, name = feature.set)
    # }

    main.exp.name <- SingleCellExperiment::mainExpName(TapestriExperiment)

    if(!feature.set %in% altExpNames(TapestriExperiment)){
        stop("feature.set not found in alternative experiments.")
    }
    message(paste0("Running PCA on Alt Experiment: ", feature.set))
    TapestriExperiment <- SingleCellExperiment::swapAltExp(TapestriExperiment, name = feature.set)

    # filter by sd.min.threshold
    feature.set.filter <- SingleCellExperiment::rowData(TapestriExperiment)$allelefreq.sd >= sd.min.threshold

    if(all(!feature.set.filter)){
        stop("All features filtered out. Reduce threshold.")
    }

    TapestriExperiment.subset <- TapestriExperiment[feature.set.filter,]

    # run PCA
    pca.assay <- t(SummarizedExperiment::assay(TapestriExperiment.subset))
    pca.result <- stats::prcomp(pca.assay, center = center, scale. = scale.)

    # save PCA to object
    SingleCellExperiment::reducedDim(TapestriExperiment, "PCA") <- pca.result$x
    S4Vectors::metadata(TapestriExperiment)$pca.proportion.of.variance <- summary(pca.result)$importance["Proportion of Variance", ]

    # switch back to main experiment
    TapestriExperiment <- SingleCellExperiment::swapAltExp(TapestriExperiment, name = main.exp.name)

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

    main.exp.name <- SingleCellExperiment::mainExpName(TapestriExperiment)

    if(!feature.set %in% altExpNames(TapestriExperiment)){
        stop("feature.set not found in alternative experiments.")
    }
    message(paste0("Running UMAP on Alt Experiment: ", feature.set))
    TapestriExperiment <- SingleCellExperiment::swapAltExp(TapestriExperiment, name = feature.set)

    # subset dims
    if(use.pca.dims){
    umap.assay <- SingleCellExperiment::reducedDim(TapestriExperiment, type = "PCA")[,input.dims]
    } else {
        stop(paste0("Set use.pca.dims to TRUE, nothing else currently supported."))
    }

    # run UMAP
    umap.result <- umap::umap(umap.assay, ... = ...)
    umap.embeddings <- as.data.frame(umap.result$layout)
    colnames(umap.embeddings) <- c("umap.1", "umap.2")

    # save UMAP to object
    SingleCellExperiment::reducedDim(TapestriExperiment, "UMAP") <- umap.embeddings

    # switch back to main experiment
    TapestriExperiment <- SingleCellExperiment::swapAltExp(TapestriExperiment, name = main.exp.name)

    return(TapestriExperiment)
}
