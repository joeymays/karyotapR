#' Simple Line Plot
#'
#' Draws a line plot with ggplot.
#'
#' @param x Numeric vector, x axis values
#' @param y Numeric vector, y axis values
#' @param labs.x Chr, x axis label
#' @param labs.y Chr, y axis label
#' @param labs.title Chr, chart title
#' @param xlim Numeric vector length = 2, min and max values of x axis
#' @param ylim Numeric vector length = 2, mix and max values of y axis
#'
#' @return ggplot line plot
#' @noRd
#'
#' @import ggplot2
#'
#' @examples
#' \dontrun{simpleLinePlot(x, y)}
simpleLinePlot <- function(x, y, labs.x ="", labs.y ="", labs.title ="", xlim, ylim){

    input.data <- data.frame(x = x, y = y)

    g1 <- ggplot2::ggplot(data = input.data, aes(x = x, y = y)) +
        geom_point() +
        geom_line(aes(group = 1)) +
        theme_bw() +
        labs(x = labs.x, y = labs.y, title = labs.title) +
        coord_cartesian(xlim = xlim, ylim = ylim) +
        scale_x_continuous(breaks = xlim[1]:xlim[2]) +
        theme(panel.grid.minor.x = element_blank())

    return(g1)
}


#' Simple Scatter Plot
#'
#' @param x Numeric vector, x axis values
#' @param y Numeric vector, y axis values
#' @param labs.x Chr, x axis label
#' @param labs.y Chr, y axis label
#' @param labs.title Chr, chart title
#'
#' @return ggplot scatter plot
#' @noRd
#'
#' @import ggplot2
#'
#' @examples
#' \dontrun{simpleScatterPlot(x, y)}
simpleScatterPlot <- function(x, y, group.label = NULL, labs.x ="", labs.y ="", labs.title ="", group.label.legend = ""){

    input.data <- data.frame(x = x, y = y)

    if(!is.null(group.label)){
        input.data$group.label <- group.label
        g1 <- ggplot2::ggplot(data = input.data, aes(x = x, y = y, color = group.label))
    } else {
        g1 <- ggplot2::ggplot(data = input.data, aes(x = x, y = y))
    }

    g2 <- g1 +
        geom_point(size = 1, alpha = 0.7) +
        labs(x = labs.x, y = labs.y, title = labs.title, color = group.label.legend) +
        theme_bw()

    return(g2)
}

#' Boxplot of Counts
#'
#' Draws boxplot of count data. Useful for visualizing altExp count data.
#'
#' @param TapestriExperiment TapestriExperiment object
#' @param alt.exp Chr string indicating altExp to use, NULL (default) uses main experiment.
#' @param assay Chr string indicating assay to use. Default "counts".
#' @param log.y Logical, if TRUE, scales counts by `log1p()`. Default TRUE.
#' @param split.features Logical, if TRUE, splits plot by rowData features.
#' @param coldata.set Chr string indicating colData column to use for X axis categories. Default NULL.
#'
#' @return ggplot object
#' @export
#'
#' @import ggplot2
#'
#' @examples
#'\dontrun{assayBoxPlot(TapestriExperiment, "chrYCounts", assay = "counts",
#'split.features = T, coldata.set = "cluster")}
assayBoxPlot <- function(TapestriExperiment, alt.exp = NULL, assay = "counts", log.y = T, split.features = F, coldata.set = NULL){

    if(is.null(alt.exp)){
        counts.to.plot <- SummarizedExperiment::assay(TapestriExperiment, assay)
    } else {
        counts.to.plot <- SummarizedExperiment::assay(SingleCellExperiment::altExp(TapestriExperiment, alt.exp), assay)
    }

    if(is.null(coldata.set)){
        to.join <- as.data.frame(SummarizedExperiment::colData(TapestriExperiment)[,c("cell.barcode"), drop = F])
    } else {
        to.join <- as.data.frame(SummarizedExperiment::colData(TapestriExperiment)[,c("cell.barcode", {{ coldata.set }})])
    }

    counts.to.plot <- counts.to.plot %>% as.data.frame() %>% tibble::rownames_to_column("probe.id") %>%
        tidyr::pivot_longer(cols = !dplyr::matches("probe.id"), values_to = "counts", names_to = "cell.barcode") %>%
        dplyr::left_join(to.join, by = "cell.barcode")

    if(log.y){
        counts.to.plot$counts <- log1p(counts.to.plot$counts)
        y.label <- "log(counts + 1)"
    } else {
        y.label <- "counts"
    }

    if(is.null(coldata.set)){
        g1 <- ggplot(counts.to.plot, aes(y = counts))
    } else {
        g1 <- ggplot(counts.to.plot, aes(x = .data[[coldata.set]], y = counts))
    }

    if(split.features){
        g1 <- g1 + geom_boxplot(aes(fill = .data$probe.id))
    } else {
        g1 <- g1 + geom_boxplot()
    }

    g1 <- g1 + labs(y = y.label, title = alt.exp, subtitle = assay) + theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.border = element_rect(colour = "black", fill=NA, size=1))

    return(g1)

}


