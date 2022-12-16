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
#' Draws boxplot of matrix data. Useful for visualizing altExp count data.
#'
#' @param TapestriExperiment TapestriExperiment object
#' @param alt.exp Chr string indicating altExp to use, NULL (default) uses main experiment.
#' @param assay Chr string indicating assay to use. NULL (default) selects first listed assay.
#' @param log.y Logical, if TRUE, scales data by `log1p()`. Default TRUE.
#' @param split.features Logical, if TRUE, splits plot by rowData features. Default FALSE.
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
assayBoxPlot <- function(TapestriExperiment, alt.exp = NULL, assay = NULL, log.y = T, split.features = F, coldata.set = NULL){

    assay <- .SelectAssay(TapestriExperiment, alt.exp = alt.exp, assay = assay)

    tidy.data <- getTidyData(TapestriExperiment, alt.exp, assay)

    if(log.y){
        tidy.data[,assay] <- log1p(tidy.data[,assay, drop = T])
        y.label <- paste0("log(", assay, ") + 1")
    } else {
        y.label <- assay
    }

    if(is.null(coldata.set)){
        g1 <- ggplot(tidy.data, aes(y = .data[[assay]]))
    } else {
        g1 <- ggplot(tidy.data, aes(x = .data[[coldata.set]], y = .data[[assay]]))
    }

    if(split.features){
        g1 <- g1 + geom_boxplot(aes(fill = .data$feature.id))
    } else {
        g1 <- g1 + geom_boxplot()
    }

    g1 <- g1 + labs(y = y.label, title = alt.exp, subtitle = assay) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              panel.border = element_rect(colour = "black", fill=NA, size=1))

    return(g1)
}


#' Generate Assay Heatmap
#'
#' Use to create a heatmap of matrix data in a TapestriObject using the `ComplexHeatmap` package.
#' Heatmaps are generated as transposed representations of the matrix.
#'
#' `color.preset` presets:
#' ploidy: `circlize::colorRamp2(c(0,1,2,3,4,8), c('#2c7bb6','#abd9e9','#ffffff','#fdae61','#d7191c', "black"))`; Blue-white-red gradient from 0-2-4. 4 to 8+ is black.
#' ploidy.denoise: `circlize::colorRamp2(c(0,1,1.5,2,2.5,3,4,8), c('#2c7bb6','#abd9e9','#ffffff','#ffffff','#ffffff','#fdae61','#d7191c', "black"))`; Similar to ploidy, but white range is from 1.5-2.5 to reduce the appearance of noise around diplod cells.
#'
#' @param TapestriExperiment TapestriExperiment object
#' @param alt.exp Chr string indicating altExp slot to pull from. `NULL` (default) pulls from top-level/main experiment.
#' @param assay Chr string indicating assay slot to pull from. `NULL` (default) pulls from first-indexed assay (often "counts").
#' @param split.col.by Chr string indicating RowData field to split columns by, usually "chr" or "arm". Default NULL.
#' @param split.row.by Chr string indicating ColData field to split rows by, usually "cluster". Default NULL.
#' @param annotate.row.by Chr string indicating ColData field to use as annotation. Default NULL.
#' @param color.preset Chr string indicating color preset to use to color heatmap, either ploidy or ploidy.denoise (see Details). Supersedes `color.custom`. `NULL` (default) uses default `ComplexHeatmap` color.
#' @param color.custom Color mapping function given by `circlize::colorRamp2()`. `color.preset` must be `NULL`.
#' @param ... Additional parameters to pass to `[ComplexHeatmap::Heatmap()]`.
#'
#' @return A ComplexHeatmap object
#' @export
#'
#' @seealso \link[ComplexHeatmap]{Heatmap}
#'
#' @examples
#' \dontrun{assayHeatmap(TapestriExperiment, alt.exp = "smoothedPloidyByArm",
#' assay = "discretePloidy", split.row.by = "cluster")}
assayHeatmap <- function(TapestriExperiment, alt.exp = NULL, assay = NULL, split.col.by = NULL, split.row.by = NULL, annotate.row.by = NULL, color.preset = NULL, color.custom = NULL, ...){

    assay <- .SelectAssay(TapestriExperiment, alt.exp, assay) #check call validity

    tidy.data <- getTidyData(TapestriExperiment, alt.exp, assay)

    hm.matrix <- tidy.data %>% dplyr::select("feature.id", "cell.barcode", {{ assay }}) %>%
        tidyr::pivot_wider(id_cols = "feature.id", names_from = "cell.barcode", values_from = {{ assay }}) %>%
        tibble::column_to_rownames("feature.id")

    if(is.null(split.col.by)){
        show.column.names <- TRUE
        column.split <- NULL
    } else {
        show.column.names <- FALSE
        column.split <- tidy.data %>% dplyr::select("feature.id", {{split.col.by}}) %>% dplyr::distinct() %>% dplyr::pull({{split.col.by}})
    }

    if(is.null(split.row.by)){
        row.split <- NULL
    } else {
        row.split <- tidy.data %>% dplyr::select("cell.barcode", {{split.row.by}}) %>% dplyr::distinct() %>% dplyr::pull({{split.row.by}})
    }

    if(is.null(annotate.row.by)){
        row.annotation <- NULL
    } else {
        row.annotation.data <- tidy.data %>% dplyr::select("cell.barcode", {{annotate.row.by}}) %>% dplyr::distinct() %>%
            dplyr::pull({{annotate.row.by}}) %>% tibble::enframe(name = NULL, value = {{annotate.row.by}}) %>% as.data.frame()

        n.colors <- length(unique(row.annotation.data[!is.na(row.annotation.data[,1]),1]))
        color.vector <- viridisLite::viridis(n.colors + 1)[1:n.colors]
        names(color.vector) <- unique(row.annotation.data[!is.na(row.annotation.data[,1]),1])
        color.list <- list(color.vector)
        names(color.list)[1] <- annotate.row.by
        row.annotation <- ComplexHeatmap::rowAnnotation(df = row.annotation.data, col = color.list, border = T, na_col = "white",
                                                        annotation_name_side = "top", annotation_name_gp = grid::gpar(fontsize = 8))
    }

    if(is.null(color.preset)){
        if(is.null(color.custom)){
            hm.col <- NULL
        } else {
            hm.col <- color.custom
        }
    } else if(color.preset == "ploidy"){
        hm.col <- circlize::colorRamp2(c(0,1,2,3,4,8), c('#2c7bb6','#abd9e9','#ffffff','#fdae61','#d7191c', "black"))
    } else if(color.preset == "ploidy.denoise"){
        hm.col <- circlize::colorRamp2(c(0,1,1.5,2,2.5,3,4,8), c('#2c7bb6','#abd9e9','#ffffff','#ffffff','#ffffff','#fdae61','#d7191c', "black"))
    } else {
        hm.col <- color.custom
    }

    hm <- ComplexHeatmap::Heatmap(matrix = t(hm.matrix),
                                  cluster_rows = TRUE,
                                  cluster_row_slices = FALSE,
                                  show_row_names = FALSE,
                                  show_row_dend = FALSE,
                                  row_split = row.split,
                                  row_title_gp = grid::gpar(fontsize = 10),
                                  #
                                  cluster_columns = FALSE,
                                  show_column_names = show.column.names,
                                  column_names_side = "top",
                                  show_column_dend = FALSE,
                                  column_split = column.split,
                                  column_title_gp = grid::gpar(fontsize = 8),
                                  column_names_gp = grid::gpar(fontsize = 10),
                                  column_title_rot = 90,
                                  column_gap = unit(0, "mm"),
                                  #
                                  left_annotation = row.annotation,
                                  name = assay,
                                  border = TRUE,
                                  col = hm.col,
                                  ...)

    return(hm)

}
