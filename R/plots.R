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
#' \dontrun{
#' simpleLinePlot(x, y)
#' }
simpleLinePlot <- function(x, y, labs.x = "", labs.y = "", labs.title = "", xlim, ylim) {
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
#' \dontrun{
#' simpleScatterPlot(x, y)
#' }
simpleScatterPlot <- function(x, y, group.label = NULL, labs.x = "", labs.y = "", labs.title = "", group.label.legend = "") {
  input.data <- data.frame(x = x, y = y)

  if (!is.null(group.label)) {
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

#' Generate a box plot from matrix data
#'
#' Draws box plot of matrix data stored in a `TapestriExperiment`.
#' This is especially useful for visualizing `altExp` count data, such as counts from
#' probes on chrY.
#'
#' @param TapestriExperiment `TapestriExperiment` object
#' @param alt.exp Character, `altExp` to plot. `NULL` (default) uses the top-level experiment in `TapestriExperiment`.
#' @param assay Character, assay to plot. `NULL` (default) selects first assay listed `TapestriExperiment`.
#' @param log.y Logical, if `TRUE`, scales data using `log1p()`. Default `TRUE.`
#' @param split.features Logical, if `TRUE`, splits plot by `rowData` features. Default `FALSE.`
#' @param split.x.by Character, `colData` column to use for X axis categories. Default `NULL`.
#'
#' @return A ggplot object using [`ggplot2::geom_boxplot()`].
#' @export
#'
#' @import ggplot2
#'
#' @examples
#' \dontrun{
#' assayBoxPlot(TapestriExperiment, "chrYCounts",
#'   assay = "counts",
#'   split.features = T, split.x.by = "cluster"
#' )
#' }
assayBoxPlot <- function(TapestriExperiment, alt.exp = NULL, assay = NULL, log.y = TRUE, split.features = FALSE, split.x.by = NULL) {
  assay <- .SelectAssay(TapestriExperiment, alt.exp = alt.exp, assay = assay)

  tidy.data <- getTidyData(TapestriExperiment, alt.exp, assay)

  if (log.y) {
    tidy.data[, assay] <- log1p(tidy.data[, assay, drop = TRUE])
    y.label <- paste0("log(", assay, "+ 1)")
  } else {
    y.label <- assay
  }

  if (is.null(split.x.by)) {
    g1 <- ggplot(tidy.data, aes(y = .data[[assay]]))
  } else {
    g1 <- ggplot(tidy.data, aes(x = .data[[split.x.by]], y = .data[[assay]]))
  }

  if (split.features) {
    g1 <- g1 + geom_boxplot(aes(fill = .data$feature.id))
  } else {
    g1 <- g1 + geom_boxplot()
  }

  g1 <- g1 + labs(y = y.label, title = alt.exp, subtitle = assay) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.border = element_rect(colour = "black", fill = NA, size = 1)
    )

  return(g1)
}


#' Generate heatmap of matrix data
#'
#' Creates a heatmap of matrix data in a `TapestriObject` using the `ComplexHeatmap` package.
#' Heatmaps are generated as transposed (i.e. x-y flipped) representations of the indicated matrix.
#'
#' @details
#' # `color.preset` Options
#' ## "ploidy"
#' Blue-white-red gradient from 0-2-4. 4 to 8+ is red-black gradient.
#' ```
#' circlize::colorRamp2(c(0,1,2,3,4,8),
#' c('#2c7bb6','#abd9e9','#ffffff','#fdae61','#d7191c', "black"))
#' ```
#' ## "ploidy.denoise"
#' Similar to ploidy, but white range is from 1.5-2.5 to reduce the appearance of noise around diploid cells.
#' ```
#' circlize::colorRamp2(c(0,1,1.5,2,2.5,3,4,8),
#' c('#2c7bb6','#abd9e9','#ffffff','#ffffff','#ffffff','#fdae61','#d7191c', "black"))
#' ````
#'
#' @param TapestriExperiment `TapestriExperiment` object
#' @param alt.exp Character, `altExp` slot to use. `NULL` (default) uses top-level/main experiment.
#' @param assay Character, `assay` slot to use. `NULL` (default) uses first-indexed assay (usually "counts").
#' @param split.col.by Character, `rowData` column to split columns by, usually "chr" or "arm". Default `NULL`.
#' @param split.row.by Character, `colData` column to split rows by, usually "cluster". Default `NULL`.
#' @param annotate.row.by Character, `colData` column to use for annotation. Default `NULL`.
#' @param color.preset Character, color preset to use to color heatmap, either "ploidy" or "ploidy.denoise" (see `Details`). Overrides `color.custom`. `NULL` (default) uses default `ComplexHeatmap` coloring.
#' @param color.custom Color mapping function given by [`circlize::colorRamp2()`]. `color.preset` must be `NULL`.
#' @param ... Additional parameters to pass to [`ComplexHeatmap::Heatmap()`].
#'
#' @return A `ComplexHeatmap` object
#' @export
#'
#' @seealso \link[ComplexHeatmap]{Heatmap}
#'
#' @examples
#' \dontrun{
#' assayHeatmap(TapestriExperiment,
#'   alt.exp = "smoothedCopyNumberByArm",
#'   assay = "discreteCopyNumber", split.row.by = "cluster"
#' )
#' }
assayHeatmap <- function(TapestriExperiment, alt.exp = NULL, assay = NULL, split.col.by = NULL, split.row.by = NULL, annotate.row.by = NULL, color.preset = NULL, color.custom = NULL, ...) {
  assay <- .SelectAssay(TapestriExperiment, alt.exp, assay) # check call validity

  tidy.data <- getTidyData(TapestriExperiment, alt.exp, assay)

  hm.matrix <- tidy.data %>%
    dplyr::select("feature.id", "cell.barcode", {{ assay }}) %>%
    tidyr::pivot_wider(id_cols = "feature.id", names_from = "cell.barcode", values_from = {{ assay }}) %>%
    tibble::column_to_rownames("feature.id")

  if (is.null(split.col.by)) {
    show.column.names <- TRUE
    column.split <- NULL
  } else {
    show.column.names <- FALSE
    column.split <- tidy.data %>%
      dplyr::select("feature.id", {{ split.col.by }}) %>%
      dplyr::distinct() %>%
      dplyr::pull({{ split.col.by }})
  }

  if (is.null(split.row.by)) {
    row.split <- NULL
  } else {
    row.split <- tidy.data %>%
      dplyr::select("cell.barcode", {{ split.row.by }}) %>%
      dplyr::distinct() %>%
      dplyr::pull({{ split.row.by }})
  }

  if (is.null(annotate.row.by)) {
    row.annotation <- NULL
  } else {
    row.annotation.data <- tidy.data %>%
      dplyr::select("cell.barcode", {{ annotate.row.by }}) %>%
      dplyr::distinct() %>%
      dplyr::pull({{ annotate.row.by }}) %>%
      tibble::enframe(name = NULL, value = {{ annotate.row.by }}) %>%
      as.data.frame()

    n.colors <- length(unique(row.annotation.data[!is.na(row.annotation.data[, 1]), 1]))
    color.vector <- viridisLite::viridis(n.colors + 1)[1:n.colors]
    names(color.vector) <- unique(row.annotation.data[!is.na(row.annotation.data[, 1]), 1])
    color.list <- list(color.vector)
    names(color.list)[1] <- annotate.row.by
    row.annotation <- ComplexHeatmap::rowAnnotation(
      df = row.annotation.data, col = color.list, border = TRUE, na_col = "white",
      annotation_name_side = "top", annotation_name_gp = grid::gpar(fontsize = 8)
    )
  }

  if (is.null(color.preset)) {
    if (is.null(color.custom)) {
      hm.col <- NULL
    } else {
      hm.col <- color.custom
    }
  } else if (color.preset == "ploidy") {
    hm.col <- circlize::colorRamp2(c(0, 1, 2, 3, 4, 8), c("#2c7bb6", "#abd9e9", "#ffffff", "#fdae61", "#d7191c", "black"))
  } else if (color.preset == "ploidy.denoise") {
    hm.col <- circlize::colorRamp2(c(0, 1, 1.5, 2, 2.5, 3, 4, 8), c("#2c7bb6", "#abd9e9", "#ffffff", "#ffffff", "#ffffff", "#fdae61", "#d7191c", "black"))
  } else {
    hm.col <- color.custom
  }

  # set default params here to allow overwriting in function call
  hm.defaults <- list("name" = assay,
                      "hm.col" = hm.col,
                      "left.annotation" = row.annotation,
                      "row.split" = row.split,
                      "show.column.names" = show.column.names,
                      "column.split" = column.split)

  hm <- .ComplexHeatmap.default(matrix = t(hm.matrix), hm.defaults = hm.defaults, ...)

  return(hm)
}

# Internal ComplexHeatmap call with reasonable default settings
.ComplexHeatmap.default <- function(matrix,
                                    cluster_rows = TRUE,
                                    cluster_row_slices = FALSE,
                                    show_row_names = FALSE,
                                    show_row_dend = FALSE,
                                    row_split = hm.defaults[["row.split"]],
                                    row_title_gp = grid::gpar(fontsize = 10),
                                    cluster_columns = FALSE,
                                    show_column_names = hm.defaults[["show.column.names"]],
                                    column_names_side = "top",
                                    show_column_dend = FALSE,
                                    column_split = hm.defaults[["column.split"]],
                                    column_title_gp = grid::gpar(fontsize = 8),
                                    column_names_gp = grid::gpar(fontsize = 10),
                                    column_names_rot = 90,
                                    column_title_rot = 90,
                                    column_gap = unit(0, "mm"),
                                    left_annotation = hm.defaults[["left.annotation"]],
                                    name = hm.defaults[["name"]],
                                    border = TRUE,
                                    col = hm.defaults[["hm.col"]],
                                    hm.defaults = hm.defaults,
                                    ...
                                    ){

    complex.hm <- ComplexHeatmap::Heatmap(

        matrix = matrix,
        cluster_rows = cluster_rows,
        cluster_row_slices = cluster_row_slices,
        show_row_names = show_row_names,
        show_row_dend = show_row_dend,
        row_split = row_split,
        row_title_gp = row_title_gp,
        cluster_columns = cluster_columns,
        show_column_names = show_column_names,
        column_names_side = column_names_side,
        show_column_dend = show_column_dend,
        column_split = column_split,
        column_title_gp = column_title_gp,
        column_names_gp = column_names_gp,
        column_title_rot = column_title_rot,
        column_gap = column_gap,
        left_annotation = left_annotation,
        name = name,
        border = border,
        col = col,
        column_names_rot = column_names_rot,
        ...
        )

    return(complex.hm)
}



