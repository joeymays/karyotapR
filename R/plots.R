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
#'@import ggplot2
#'
#' @examples
#' \dontrun{simpleScatterPlot(x, y)}
simpleScatterPlot <- function(x, y, group.label = NULL, labs.x ="", labs.y ="", labs.title =""){

    input.data <- data.frame(x = x, y = y)

    if(!is.null(group.label)){
        input.data$group.label <- group.label
        g1 <- ggplot2::ggplot(data = input.data, aes(x = x, y = y, color = group.label))
    } else {
        g1 <- ggplot2::ggplot(data = input.data, aes(x = x, y = y))
    }

    g2 <- g1 +
        geom_point(size = 1, alpha = 0.7) +
        labs(x = labs.x, y = labs.y, title = labs.title) +
        theme_bw()

    return(g2)
}
