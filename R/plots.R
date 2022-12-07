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
#' @return ggplot
#' @noRd
#'
#' @import ggplot2
#'
#' @examples
#' \dontrun{simpleLinePlot(x, y)}
simpleLinePlot <- function(x, y, labs.x, labs.y, labs.title, xlim, ylim){

input.data <- data.frame(x = x, y = y)

ggplot2::ggplot(data = input.data, aes(x = x, y = y)) +
    geom_point() +
    geom_line(aes(group = 1)) +
    theme_bw() +
    labs(x = labs.x, y = labs.y, title = labs.title) +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    scale_x_continuous(breaks = xlim[1]:xlim[2]) +
    theme(panel.grid.minor.x = element_blank())
}
