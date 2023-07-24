#' Get chromosome order from a string of chromosome/contig names
#'
#' `getChrOrder()` takes a string of chromosome or contig names and returns the indices of the string in typical chromosome order, i.e. 1 through 22, X, Y.
#' Contig names that do not match 1:22, X, or Y are sorted numerically and alphabetically (with numbers coming first), and added to the end of the order.
#' The output string can then be used to sort the input string into typical chromosome order.
#'
#' @param chr.vector Character vector of chromosome or contig names.
#'
#' @return A numerical vector of the input vectors indices in chromosome order.
#' @export
#'
#' @concept build experiment
#'
#' @examples
#' chr.order <- getChrOrder(c(1,"virus",5,"X",22,"plasmid","Y"))
#' ordered.vector <- c(1,"virus",5,"X",22,"plasmid","Y")[chr.order]

getChrOrder <- function(chr.vector){

    chr.df.input <- data.frame(index = seq_along(chr.vector), id = chr.vector)

    chr.df <- chr.df.input[which(chr.df.input$id %in% c(1:22,"X", "Y")),]

    non.chr.df <- chr.df.input[which(!chr.df.input$id %in% c(1:22,"X", "Y")),]

    chr.df.ordered <- chr.df[gtools::mixedorder(chr.df$id),]
    non.chr.df.ordered <- non.chr.df[gtools::mixedorder(non.chr.df$id),]

    chr.vector.order <- c(chr.df.ordered$index, non.chr.df.ordered$index)

    return(chr.vector.order)
}
