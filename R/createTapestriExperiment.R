#' Create TapestriExperiment object from Tapestri Pipeline output
#'
#' @param h5.filename file path for .h5 file from Tapestri Pipeline output.
#' @param panel.id Tapestri panel name, CO261 and CO293 supported only. Default NULL.
#'
#' @return TapestriExperiment
#' @export
#'
#' @import SingleCellExperiment
#'
#' @examples
#' \dontrun{x <- createTapestriExperiment("myh5file.h5", "CO293")}
createTapestriExperiment <- function(h5.filename, panel.id = NULL){

    # read panel ID
    if(is.null(panel.id)){
        barcodeProbe = "not specified"
        grnaProbe = "not specified"
    } else if(panel.id == "CO293"){
        barcodeProbe = "AMPL205334"
        grnaProbe = "AMPL205666"
    } else if(panel.id == "CO261") {
        barcodeProbe = "not specified"
        grnaProbe = "not specified"
    } else {
        stop(paste("panel.id", panel.id, "is not recognized. Please specify CO261 or CO293, or set speciality probes manually."))
    }

    #construct object
    tapestri.h5 <- rhdf5::H5Fopen(file.path(h5.filename))

    read.counts.raw <- t(matrix(data = tapestri.h5$"/assays/dna_read_counts/layers/read_counts",
                                ncol=length(tapestri.h5$'/assays/dna_read_counts/ca/id'),
                                byrow = T))

    read.counts.raw.colData <- S4Vectors::DataFrame(cell.barcode = tapestri.h5$'/assays/dna_read_counts/ra/barcode',
                                                    row.names = tapestri.h5$'/assays/dna_read_counts/ra/barcode')

    read.counts.raw.rowData <- S4Vectors::DataFrame(probe.id = tapestri.h5$'/assays/dna_read_counts/ca/id',
                                                    chr = tapestri.h5$'/assays/dna_read_counts/ca/CHROM',
                                                    start.pos = tapestri.h5$'/assays/dna_read_counts/ca/start_pos',
                                                    end.pos = tapestri.h5$'/assays/dna_read_counts/ca/end_pos',
                                                    row.names = tapestri.h5$'/assays/dna_read_counts/ca/id')

    chr.order <- getChrOrder(read.counts.raw.rowData$chr) #reorder amplicon metadata by chromosome order

    read.counts.raw.rowData <- read.counts.raw.rowData[chr.order,]

    read.counts.raw <- read.counts.raw[chr.order,]

    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = read.counts.raw),
                                                      colData = read.counts.raw.colData,
                                                      rowData = read.counts.raw.rowData)

    tapestri.object <- .TapestriExperiment(sce)

    SingleCellExperiment::mainExpName(tapestri.object) <- "gwCNV"

    # apply panel ID probe shortcuts
    tapestri.object@barcodeProbe = barcodeProbe
    tapestri.object@grnaProbe = grnaProbe

    #variant allele frequency data
    variant.metadata <- data.frame(
        filtered = as.logical(tapestri.h5$'/assays/dna_variants/ca/filtered'),
        chr = as.character(tapestri.h5$'/assays/dna_variants/ca/CHROM'),
        position = as.numeric(tapestri.h5$'/assays/dna_variants/ca/POS'),
        quality = as.numeric(tapestri.h5$'/assays/dna_variants/ca/QUAL'),
        amplicon.id = as.character(tapestri.h5$'/assays/dna_variants/ca/amplicon'),
        variant.id = as.character(tapestri.h5$'/assays/dna_variants/ca/id'),
        ado.rate = as.numeric(tapestri.h5$'/assays/dna_variants/ca/ado_rate'),
        reference.allele = as.character(tapestri.h5$'/assays/dna_variants/ca/REF'),
        alternate.allele = as.character(tapestri.h5$'/assays/dna_variants/ca/ALT'),
        ado.gt.cells = as.numeric(tapestri.h5$'/assays/dna_variants/ca/ado_gt_cells'))

    variant.metadata$chr <- factor(variant.metadata$chr, unique(variant.metadata$chr))

    af.matrix <- tapestri.h5$'/assays/dna_variants/layers/AF'

    af.sd <- apply(af.matrix, 1, stats::sd)

    allele.frequency <- SummarizedExperiment::SummarizedExperiment(list(alleleFrequency = af.matrix),
                                                                   rowData = S4Vectors::DataFrame(variant.metadata,
                                                                                                  allelefreq.sd = af.sd,
                                                                                                  row.names = variant.metadata$variant.id))
    colnames(allele.frequency) <- tapestri.h5$'/assays/dna_read_counts/ra/barcode'
    SingleCellExperiment::altExp(tapestri.object, "alleleFrequency") <- allele.frequency

    #close h5
    rhdf5::H5Fclose(tapestri.h5)

    return(tapestri.object)
}
