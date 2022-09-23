#' Create TapestriExperiment object from Tapestri Pipeline output
#'
#' @param h5.filename file path for .h5 file form Tapestri Pipeline output
#' @param panel.id Tapestri panel name, CO261 and CO293 supported only.
#'
#' @return TapestriExperiment
#' @export
#'
#' @import rhdf5
#' @import S4Vectors
#' @importFrom SingleCellExperiment SingleCellExperiment
#'
#' @examples
#' \dontrun{x <- CreateTapestriExperiment("myh5file.h5", "CO293")}
createTapestriExperiment <- function(h5.filename, panel.id = ""){

    if(!panel.id %in% c("CO261", "CO293")){
        stop(paste("panel.id", panel.id, "is not recognized. Please use CO261 or CO293."))
    }

    tapestri.h5 <- rhdf5::H5Fopen(file.path(h5.filename))

    read.counts.raw <- t(matrix(data = tapestri.h5$"/assays/dna_read_counts/layers/read_counts",
                                ncol=length(tapestri.h5$'/assays/dna_read_counts/ca/id'),
                                byrow = T))

    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = read.counts.raw),
                                                      colData = S4Vectors::DataFrame(cell.barcode = tapestri.h5$'/assays/dna_read_counts/ra/barcode',
                                                                                     row.names = tapestri.h5$'/assays/dna_read_counts/ra/barcode'),
                                                      rowData = S4Vectors::DataFrame(probe.id = tapestri.h5$'/assays/dna_read_counts/ca/id',
                                                                                     chr = tapestri.h5$'/assays/dna_read_counts/ca/CHROM',
                                                                                     start.pos = tapestri.h5$'/assays/dna_read_counts/ca/start_pos',
                                                                                     end.pos = tapestri.h5$'/assays/dna_read_counts/ca/end_pos',
                                                                                     row.names = tapestri.h5$'/assays/dna_read_counts/ca/id'))

    tapestri.object <- .TapestriExperiment(sce)

    SingleCellExperiment::mainExpName(tapestri.object) <- "gwCNV"

    if(panel.id == "CO293"){
        tapestri.object@barcodeProbe = "AMPL205334"
        tapestri.object@grnaProbe = "AMPL205666"
    } else if(panel.id == "CO261") {
        tapestri.object@barcodeProbe = "CO261_AMP1"
        tapestri.object@grnaProbe = NULL
    }

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

    af.sd <- apply(af.matrix, 1, sd)

    allele.frequency <- SummarizedExperiment::SummarizedExperiment(list(alleleFrequency = af.matrix),
                                                                   rowData = S4Vectors::DataFrame(variant.metadata,
                                                                                                  allelefreq.sd = af.sd,
                                                                                                  row.names = variant.metadata$variant.id))
    colnames(allele.frequency) <- tapestri.h5$'/assays/dna_read_counts/ra/barcode'
    SingleCellExperiment::altExp(tapestri.object, "alleleFrequency") <- allele.frequency


    #close
    rhdf5::H5Fclose(tapestri.h5)

    return(tapestri.object)
}
