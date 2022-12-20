#' Create TapestriExperiment object from Tapestri Pipeline output
#'
#' `createTapestriExperiment()` constructs a `TapestriExperiment` container object from data in the `.h5` file output by the Tapestri Pipeline.
#' Read count matrix (probe x cell barcode) is stored in the "counts" `assay` slot of the top-level experiment.
#' Allele frequency matrix (variant x cell barcode) is stored in the "alleleFrequency" `assay` slot of the "alleleFrequency" `altExp` (alternative experiment) slot.
#' `panel.id` is an optional shortcut to set special probe identities for specific custom panels.
#' See details for full explanation of construction.
#'
#' # Automatic Operations
#' ## Raw Data
#' Read count and allele frequency matrices are imported to their appropriate slots as described above.
#' `filter.variants == TRUE` (default) only loads allele frequency variants that have passed internal filters in the Tapestri Pipeline.
#' This greatly reduces the number of variants from tens of thousands to hundreds of likely more consequential variants,
#' saving RAM and reducing operation time.
#'
#' ## Metadata
#' Several metadata sets are copied or generated and then stored in the appropriate `TapestriExperiment` slot during construction.
#' - Probe panel metadata stored in the `.h5` file are copied to `rowData`.
#' - Basic QC stats (e.g. total number of reads per probe) are added to `rowData.`
#' - Basic QC stats (e.g. total number of reads per cell barcode) are added to `colData.`
#' - Experiment-level metadata is stored in `metadata`.
#'
#' # Optional Operations
#' Two additional major operations are called by default during `TapestriExperiment` construction for convenience.
#' `get.cytobands == TRUE` (default) calls [getCytobands()], which retrieves the chromosome arm and cytoband for each probe based on stored positional data and saves them in `rowData`.
#' Some downstream smoothing and plotting functions may fail if chromosome arms are not present in `rowData`, so this generally should always be run.
#' `move.non.genome.probes` calls [moveNonGenomeProbes()], which moves probes corresponding to the specified tags to `altExp` (alternative experiment) slots in the `TapestriExperiment` object.
#' These probes should be those which do not correspond to a chromosome and therefore would not be used to call copy number variants.
#' The exception is probes on chromosome Y; chrY does not generally experience CNVs, so we move it to an `altExp` for separate analysis.
#' Probes corresponding to the `barcodeProbe` and `grnaProbe` slots, which are specified by the `panel.id` shortcut or manually (see [Custom Slot Getters and Setters]),
#' are automatically moved to `altExp` by this operation as well.
#' If such probes are not present, the function will generate a message but not throw an error, so it is always safe to run by default.
#'
#' @param h5.filename File path for `.h5` file from Tapestri Pipeline output.
#' @param panel.id Tapestri panel name, CO261 and CO293 only supported. Acts as shortcut for setting custom slots. Default `NULL`.
#' @param get.cytobands Logical, whether to retrieve and add chromosome cytobands and chromosome arms to probe metadata. Default `TRUE`.
#' @param genome Character, reference genome to pull cytobands and arms from. Only "hg19" (default) is currently supported.
#' @param move.non.genome.probes Character vector,  non-genomic probes to move counts and metadata to `altExp` slots, if available. Default `c("grna", "sample.barcode", "Y")`.
#' @param filter.variants Logical, if `TRUE`, only loads variants that have passed Tapestri Pipeline filters. Default `TRUE`.
#'
#' @return Constructed `TapestriExperiment` object.
#' @export
#'
#' @import SingleCellExperiment
#'
#' @seealso [moveNonGenomeProbes()], [getCytobands()], which are run as part of this function by default for convenience.
#'
#' @examples
#' \dontrun{
#' tapExperiment <- createTapestriExperiment("myh5file.h5", "CO293")
#' }
createTapestriExperiment <- function(h5.filename, panel.id = NULL, get.cytobands = TRUE, genome = "hg19", move.non.genome.probes = c("grna", "sample.barcode", "Y"),
                                     filter.variants = T) {
  # read panel ID
  if (is.null(panel.id)) {
    barcodeProbe <- "not specified"
    grnaProbe <- "not specified"
  } else if (panel.id == "CO293") {
    barcodeProbe <- "AMPL205334"
    grnaProbe <- "AMPL205666"
  } else if (panel.id == "CO261") {
    barcodeProbe <- "not specified"
    grnaProbe <- "not specified"
  } else {
    stop(paste("panel.id", panel.id, "is not recognized. Please specify CO261 or CO293, or NULL to set speciality probes manually."))
  }

  # import data
  tapestri.h5 <- rhdf5::H5Fopen(file.path(h5.filename))

  # report data
  message(paste("sample name:", tapestri.h5$"/assays/dna_read_counts/metadata/sample_name"))
  message(paste("pipeline panel name:", tapestri.h5$"/assays/dna_read_counts/metadata/panel_name"))
  message(paste("pipeline version:", tapestri.h5$"/assays/dna_read_counts/metadata/pipeline_version"))
  message(paste("number of cells:", tapestri.h5$"/assays/dna_read_counts/metadata/n_cells"))
  message(paste("number of probes:", tapestri.h5$"/assays/dna_read_counts/metadata/n_amplicons"))
  message(paste("date created:", tapestri.h5$"/metadata/date_created"))

  # construct object
  read.counts.raw <- t(matrix(
    data = tapestri.h5$"/assays/dna_read_counts/layers/read_counts",
    ncol = length(tapestri.h5$"/assays/dna_read_counts/ca/id"),
    byrow = T
  ))

  read.counts.raw.colData <- S4Vectors::DataFrame(
    cell.barcode = as.character(tapestri.h5$"/assays/dna_read_counts/ra/barcode"),
    row.names = as.character(tapestri.h5$"/assays/dna_read_counts/ra/barcode")
  )

  read.counts.raw.rowData <- S4Vectors::DataFrame(
    probe.id = as.character(tapestri.h5$"/assays/dna_read_counts/ca/id"),
    chr = tapestri.h5$"/assays/dna_read_counts/ca/CHROM",
    start.pos = as.numeric(tapestri.h5$"/assays/dna_read_counts/ca/start_pos"),
    end.pos = as.numeric(tapestri.h5$"/assays/dna_read_counts/ca/end_pos"),
    row.names = as.character(tapestri.h5$"/assays/dna_read_counts/ca/id")
  )

  chr.order <- getChrOrder(read.counts.raw.rowData$chr) # reorder amplicon metadata by chromosome order

  read.counts.raw.rowData <- read.counts.raw.rowData[chr.order, ]

  read.counts.raw <- read.counts.raw[chr.order, ]

  read.counts.raw.rowData$chr <- factor(read.counts.raw.rowData$chr, levels = unique(read.counts.raw.rowData$chr))

  # basic metrics
  read.counts.raw.colData$total.reads <- colSums(read.counts.raw) # total reads per cell
  read.counts.raw.rowData$total.reads <- rowSums(read.counts.raw) # total reads per probe
  read.counts.raw.rowData$median.reads <- apply(read.counts.raw, 1, median) # median reads per probe
  mean.reads <- round(mean(colMeans(read.counts.raw)), 2) # mean reads/cell/probe(or amplicon)
  message(paste("mean reads per cell per probe:", mean.reads))

  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = read.counts.raw),
    colData = read.counts.raw.colData,
    rowData = read.counts.raw.rowData
  )

  tapestri.object <- .TapestriExperiment(sce)

  SingleCellExperiment::mainExpName(tapestri.object) <- "CNV"

  # save experiment metadata
  S4Vectors::metadata(tapestri.object)$sample.name <- as.character(tapestri.h5$"/assays/dna_read_counts/metadata/sample_name")
  S4Vectors::metadata(tapestri.object)$pipeline.panel.name <- as.character(tapestri.h5$"/assays/dna_read_counts/metadata/panel_name")
  S4Vectors::metadata(tapestri.object)$pipeline.version <- as.character(tapestri.h5$"/assays/dna_read_counts/metadata/pipeline_version")
  S4Vectors::metadata(tapestri.object)$number.of.cells <- as.character(tapestri.h5$"/assays/dna_read_counts/metadata/n_cells")
  S4Vectors::metadata(tapestri.object)$number.of.probes <- as.character(tapestri.h5$"/assays/dna_read_counts/metadata/n_amplicons")
  S4Vectors::metadata(tapestri.object)$date.h5.created <- as.character(tapestri.h5$"/metadata/date_created")
  S4Vectors::metadata(tapestri.object)$mean.reads.per.cell.per.probe <- as.character(mean.reads)

  # apply panel ID probe shortcuts
  tapestri.object@barcodeProbe <- barcodeProbe
  tapestri.object@grnaProbe <- grnaProbe

  # if keeping only variants that passed thresholds
  if (filter.variants == T) {
    # variant filtering; filtered == TRUE in h5 object means "filtered out".
    filtered.variants <- as.logical(tapestri.h5$"/assays/dna_variants/ca/filtered")
    filtered.variants <- !filtered.variants

    # variant allele frequency data
    variant.metadata <- data.frame(
      filtered = as.logical(tapestri.h5$"/assays/dna_variants/ca/filtered")[filtered.variants],
      chr = as.character(tapestri.h5$"/assays/dna_variants/ca/CHROM")[filtered.variants],
      position = as.numeric(tapestri.h5$"/assays/dna_variants/ca/POS")[filtered.variants],
      quality = as.numeric(tapestri.h5$"/assays/dna_variants/ca/QUAL")[filtered.variants],
      amplicon.id = as.character(tapestri.h5$"/assays/dna_variants/ca/amplicon")[filtered.variants],
      variant.id = as.character(tapestri.h5$"/assays/dna_variants/ca/id")[filtered.variants],
      ado.rate = as.numeric(tapestri.h5$"/assays/dna_variants/ca/ado_rate")[filtered.variants],
      reference.allele = as.character(tapestri.h5$"/assays/dna_variants/ca/REF")[filtered.variants],
      alternate.allele = as.character(tapestri.h5$"/assays/dna_variants/ca/ALT")[filtered.variants],
      ado.gt.cells = as.numeric(tapestri.h5$"/assays/dna_variants/ca/ado_gt_cells")[filtered.variants]
    )

    af.matrix <- tapestri.h5$"/assays/dna_variants/layers/AF"[filtered.variants, ]
    dimnames(af.matrix) <- list(variant.metadata$variant.id, tapestri.h5$"/assays/dna_read_counts/ra/barcode")
  } else {
    # variant allele frequency data, no filtering, keep all variants
    variant.metadata <- data.frame(
      filtered = as.logical(tapestri.h5$"/assays/dna_variants/ca/filtered"),
      chr = as.character(tapestri.h5$"/assays/dna_variants/ca/CHROM"),
      position = as.numeric(tapestri.h5$"/assays/dna_variants/ca/POS"),
      quality = as.numeric(tapestri.h5$"/assays/dna_variants/ca/QUAL"),
      amplicon.id = as.character(tapestri.h5$"/assays/dna_variants/ca/amplicon"),
      variant.id = as.character(tapestri.h5$"/assays/dna_variants/ca/id"),
      ado.rate = as.numeric(tapestri.h5$"/assays/dna_variants/ca/ado_rate"),
      reference.allele = as.character(tapestri.h5$"/assays/dna_variants/ca/REF"),
      alternate.allele = as.character(tapestri.h5$"/assays/dna_variants/ca/ALT"),
      ado.gt.cells = as.numeric(tapestri.h5$"/assays/dna_variants/ca/ado_gt_cells")
    )

    af.matrix <- tapestri.h5$"/assays/dna_variants/layers/AF"
    dimnames(af.matrix) <- list(variant.metadata$variant.id, tapestri.h5$"/assays/dna_read_counts/ra/barcode")
  }

  variant.metadata$chr <- factor(variant.metadata$chr, unique(variant.metadata$chr))
  variant.metadata$allelefreq.sd <- apply(af.matrix, 1, stats::sd)
  rownames(variant.metadata) <- variant.metadata$variant.id

  allele.frequency <- SingleCellExperiment::SingleCellExperiment(list(alleleFrequency = af.matrix),
    rowData = S4Vectors::DataFrame(variant.metadata),
    colData = S4Vectors::DataFrame(cell.barcode = colnames(af.matrix))
  )

  allele.frequency <- .TapestriExperiment(allele.frequency)
  allele.frequency@barcodeProbe <- barcodeProbe
  allele.frequency@grnaProbe <- grnaProbe

  SingleCellExperiment::altExp(tapestri.object, "alleleFrequency", withDimnames = T) <- allele.frequency

  # close h5
  rhdf5::H5Fclose(tapestri.h5)

  # get cytobands
  if (get.cytobands) {
    tapestri.object <- getCytobands(tapestri.object)
  }

  # move non-genomic probes to altExp slots
  if (!identical(move.non.genome.probes, F)) {
    tapestri.object <- moveNonGenomeProbes(tapestri.object, move.non.genome.probes)
  }

  show(tapestri.object)

  return(tapestri.object)
}
