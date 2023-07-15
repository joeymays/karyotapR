
calcGMMCopyNumber <- function(TapestriExperiment, cell.barcodes, control.copy.number, model.components = 1:5, model.priors = NULL, ...){

    #control.copy.number <- generateControlCopyNumberTemplate(TapestriExperiment, copy.number = 2)
    #control.copy.number["chr10q","copy.number"] <- 3

    if(is.null(model.priors)){
        model.priors <- rep(1, length(model.components))
    } else {
        if(length(model.priors) != length(model.components)){
            cli::cli_abort("model.priors must be same length as model.components, or `NULL` for equal priors.")
        }
    }

    filtered.tapestri.exp <- TapestriExperiment[,cell.barcodes]

    #simulate probe counts
    simulated.norm.counts <- .generateSimulatedCNVCells(TapestriExperiment = filtered.tapestri.exp, control.copy.number = control.copy.number, ...)

    #fit Gaussian distributions to simulated cells
    cli::cli_progress_step("Fitting Gaussian distributions to simulated cells...")
    cn.model.params <- .fitGaussianDistributions(normalized.counts = simulated.norm.counts, probe.metadata = rowData(TapestriExperiment), ...)

    #calculate posterior probabilities for each data point under each model component
    cli::cli_progress_step("Calculating posterior probabilities...")
    cn.model.table <- .calcClassPosteriors(TapestriExperiment = TapestriExperiment, cn.model.params = cn.model.params, model.components = model.components, model.priors = model.priors)

    #call copy number values from posterior probabilities
    cli::cli_progress_step("Calling copy number from posterior probabilities...")
    cn.model.table <- .callCopyNumberClasses(cn.model.table)

    #transform copy number calls to matrix
    cli::cli_progress_step("Saving copy number calls to altExp: smoothedCopyNumberByChr, assay: gmmCopyNumber...")
    class.labels.df <- cn.model.table %>% dplyr::pull("cn.class") %>%
        purrr::map(\(x) tidyr::pivot_wider(x, names_from = "cell.barcode", values_from = "cn.class")) %>% purrr::list_rbind() %>%
        as.data.frame()
    rownames(class.labels.df) <- cn.model.table$feature.id

    #add copy number calls and model metadata to TapestriExperiment
    SummarizedExperiment::assay(altExp(TapestriExperiment, "smoothedCopyNumberByChr"), "gmmCopyNumber") <- class.labels.df
    #assay(altExp(TapestriExperiment, "smoothedCopyNumberByArm"), "gmmCopyNumber") <- class.labels.df

    S4Vectors::metadata(TapestriExperiment)$gmmParametersByChr <- cn.model.table
    #metadata(TapestriExperiment)$gmmParametersByArm <- cn.model.table

    return(TapestriExperiment)
}

#model probes and generate simulated cells with normalized counts
.generateSimulatedCNVCells <- function(TapestriExperiment, control.copy.number, n.simulated.cells = 500){

    cli::cli_progress_step("Simulating probes for {n.simulated.cells} cells...")

    raw.counts <- SummarizedExperiment::assay(TapestriExperiment, "counts")
    norm.counts <-  sweep(raw.counts, 2, colMeans(raw.counts), "/") * 100 # set cell means to 100
    norm.counts[norm.counts == 0] <- 1 #pseudocount zeros to ones
    norm.counts <- as.list(as.data.frame(t(norm.counts))) #convert to list of probes

    #fit probe parameters
    probe.model.fit <- suppressWarnings(norm.counts %>% purrr::map(\(x) fitdistrplus::fitdist(data = x, distr = "weibull")))
    probe.model.fit <- probe.model.fit %>% purrr::map(~.x$estimate) %>% purrr::list_transpose() %>% as.data.frame() %>%
        tibble::rownames_to_column("probe.id")

    #get control copy number, calculate scalar values for copy number = 1:6,
    #multiply scale parameter by scalar value
    probe.transform <- SingleCellExperiment::rowData(TapestriExperiment) %>% tibble::as_tibble() %>%
        dplyr::select("probe.id", "arm") %>% dplyr::inner_join(control.copy.number %>% dplyr::select(!c("sample.label")), by = "arm")
    probe.transform <- probe.transform %>% dplyr::mutate("scalar" = purrr::map(.data$copy.number, \(x) 1:6/x))
    probe.transform <- dplyr::inner_join(probe.transform, probe.model.fit, by = "probe.id")
    probe.transform <- probe.transform %>% dplyr::mutate("scale.transform" = purrr::map2(.data$scalar, .data$scale, \(x, y) x * y))
    probe.transform$probe.id <- factor(probe.transform$probe.id, levels = unique(probe.transform$probe.id))

    #simulate normalized counts for each parameter set
    simulated.counts <- probe.transform %>% tidyr::unnest(.data$scale.transform) %>%
        dplyr::mutate("sim.counts" = purrr::map2(.data$shape, .data$scale.transform, function(shape, scale.transform){
            stats::rweibull(n=n.simulated.cells, shape = shape, scale = scale.transform)
        }))

    #combine counts into matrix
    simulated.counts <- simulated.counts %>% dplyr::select(.data$probe.id, .data$sim.counts)
    simulated.counts <- split(simulated.counts, simulated.counts$probe.id) %>%
        purrr::map(\(x) unlist(x$sim.counts)) %>% as.data.frame() %>% t()
    colnames(simulated.counts) <- paste0("sim_cn", 1:6) %>% purrr::map(\(x) paste(x, 1:n.simulated.cells, sep = "_")) %>% unlist()

    return(simulated.counts)
}

#fit Gaussian distributions to simulated cells
.fitGaussianDistributions <- function(normalized.counts, probe.metadata, genome = "hg19"){

    #generate tapestri experiment object and copy normalized counts normcounts slot
    tapestri.sim <- .createTapestriExperiment.sim(counts = normalized.counts, probe.metadata = probe.metadata, genome = genome)
    SummarizedExperiment::assay(tapestri.sim, "normcounts") <- normalized.counts

    #delimit cell barcodes to get copy number classes
    cn.sim.class <- stringr::str_split_fixed(colData(tapestri.sim)$cell.barcode, pattern = "_", n = 3)
    cn.sim.class <- stringr::str_c(cn.sim.class[,1], cn.sim.class[,2], sep = "_")
    SingleCellExperiment::colData(tapestri.sim)$cn.sim.class <- as.factor(cn.sim.class)

    control.cn <- generateControlCopyNumberTemplate(tapestri.sim, copy.number = 2, sample.feature.label = "sim_cn2")
    tapestri.sim <- calcCopyNumber(tapestri.sim, control.copy.number = control.cn, sample.feature = "cn.sim.class")
    tapestri.sim <- suppressMessages(calcSmoothCopyNumber(tapestri.sim))

    sim.data.tidy <- getTidyData(tapestri.sim, alt.exp = "smoothedCopyNumberByChr", assay = "smoothedCopyNumber")

    #fit Gaussian distributions; group_by makes cluster and feature columns persist when using transmute
    cn.model.params <- sim.data.tidy %>% tidyr::nest(.by = c("feature.id", "cn.sim.class"), .key = "smoothed.counts") %>%
      dplyr::mutate("fit" = purrr::map(.data$smoothed.counts, \(x) suppressWarnings(fitdistrplus::fitdist(x$smoothedCopyNumber, distr = "norm")), .progress = T)) %>%
      dplyr::mutate("params" = purrr::map(.data$fit, \(x) x$estimate)) %>% tidyr::unnest_wider(.data$params) %>% dplyr::select("feature.id", "cn.sim.class", "mean", "sd")

    return(cn.model.params)
}

#get probabilities for belonging to a given component of the Gaussian mixture
.calcClassPosteriors <- function(TapestriExperiment, cn.model.params, model.components, model.priors){

    components.filtered <- paste0("sim_cn", model.components)

    #get smoothed copy number and combine in tibble with copy number model parameters
    smoothed.cn.df <- getTidyData(TapestriExperiment, alt.exp = "smoothedCopyNumberByChr", assay = "smoothedCopyNumber") %>%
        dplyr::select("feature.id", "cell.barcode", "smoothedCopyNumber") %>% tidyr::nest(.by = "feature.id", .key = "smoothed.cn")
    cn.params.df <- cn.model.params %>% dplyr::filter(.data$cn.sim.class %in% components.filtered) %>% tidyr::nest(.by = "feature.id", .key = "model")
    gmm.table <- dplyr::inner_join(smoothed.cn.df, cn.params.df, by = "feature.id")

    #iterate over both lists in parallel to get probability density function values across all GMM components for all data points
    gmm.table <- gmm.table %>% dplyr::mutate("pdf" = purrr::map2(.data$smoothed.cn, .data$model, function(smoothed.cn, model){
        df <- purrr::pmap(model, function(cn.sim.class, mean, sd){stats::setNames(data.frame(stats::dnorm(smoothed.cn[["smoothedCopyNumber"]], mean = mean, sd = sd)), cn.sim.class) }) %>%
            purrr::list_cbind()
        rownames(df) <- smoothed.cn[["cell.barcode"]]
        df <- tibble::as_tibble(df)
        return(df)
    }))

    # get Bayes theorem denominator for each GMM (aka evidence or marginal likelihood)
    gmm.table <- gmm.table %>% dplyr::mutate("model.evidence" = purrr::map(.data$pdf, function(pdf){
        as.matrix(pdf) %*% model.priors %>% drop()
    }))

    #calculate probability of a data point belonging to each copy number class (i.e. copy number class given a data point)
    gmm.table <- gmm.table %>% dplyr::mutate("cn.probability" = purrr::map2(.data$pdf, .data$model.evidence, function(pdf, model.evidence){
        probs <- sweep(pdf, 2, model.priors, '*') #priors multiplied through each column
        probs <- sweep(probs, 1, model.evidence, "/")  #model.evidence divided through each row
        probs <- tibble::as_tibble(probs)
        return(probs)
    }))

    return(gmm.table)
}

#assign copy numbers to data points by highest posterior probability
.callCopyNumberClasses <- function(cn.model.table){

    #iterate over cn.probability and smoothed.cn lists in parallel
    cn.model.table <- cn.model.table %>%
        dplyr::mutate(cn.class = purrr::map2(.data$cn.probability, .data$smoothed.cn, function(cn.probability, smoothed.cn){
            result <- max.col(cn.probability, "first") #pick class with highest probability
            result <- colnames(cn.probability)[result] #label selections
            result <- as.numeric(stringr::str_split_fixed(result, pattern = "cn", n = 2)[,2]) #pull copy number value from class label
            result <- tibble::as_tibble(data.frame(cell.barcode = purrr::pluck(smoothed.cn, "cell.barcode"), cn.class = result))
        }))

    return(cn.model.table)
}































