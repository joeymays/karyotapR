
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
    class.labels.df <- cn.model.table %>% pull(cn.class) %>%
        map(\(x) pivot_wider(x, names_from = cell.barcode, values_from = cn.class)) %>% list_rbind() %>%
        as.data.frame()
    rownames(class.labels.df) <- cn.model.table$feature.id

    #add copy number calls and model metadata to TapestriExperiment
    assay(altExp(TapestriExperiment, "smoothedCopyNumberByChr"), "gmmCopyNumber") <- class.labels.df
    #assay(altExp(TapestriExperiment, "smoothedCopyNumberByArm"), "gmmCopyNumber") <- class.labels.df

    metadata(TapestriExperiment)$gmmParametersByChr <- cn.model.table
    #metadata(TapestriExperiment)$gmmParametersByArm <- cn.model.table

    return(TapestriExperiment)
}

#model probes and generate simulated cells with normalized counts
.generateSimulatedCNVCells <- function(TapestriExperiment, control.copy.number, n.simulated.cells = 500){

    cli::cli_progress_step("Simulating probes for {n.simulated.cells} cells...")

    raw.counts <- assay(TapestriExperiment, "counts")
    norm.counts <-  sweep(raw.counts, 2, colMeans(raw.counts), "/") * 100 # set cell means to 100
    norm.counts[norm.counts == 0] <- 1 #pseudocount zeros to ones
    norm.counts <- as.list(as.data.frame(t(norm.counts))) #convert to list of probes

    #fit probe parameters
    probe.model.fit <- suppressWarnings(norm.data %>% map(~fitdist(data = .x, distr = probe.dist)))
    probe.model.fit <- probe.model.fit %>% map(~.x$estimate) %>% list_transpose() %>% as.data.frame() %>%
        rownames_to_column("probe.id")

    #get control copy number, calculate scalar values for copy number = 1:6,
    #multiply scale parameter by scalar value
    probe.transform <- rowData(tapestri.experiment) %>% as_tibble() %>%
        select(probe.id, arm) %>% inner_join(control.copy.number %>% select(!sample.label), by = "arm")
    probe.transform <- probe.transform %>% mutate(scalar = map(copy.number, \(x) 1:6/x))
    probe.transform <- inner_join(probe.transform, probe.model.fit, by = "probe.id")
    probe.transform <- probe.transform %>% mutate(scale.transform = map2(scalar, scale, \(x, y) x * y))
    probe.transform$probe.id <- factor(probe.transform$probe.id, levels = unique(probe.transform$probe.id))

    #simulate normalized counts for each parameter set
    simulated.counts <- probe.transform %>% unnest(scale.transform) %>%
        mutate(sim.counts = map2(shape, scale.transform, function(shape, scale.transform){
            rweibull(n=n.simulated.cells, shape = shape, scale = scale.transform)
        }))

    #combine counts into matrix
    simulated.counts <- simulated.counts %>% select(probe.id, sim.counts) %>% split(., .$probe.id) %>%
        map(\(x) unlist(x$sim.counts)) %>% as.data.frame() %>% t()
    colnames(simulated.counts) <- paste0("sim_cn", 1:6) %>% map(\(x) paste(x, 1:n.simulated.cells, sep = "_")) %>% unlist()

    return(simulated.counts)
}

#fit Gaussian distributions to simulated cells
.fitGaussianDistributions <- function(normalized.counts, probe.metadata, genome = "hg19"){

    #generate tapestri experiment object and copy normalized counts normcounts slot
    tapestri.sim <- CNweaveR:::.createTapestriExperiment.sim(counts = normalized.counts, probe.metadata = probe.metadata, genome = genome)
    assay(tapestri.sim, "normcounts") <- normalized.counts

    #delimit cell barcodes to get copy number classes
    cn.sim.class <- str_split_fixed(colData(tap.sim)$cell.barcode, pattern = "_", n = 3)
    cn.sim.class <- str_c(cn.sim.class[,1], cn.sim.class[,2], sep = "_")
    colData(tapestri.sim)$cn.sim.class <- as.factor(cn.sim.class)

    control.cn <- generateControlCopyNumberTemplate(tapestri.sim, copy.number = 2, sample.feature.label = "sim_cn2")
    tapestri.sim <- calcCopyNumber(tapestri.sim, control.copy.number = control.cn, sample.feature = "cn.sim.class")
    tapestri.sim <- suppressMessages(calcSmoothCopyNumber(tapestri.sim))

    sim.data.tidy <- getTidyData(tapestri.sim, alt.exp = "smoothedCopyNumberByChr", assay = "smoothedCopyNumber")

    #fit Gaussian distributions; group_by makes cluster and feature columns persist when using transmute
    cn.model.params <- sim.data.tidy %>% nest(.by = c(feature.id, cn.sim.class))%>%
        mutate(fit = map(data, \(x) suppressWarnings(fitdist(x$smoothedCopyNumber, "norm")), .progress = T)) %>%
        mutate(params = map(fit, \(x) x$estimate)) %>% unnest_wider(params) %>% select(feature.id, cn.sim.class, mean, sd)

    return(cn.model.params)
}

#get probabilities for belonging to a given component of the Gaussian mixture
.calcClassPosteriors <- function(TapestriExperiment, cn.model.params, model.components, model.priors){

    components.filtered <- paste0("sim_cn", model.components)

    #get smoothed copy number and combine in tibble with copy number model parameters
    smoothed.cn.df <- getTidyData(TapestriExperiment, alt.exp = "smoothedCopyNumberByChr", assay = "smoothedCopyNumber") %>%
        dplyr::select(feature.id, cell.barcode, smoothedCopyNumber) %>% nest(.by = feature.id, .key = "smoothed.cn")
    cn.params.df <- cn.model.params %>% filter(cn.sim.class %in% components.filtered) %>% nest(.by = feature.id, .key = "model")
    gmm.table <- inner_join(smoothed.cn.df, cn.params.df, by = "feature.id")

    #iterate over both lists in parallel to get probability density function values across all GMM components for all data points
    gmm.table <- gmm.table %>% mutate(pdf = map2(smoothed.cn, model, function(smoothed.cn, model){
        df <- pmap(model, function(cn.sim.class, mean, sd){setNames(data.frame(dnorm(smoothed.cn[["smoothedCopyNumber"]], mean = mean, sd = sd)), cn.sim.class) }) %>%
            list_cbind()
        rownames(df) <- smoothed.cn[["cell.barcode"]]
        df <- as_tibble(df)
        return(df)
    }))

    # get Bayes theorem denominator for each GMM (aka evidence or marginal likelihood)
    gmm.table <- gmm.table %>% mutate(model.evidence = map(pdf, function(pdf){
        as.matrix(pdf) %*% model.priors %>% drop()
    }))

    #calculate probability of a data point belonging to each copy number class (i.e. copy number class given a data point)
    gmm.table <- gmm.table %>% mutate(cn.probability = map2(pdf, model.evidence, function(pdf, model.evidence){
        probs <- sweep(pdf, 2, model.priors, '*') #priors multiplied through each column
        probs <- sweep(probs, 1, model.evidence, "/")  #model.evidence divided through each row
        probs <- as_tibble(probs)
        return(probs)
    }))

    return(gmm.table)
}

#assign copy numbers to data points by highest posterior probability
.callCopyNumberClasses <- function(cn.model.table){

    #iterate over cn.probability and smoothed.cn lists in parallel
    cn.model.table <- cn.model.table %>%
        mutate(cn.class = map2(cn.probability, smoothed.cn, function(cn.probability, smoothed.cn){
            result <- max.col(cn.probability, "first") #pick class with highest probability
            result <- colnames(cn.probability)[result] #label selections
            result <- as.numeric(str_split_fixed(result, pattern = "cn", n = 2)[,2]) #pull copy number value from class label
            result <- as_tibble(data.frame(cell.barcode = pluck(smoothed.cn, "cell.barcode"), cn.class = result))
        }))

    return(cn.model.table)
}































