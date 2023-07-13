
calcGMMCopyNumber <- function(TapestriExperiment, cell.barcodes, model.components = 1:5, model.priors = NULL, ...){

    if(is.null(model.priors)){
        model.priors <- rep(1, length(model.components))
    }

    filtered.tapestri.exp <- TapestriExperiment[,cell.barcodes]

    #simulate probe counts
    cli::cli_progress_step("Simulating Probes...")
    simulated.norm.counts <- .generateSimulatedCNVCells(filtered.tapestri.exp, ...)




    #fitGaussianParameters()
    #calcClassPosteriors()
    #callCopyNumberClasses()

    return(simulated.norm.counts)

}

#model probes and generate simulated cells with normalized counts
.generateSimulatedCNVCells <- function(TapestriExperiment, n.simulated.cells = 500){

    raw.counts <- assay(TapestriExperiment, "counts")
    norm.counts <-  sweep(raw.counts, 2, colMeans(raw.counts), "/") * 100 # set cell means to 100
    norm.counts[norm.counts == 0] <- 1 #pseudocount zeros to ones
    norm.counts <- as.list(as.data.frame(t(norm.counts))) #convert to list of probes

    #get probe parameters
    probe.model.fit <- suppressWarnings(norm.counts %>% map(\(x) fitdist(data = x, distr = "weibull")))
    probe.model.fit <- probe.model.fit %>% map(\(x) x$estimate) %>% as.data.frame()

    #generate scaled scale parameters for Copy Number = 1:6
    scale.multiplier <- c(0.5, 1.0, 1.5, 2.0, 2.5, 3.0)
    probe.model.all.params <- list()

    for(i in seq_len(length(scale.multiplier))){
        temp.params <- probe.model.fit
        temp.params["scale", ] <- temp.params["scale", ] * scale.multiplier[i]
        probe.model.all.params[[i]] <- as.list(temp.params)
    }
    names(probe.model.all.params) <- paste0("X", scale.multiplier)

    #simulate normalized counts for each probe in each scale set
    #map_depth applies function to 2nd level in list (each probe instead of the list of probes)
    cnv.norm.counts <- probe.model.all.params %>% map_depth(.depth = 2, function(y) rweibull(n = n.simulated.cells, shape = y[1], scale = y[2]))

    #combine cell vectors into matrices and transpose
    cnv.norm.counts <- cnv.norm.counts %>% map(\(x) as.data.frame(x)) %>% list_rbind() %>% t()
    colnames(cnv.norm.counts) <- paste0("sim_cn", 1:6) %>% map(\(x) paste(x, 1:n.simulated.cells, sep = "_")) %>% unlist()

    return(cnv.norm.counts)
}
