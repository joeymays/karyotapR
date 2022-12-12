corner <- function(input.mat){
    if(nrow(input.mat) > 4){
        row.out = 5
    } else { row.out = nrow(input.mat)}

    if(ncol(input.mat) > 7){
        col.out = 8
    } else { col.out = ncol(input.mat)}

    print(input.mat[1:row.out, 1:col.out])
}

getTidyPloidy <- function(TapestriExperiment){

    ploidy.tidy <- assay(TapestriExperiment, "ploidy")

    ploidy.tidy <- ploidy.tidy %>% as.data.frame() %>% rownames_to_column("probe.id") %>% as_tibble() %>%
        pivot_longer(cols = !matches("probe.id"), names_to = "cell.barcode", values_to = "ploidy")

    ploidy.tidy <- ploidy.tidy %>% left_join(as.data.frame(colData(TapestriExperiment)), by = "cell.barcode")

    ploidy.tidy <- ploidy.tidy %>% left_join(as.data.frame(rowData(TapestriExperiment)), by = "probe.id")

    return(ploidy.tidy)
}
