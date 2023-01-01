tap.file <- "~/Documents/data-and-analyses/tapestri/data/20220331-CO293/20220331-CO293-reference20220419.dna.h5"
h5.filename <- tap.file

bam.file <- "~/Documents/data-and-analyses/tapestri/data/20220331-CO293/20220331-CO293-reference20220419.tube1.cells.bam"

barcode.lookup <- data.frame(ids = c("g7x1", "gNC", "g9p21L2", "g9p21R2"),
                             sequences = c("TGGATATATGGACCGCATTG", "ACGGAGGCTAAGCGTCGCAA", "TATTTACAGGGACAATACCG", "CGGTAGAATAAGCTGTACCG"))

#full set
barcode.lookup <- data.frame(ids = c("g7x1", "gNC", "g9p21L2", "g9p21R2", "g7x3"),
                             sequences = c("TGGATATATGGACCGCATTG", "ACGGAGGCTAAGCGTCGCAA", "TATTTACAGGGACAATACCG", "CGGTAGAATAAGCTGTACCG", "ACTCTTGCTGTGGCATTTTC"))

#Troubleshooting object, gNC is random sequence
barcode.lookup <- data.frame(ids = c("g7x1", "gNC", "g9p21L2", "g9p21R2"),
                             sequences = c("TGGATATATGGACCGCATTG", "CGACTGTTCCCAAATTGTAA", "TATTTACAGGGACAATACCG", "CGGTAGAATAAGCTGTACCG"))

#Troubleshooting object, all bcs are random sequences
barcode.lookup <- data.frame(ids = c("g7x1", "gNC", "g9p21L2", "g9p21R2"),
           sequences = c("GATGTAACCATATACTTACG", "CTGGATCTTCTCCCGCGAAT", "TTTAACCCTCACCAACTACG", "AGATTTGAGGTAAACCAAAT"))

parsed.bc <- parseBarcodedReads(tap, bam.file = bam.file, barcode.lookup = barcode.lookup, probe = "GRNA")

parseBarcodedReads(tap, bam.file = bam.file, barcode.lookup = barcode.lookup, probe = "sample.barcode")


bam.matches <- bam.filter$tag[sequence.matches[[x]]]

coldata.labels <- c("g7x1", "gNC")

library(devtools)
load_all()
tap.file <- "~/Documents/data-and-analyses/tapestri/data/20220331-CO293/20220331-CO293-reference20220419.dna.h5"

tap <- createTapestriExperiment(tap.file, "CO293")
tap <- parseBarcodedReads(tap, bam.file = bam.file, barcode.lookup = barcode.lookup, probe.tag = "grna")
colData(tap)$sample.grna <- callSampleLables(tap, coldata.labels = c(coldata.labels, "g9p21L2", "g9p21R2"))
tap1 <- runPCA(tap, sd.min.threshold = 35)
PCAKneePlot(tap1)
tap2 <- runUMAP(tap1, input.dims = 1:3)
reducedDimPlot(tap2, "umap", group.label = "sample.grna")
umap.dbscan <- dbscan::dbscan(reducedDim(altExp(tap2, "alleleFrequency"), "UMAP"), eps = 0.8)
colData(tap2)$cluster <- as.factor(umap.dbscan$cluster)
reducedDimPlot(tap2, "umap", group.label = "cluster")
reducedDimPlot(tap2, "umap", group.label = "sample.grna")


library(tidyverse)
load_all()
tap.file <- "~/Documents/data-and-analyses/tapestri/data/Exp003-20220808-CO293/run20220808-02-panelCO293-ref20220419.dna.h5"
bam.file <- "~/Documents/data-and-analyses/tapestri/data/Exp003-20220808-CO293/run20220808-02-panelCO293-ref20220419.tube1.cells.bam"
barcode.lookup <- data.frame(ids = c("g7x1", "gNC", "g9p21L2", "g9p21R2", "g7x3"), sequences = c("TGGATATATGGACCGCATTG", "ACGGAGGCTAAGCGTCGCAA", "TATTTACAGGGACAATACCG", "CGGTAGAATAAGCTGTACCG", "ACTCTTGCTGTGGCATTTTC"))
label.features <- c("g7x1", "gNC", "g9p21L2", "g9p21R2", "g7x3")

exp3 <- createTapestriExperiment(tap.file, "CO293")
exp3 <- parseBarcodedReads(exp3, bam.file = bam.file, barcode.lookup = barcode.lookup, probe.tag = "grna")
exp3 <- callSampleLables(exp3, label.features = c("g7x3", "gNC"), sample.label = "sample.grna")
exp3 <- runPCA(exp3, sd.min.threshold = 40)
PCAKneePlot(exp3)
exp3 <- runUMAP(exp3, pca.dims = 1:2)
reducedDimPlot(exp3, dim.reduction = "umap")
exp3 <- getClusters(exp3, eps = 0.9)
reducedDimPlot(exp3, dim.reduction = "umap", group.label = "cluster")
reducedDimPlot(exp3, dim.reduction = "umap", group.label = "sample.grna")
fct_count(colData(exp3)$cluster)
exp3.subset <- exp3[,colData(exp3)$cluster %in% 1:3]
reducedDimPlot(exp3.subset, dim.reduction = "umap", group.label = "cluster")
colData(exp3.subset)$cluster <- fct_recode(colData(exp3.subset)$cluster, RPE1 = "1", cellline1 = "2", cellline2 = "3")
colData(exp3.subset)$cluster <- fct_drop(colData(exp3.subset)$cluster)
fct_count(colData(exp3.subset)$cluster)
summary(rowData(exp3.subset)$median.reads > 0)
exp3.subset <- exp3.subset[rowData(exp3.subset)$median.reads > 0,] #1 probe filtered for exp 3
exp3.subset <- normalizeCounts(exp3.subset)
control.ploidy <- generateControlPloidyTemplate(sample.label.all = "RPE1")
control.ploidy["chr10q", "ploidy"] <- 3
exp3.subset <- getPloidy(exp3.subset, control.ploidy = control.ploidy, sample.category = "cluster")
exp3.subset <- smoothPloidy(exp3.subset)
assayHeatmap(exp3.subset, assay = "ploidy", split.col.by = "chr", split.row.by = "cluster", annotate.row.by = "sample.grna", color.preset = "ploidy")
assayHeatmap(exp3.subset, alt.exp = "smoothedPloidyByChrom", assay = "discretePloidy", split.row.by = "cluster", annotate.row.by = "sample.grna", color.preset = "ploidy")
assayBoxPlot(exp3.subset, alt.exp = "chrYCounts", split.features = T, split.x.by = "cluster")

getTidyData(exp3.subset)
getTidyData(exp3.subset, alt.exp = "alleleFrequency")
getTidyData(exp3.subset, alt.exp = "chrYCounts")
getTidyData(exp3.subset, alt.exp = "smoothedPloidyByChrom")


# readme test

example.exp <- newTapestriExperimentExample()
assayHeatmap(example.exp, split.col.by = "arm", split.row.by = "test.cluster", annotate.row.by = "test.cluster")

example.exp <- runPCA(example.exp)
PCAKneePlot(example.exp)
example.exp <- runUMAP(example.exp, pca.dims = 1:2)
reducedDimPlot(example.exp, dim.reduction = "umap")
example.exp <- getClusters(example.exp, eps = 0.9)
reducedDimPlot(example.exp, dim.reduction = "umap", group.label = "cluster")
colData(example.exp)$cluster <- forcats::fct_recode(colData(example.exp)$cluster, cellline1 = "1", cellline2 = "2", cellline3 = "3")

example.exp <- normalizeCounts(example.exp)
control.ploidy <- generateControlPloidyTemplate(example.exp, sample.label.all = "cellline3", ploidy.all = 2)
example.exp <- getPloidy(example.exp, control.ploidy = control.ploidy, sample.category = "cluster")
example.exp <- smoothPloidy(example.exp)

assayHeatmap(example.exp, assay = "ploidy", split.col.by = "arm", split.row.by = "test.cluster", annotate.row.by = "test.cluster", color.preset = "ploidy")
assayHeatmap(example.exp, alt.exp = "smoothedPloidyByArm", assay = "discretePloidy", split.row.by = "test.cluster", annotate.row.by = "test.cluster", color.preset = "ploidy")
assayBoxPlot(example.exp, alt.exp = "chrYCounts", split.features = T, split.x.by = "test.cluster")



