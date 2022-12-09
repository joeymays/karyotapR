tap.file <- "~/Documents/data-and-analyses/tapestri/data/20220331-CO293/20220331-CO293-reference20220419.dna.h5"
h5.filename <- tap.file

bam.file <- "~/Documents/data-and-analyses/tapestri/data/20220331-CO293/20220331-CO293-reference20220419.tube1.cells.bam"

barcode.lookup <- data.frame(ids = c("g7x1", "gNC", "g9p21L2", "g9p21R2"),
                             sequences = c("TGGATATATGGACCGCATTG", "ACGGAGGCTAAGCGTCGCAA", "TATTTACAGGGACAATACCG", "CGGTAGAATAAGCTGTACCG"))

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



tap.file <- "~/Documents/data-and-analyses/tapestri/data/Exp003-20220808-CO293/run20220808-02-panelCO293-ref20220419.dna.h5"
bam.file <- "~/Documents/data-and-analyses/tapestri/data/Exp003-20220808-CO293/run20220808-02-panelCO293-ref20220419.tube1.cells.bam"
barcode.lookup <- data.frame(ids = c("g7x1", "gNC", "g9p21L2", "g9p21R2", "g7x3"), sequences = c("TGGATATATGGACCGCATTG", "ACGGAGGCTAAGCGTCGCAA", "TATTTACAGGGACAATACCG", "CGGTAGAATAAGCTGTACCG", "ACTCTTGCTGTGGCATTTTC"))
coldata.labels <- c("g7x1", "gNC", "g9p21L2", "g9p21R2", "g7x3")

tap <- createTapestriExperiment(tap.file, "CO293")
tap <- parseBarcodedReads(tap, bam.file = bam.file, barcode.lookup = barcode.lookup, probe.tag = "grna")
tap <- callSampleLables(tap, coldata.labels = c("g7x3", "gNC"), sample.label = "sample.grna")
tap <- runPCA(tap, sd.min.threshold = 40)
PCAKneePlot(tap)
tap <- runUMAP(tap, input.dims = 1:2)
reducedDimPlot(tap, "umap")
tap <- getClusters(tap, eps = 0.9)
reducedDimPlot(tap, "umap", group.label = "cluster")
reducedDimPlot(tap, "umap", group.label = "sample.grna")


