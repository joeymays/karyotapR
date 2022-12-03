tap.file <- "~/Documents/data-and-analyses/tapestri/data/20220331-CO293/20220331-CO293-reference20220419.dna.h5"

bam.file <- "~/Documents/data-and-analyses/tapestri/data/20220331-CO293/20220331-CO293-reference20220419.tube1.cells.bam"

barcode.lookup <- data.frame(ids = c("g7x1", "gNC", "g9p21L2", "g9p21R2"),
                             sequences = c("TGGATATATGGACCGCATTG", "ACGGAGGCTAAGCGTCGCAA", "TATTTACAGGGACAATACCG", "CGGTAGAATAAGCTGTACCG"))

barcode.lookup <- data.frame(ids = c("g7x1", "gNC", "g9p21L2", "g9p21R2"),
                             sequences = c("TGGATATATGGACCGCATTG", "ACGGAGGCTAAGCGTCGCAA", "TATTTACAGGGACAATACCG", "CGGTAGAATAAGCTGTACCG"))

#Troubleshooting object, gNC is random sequence
barcode.lookup <- data.frame(ids = c("g7x1", "gNC", "g9p21L2", "g9p21R2"),
                             sequences = c("TGGATATATGGACCGCATTG", "CGACTGTTCCCAAATTGTAA", "TATTTACAGGGACAATACCG", "CGGTAGAATAAGCTGTACCG"))

#Troubleshooting object, all bcs are random sequences
data.frame(ids = c("g7x1", "gNC", "g9p21L2", "g9p21R2"),
           sequences = c("GATGTAACCATATACTTACG", "CTGGATCTTCTCCCGCGAAT", "TTTAACCCTCACCAACTACG", "AGATTTGAGGTAAACCAAAT"))

parsed.bc <- parseBarcodedReads(tap, bam.file = bam.file, barcode.lookup = barcode.lookup, probe = "GRNA")

parseBarcodedReads(tap, bam.file = bam.file, barcode.lookup = barcode.lookup, probe = "sample.barcode")


bam.matches <- bam.filter$tag[sequence.matches[[x]]]

coldata.labels <- c("g7x1", "gNC")
