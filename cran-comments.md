## R CMD check results

0 errors | 0 warnings | 1 note

* Namespace in Imports field not imported from: 'GenomeInfoDb'
    All declared Imports should be used.

The NOTE regarding GenomeInfoDb reflects that it is a transitive dependency 
required by Bioconductor packages used internally (e.g. GenomicRanges).
It is explicitly listed in Imports to ensure correct installation order on 
systems where Bioconductor dependencies may be incomplete 
(e.g. R-devel during Bioconductor transitions). Such a transition was the 
source of an error for automated CRAN checks with devel environments.