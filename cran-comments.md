## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## Resubmission
This is a resubmission. In this version I have:

* Changed \dontrun{} to \donttest for the `calcGMMCopyNumber`, `plotCopyNumberGMM`, 
and `getGMMBoundaries` functions, which otherwise would run longer than 5 seconds.
  
* I also note that there are no references describing the methods in this package.