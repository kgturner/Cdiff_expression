pdinfo <- "pd.110405.cdiffusa.lz.exp"

# the command
# eg. library(pd.mapping50k.hind240)
# will have been run in the harness

## test_crlmm <- function() {
## # this stuff will work for the 4 chips we know about now
## # 100k hind240, xba240; 250k sty, nsp
## size = as.character(2*as.numeric(gsub(".*g(.*)k.*", "\\1", pdinfo)))
## enz = gsub(".*k.(.*$)", "\\1", pdinfo)
## enz = gsub("240", "", enz)
## hapmPackName = paste("hapmap", size, "k", enz, sep="")
## library(oligo)
## library(hapmPackName, character.only=TRUE)
## xxr = justSNPRMA(dir(system.file( "celFiles", package=hapmPackName), full=TRUE))
## xxc = crlmm(xxr, correctionFile="corr.rda")
## if (exists("xxc")) return(TRUE)
## return(FALSE)
## }

test_crlmm <- function(){
  pkgname <- switch(pdinfo,
                    pd.mapping50k.xba240 = "hapmap100kxba",
                    pd.mapping50k.hind240 = "hapmap100khind",
                    pd.mapping250k.nsp = "hapmap500knsp",
                    pd.mapping250k.sty = "hapmap500ksty",
                    pd.genomewidesnp.5 = "hapmapsnp5",
                    pd.genomewidesnp.6 = "hapmapsnp6")
  library(oligo)
  library(pkgname, character.only=TRUE)
  outdir <- file.path(getwd(), paste(pkgname, "test", sep="_"))
  celfiles <- list.celfiles(system.file("celFiles", package=pkgname), full=TRUE)
  result <- crlmm(celfiles, outdir)
  unlink(outdir, recursive=TRUE)
  return(result)
}
