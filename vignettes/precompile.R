# Precompute the vignettes whose models take too long to fit on every
# R CMD check.
#
# For each name below, vignettes/<name>.Rmd.orig is the real source. Knitting
# it here runs every chunk once and writes vignettes/<name>.Rmd with the code
# and its output already in it, plus the figures as vignettes/<name>-*.png.
# R CMD check then only renders that file to HTML; no model is fitted.
#
# The output is only as current as the last knit. Re-run this on the code
# being pushed before every push, and commit whatever it changes:
#
#   Rscript vignettes/precompile.R           # every precomputed vignette
#   Rscript vignettes/precompile.R rilta     # just one
#
# from the package root. Needs pandoc only for the later render, not here.

precomputed <- c("lta", "mglca_yrbs", "rilta")

knit_precomputed <- function(name, src_dir = "vignettes", out_dir = src_dir) {
  src <- normalizePath(file.path(src_dir, paste0(name, ".Rmd.orig")), mustWork = TRUE)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  old <- setwd(out_dir)
  on.exit(setwd(old), add = TRUE)
  unlink(Sys.glob(paste0(name, "-*.png")))
  knitr::opts_chunk$set(fig.path = paste0(name, "-"))
  knitr::knit(src, output = paste0(name, ".Rmd"), envir = new.env(), quiet = TRUE)
  invisible(file.path(out_dir, paste0(name, ".Rmd")))
}

if (sys.nframe() == 0L) {
  pkgload::load_all(".", quiet = TRUE)
  args <- commandArgs(trailingOnly = TRUE)
  for (name in if (length(args)) args else precomputed) {
    t0 <- Sys.time()
    knit_precomputed(name)
    cat(sprintf("%s: %.1f min\n", name,
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
}
