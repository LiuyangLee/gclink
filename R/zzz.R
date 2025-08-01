# R/zzz.R

.onLoad <- function(libname, pkgname) {
  op <- options()
  op.gclink <- list(
    gclink.path = tempdir(),
    gclink.verbose = FALSE
  )
  toset <- !(names(op.gclink) %in% names(op))
  if (any(toset)) options(op.gclink[toset])
  ns <- asNamespace(pkgname)
  tryCatch({
    utils::data(list = c("blastp_df", "PGC_group", "photosynthesis_gene_list", "seq_data"),
                package = pkgname,
                envir = ns)
  }, error = function(e) {
    warning("Data load failed: ", e$message, call. = FALSE)
  })
  invisible()
}

.onAttach <- function(libname, pkgname) {
  version <- utils::packageVersion(pkgname)
  packageStartupMessage("Welcome to gclink v", version)
}
