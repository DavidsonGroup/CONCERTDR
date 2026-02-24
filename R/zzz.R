.onLoad <- function(libname, pkgname) {
  # Check if necessary packages are available
  required_packages <- c("BiocManager", "data.table", "cmapR")
  missing_packages <- setdiff(required_packages, installed.packages()[,"Package"])

  if (length(missing_packages) > 0) {
    missing_bioc <- intersect(missing_packages, c("cmapR"))
    missing_cran <- setdiff(missing_packages, missing_bioc)

    warning_msg <- "CONCERTDR depends on the following packages that are not installed: "

    if (length(missing_cran) > 0) {
      warning_msg <- paste0(warning_msg,
                            "\n  From CRAN: ", paste(missing_cran, collapse = ", "),
                            ". Install using: install.packages(c('",
                            paste(missing_cran, collapse = "', '"), "'))")
    }

    if (length(missing_bioc) > 0) {
      warning_msg <- paste0(warning_msg,
                            "\n  From Bioconductor: ", paste(missing_bioc, collapse = ", "),
                            ". Install using: BiocManager::install(c('",
                            paste(missing_bioc, collapse = "', '"), "'))")
    }

    warning(warning_msg)
  }
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("CONCERTDR: Drug Response Data Analysis Tools")
  packageStartupMessage("Version: ", utils::packageVersion("CONCERTDR"))

  # Add recommendation to install data.table if not present
  if (!requireNamespace("data.table", quietly = TRUE)) {
    packageStartupMessage("\nNOTE: For optimal performance and full functionality, it's recommended to install 'data.table':")
    packageStartupMessage("      install.packages('data.table')")
  }
}
