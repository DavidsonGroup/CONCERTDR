.onAttach <- function(libname, pkgname) {
  packageStartupMessage("CONCERTDR: Drug Response Data Analysis Tools")
  packageStartupMessage("Version: ", utils::packageVersion("CONCERTDR"))
}
