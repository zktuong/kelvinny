#' documenting kelvinny shortcut
#'
#' @return quickly redocumenting kelvinny
#' @examples
#' init_Kelvinny() # leaves folder, install, and change back to folder
#' @import devtools
#' @import roxygen2
#' @export
init_Kelvinny <- function() {
roxygen2::document()
setwd('..')
roxygen2::install('kelvinny')
setwd('kelvinny')
}
