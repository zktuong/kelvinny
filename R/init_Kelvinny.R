#' documenting kelvinny shortcut
#'
#' @return quickly redocumenting kelvinny
#' @import devtools
#' @import roxygen2
#' @export
init_Kelvinny <- function() {
document()
setwd('..')
install('kelvinny')
setwd('kelvinny')
}
