#' documenting kelvinny shortcut
#'
#' @return quickly redocumenting kelvinny
#' @examples
#' init_Kelvinny() # leaves folder, install, and change back to folder
#' @export
init_Kelvinny <- function() {
devtools::document()
setwd('..')
devtools::install('kelvinny')
setwd('kelvinny')
}
