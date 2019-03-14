#' Generate default ggplot colour palette depending on the number of colours required.
#' 
#' @param n Number.
#' @return Generate default ggplot colour palette depending on the number of colours required.
#' @examples
#' col <- gg_color_hue(5)
#' col
#' [1] "#F8766D" "#A3A500" "#00BF7D" "#00B0F6" "#E76BF3"
#' @export
gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }