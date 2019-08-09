#' parse gmt file format as gmx
#' 
#' @param file file
#' @param header default TRUE
#' @param sep default "\t"
#' @param ... passes to read.csv
#' @return reads a .gmt file like a .gmx file
#' @export
#'      

parse_gmt <- function(file, header = TRUE, sep = "\t", ...) {

  n = max(count.fields(file, sep=sep), na.rm=TRUE)
  x = readLines(file)

  .splitvar = function(x, sep, n) {
    var = unlist(strsplit(x, split=sep))
    length(var) = n
    return(var)
  }

  x = do.call(cbind, lapply(x, .splitvar, sep=sep, n=n))
  x = apply(x, 1, paste, collapse=sep) 
  out = read.csv(text=x, sep=sep, header=header, ...)
  return(out)
}
