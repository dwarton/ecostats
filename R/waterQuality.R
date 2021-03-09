#' Water Quality data
#'
#' Data from a study of the association between water quality and catchment area. Fish composition were used to construct an
#' index of biotic integrity (IBI) and relate it to catchment are in the Seine river basin in France.
#' 
#' @docType data
#'
#' @usage data(waterQuality)
#'
#' @format A dataframe containing:\describe{
#' \item{catchment}{The catchment area in square kilometres}
#' \item{quality}{Index of Biotic Integrity (IBI)}
#' \item{logCatchment}{log base 10 of catchment area}
#' }
#'
#' @keywords datasets
#'
#' @references Oberdorff, T. & Hughes, R. M. (1992). Modification of an index of biotic integrity based
#' on fish assemblages to characterize rivers of the Seine Basin, France". Hydrobiologia, \bold{228}, 117-130.
#'
#' @examples
#' data(waterQuality)
#' plot(quality~logCatchment, data=waterQuality)
"waterQuality"