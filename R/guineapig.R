#' Guineapig data
#'
#' Data derived from a study of the effects of nicotine on the development of guineapig offspring (Johns et al 1993).
#' Ten pregnant guineapigs were injected with nicotine hydrogen tartite solution, and
#' ten control guinea pigs were injected with saline solution as a control. Learning capabilities
#' of offspring were then measured as the number of errors made when trying to remember
#' where to find food in a simple maze. Data presented  are not from the original paper- they have been generated to match the summary statistics reported in Johns et al (1993).
#'
#' @docType data
#'
#' @usage data(guineapig)
#'
#' @format A list containing two dataframes, \code{oat} and \code{wheat}, depending on the crop. Each
#' dataframe contains:\describe{
#' \item{errors}{Number of errors made}
#' \item{treatment}{N=Nicotine, C=Control}
#' }
#'
#' @keywords datasets
#'
#' @references Johns \emph{et al.} (1993) The effects of chronic prenatal exposure to nicotine on the behavior of guinea pigs (\emph{Cavia porcellus}). The Journal of
#' General Psychology \bold{120}, 49-63.
#'
#' @examples
#' data(guineapig)
#' plot(errors~treatment,data=guineapig)
"guineapig"