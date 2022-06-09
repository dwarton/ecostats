#' Germination rates of \emph{Abutilon angulatum} at different temperatures
#'
#' Germination rates of \emph{Abutilon angulatum} from 29 different studies, undertaken
#' at different ambient temperatures. We would like to know how germination rate
#' varies with temperature, in particular, the range of temperatures at which \emph{Abutilon angulatum} 
#' will germinate. These data come from a larger study across species and environments
#' to look for a latitudinal signal in tolerance to changing temperature (Sentinella et al 2020).

#' @docType data
#'
#' @usage data(seedsTemp)
#'
#' @format A dataframe containing:\describe{
#' \item{NumSown}{The number of seeds sown.}
#' \item{NumGerm}{The number of seeds that germinated.}
#' \item{Test.Temp}{The ambient temperature (in degrees Celsius) of the location at which seeds were sown.}
#' }
#'
#' @keywords datasets
#'
#' @references Sentinella, AT, Warton, DI, Sherwin, WB, Offord, CA, Moles, AT. (2020) Tropical plants do not have narrower temperature tolerances, but are more at risk from warming because they are close to their upper thermal limits. Global Ecol Biogeogr. \bold{29}, 1387-1398.
#' 
#' @examples
#' data(seedsTemp)
#' seedsTemp$propGerm = seedsTemp$NumGerm / seedsTemp$NumSown
#' plot(propGerm/(1-propGerm)~Test.Temp,data=seedsTemp,log="y",
#'  ylab="Germination rate [logit scale]", xlab="Temperature (Celsius)")
"seedsTemp"