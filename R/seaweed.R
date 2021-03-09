#' Habitat Configuration data from seaweed experiment
#'
#' Data from a study of habitat configuration, specifically, does density of invertebrate epifauna on seaweed
#' vary across sites with different levels of isolation from each other.
#' 
#' @docType data
#'
#' @usage data(seaweed)
#'
#' @format A dataframe containing:\describe{
#' \item{Size}{A character vector describing size of experimental plots as "SMALL" or "LARGE"}
#' \item{Dist}{Distance of isolation -- 0, 2 or 10 metres from other algal beds}
#' \item{Time}{Sampling time - either 5 or 10 weeks from the start of the experiment.}
#' \item{Rep}{The replicate number (1 to 5).}
#' \item{Wmass}{Wet mass of the algal bed for that plot.}
#' \item{Total}{Total invertebrate density in the plot, calculated as nuber of individuals divided by \code{Wmass}.}
#' }
#' Other variables in the dataset give invertebrate counts separately for different taxa.
#'
#' @keywords datasets
#'
#' @references Roberts, D. A. & Poore, A. G. (2006). Habitat configuration affects colonisation of
#' epifauna in a marine algal bed. Biological Conservation \bold{127}, 18-26.
#'
#' @examples
#' data(seaweed)
#' boxplot(Total~Dist, data=seaweed)
"seaweed"