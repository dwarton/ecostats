#' Effect of pollution on marine microinvertebrates in estuaries in different zones
#'
#' Data from an observational study of whether there is a different in microinvertebrate
#' communities between estuaries that have been heavily modified by human activity and those
#' that have not, across seven estuaries along the coast of New South Wales, Australia (Clark et al. 2015).
#' Sampling was undertaken in both Inner and Outer zones of the estuary. 
#' 
#' @docType data
#'
#' @usage data(estuaryZone)
#'
#' @format A dataframe containing (amongst other things):\describe{
#' \item{Mod}{A factor describing whether the sample was taken from a 'Modified' or 'Pristine' estuary.}
#' \item{Zone}{Whether the sample was taken from Inner (upstream) or Outer (downstream) zone of the estuary.}
#' \item{Estuary}{A factor with seven levels identifying which estuary the sample was taken from.}
#' \item{Total}{Total abundance of all invertebrates in the sample} 
#' \item{Richness}{Richness of taxa in the sample -- the number of responses (of those in columns 8-94) taking a non-zero value} 
#' }
#' Other variables in the dataset give invertebrate counts separately for different taxa.
#'
#' @keywords datasets
#'
#' @references Clark, G. F., Kelaher, B. P., Dafforn, K. A., Coleman, M. A., Knott, N. A., 
#' Marzinelli, E. M., & Johnston, E. L. (2015). What does impacted look like? high diversity and
#' abundance of epibiota in modified estuaries. Environmental Pollution 196, 12-20.
#'
#' @examples
#' data(estuaryZone)
#' cols=c("blue","red","lightblue","pink")
#' plot(Total~interaction(Estuary,Zone),data=estuaryZone,col=cols[c(1,2,2,1,2,1,2,3,4,4,3,4,3,4)])
"estuaryZone"