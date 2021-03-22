#' Atmospheric carbon dioxide concentration from the Mauna Loa Observatory
#'
#' Monthly average measurements of carbon dioxide concentration from the Mauna Loa Observatory in Hawaii,
#' from March 1958 to February 2021. Data available courtesy of the Global Monitoring Laboratory at
#' the National Oceanic and Atmospheric Administration (NOAA) in the United States
#' (https://www.esrl.noaa.gov/gmd/ccgg/trends/data.html). 

#' @docType data
#'
#' @usage data(maunaloa)
#'
#' @format A dataframe containing:\describe{
#' \item{Date}{The data of the measurement in date format. One measurement is available for each month, the first day of the month is assumed here.}
#' \item{year}{The year of the measurement.}
#' \item{month}{The month of the measurement.}
#' \item{DateNum}{The date in numerical format, as \code{year+month/12}.}
#' \item{co2}{Carbon dioxide measurement in parts per million. Calculated as the average of all daily measurements for the month.}
#' }
#'
#' @keywords datasets
#'
#' @examples
#' data(maunaloa)
#' plot(co2~Date, type="l", data=maunaloa)
"maunaloa"