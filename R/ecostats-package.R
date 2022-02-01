#' @keywords internal
"_PACKAGE"

#' Code and Data Accompanying the Eco-Stats Text
#'
#' Functions and data supporting the Ecostats text (Warton 2022, Springer), and solutions to 
#' exercises. Functions include tools for using simulation envelopes in diagnostic plots, and a
#' function for diagnostic plots of multivariate linear models. Datasets mentioned in the text
#' are included here (where not available elsewhere) and vignette solutions to textbook exercises.
#'
#' @docType package
#' 
#' @examples
#' # spaghetti plot of longitudinal data from bird exclusion field experiment:
#' data(aphids)
#' cols=c(rgb(1,0,0,alpha=0.5),rgb(0,0,1,alpha=0.5)) #transparent colours
#' with(aphids$oat, interaction.plot(Time,Plot,logcount,legend=FALSE,
#'                                col=cols[Treatment], lty=1, ylab="Counts [log(y+1) scale]",
#'                                xlab="Time (days since treatment)") )
#'                                legend("bottomleft",c("Excluded","Present"),col=cols,lty=1)
#'                                
#' # diagnostic plots for multivariate linear models:
#' data(iris)
#' iris.mlm=lm(cbind(Sepal.Length,Sepal.Width,Petal.Length,Petal.Width)~Species,data=iris)
#' # construct residual vs fits and QQ plot
#' plot(iris.mlm, which=1:2)

#' @export

