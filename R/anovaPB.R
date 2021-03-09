#' Parametric bootstrap testing to compare two models by analysis of variance (deviance)
#'
#' Compute analysis of variance (or deviance) tables for two fitted model objects,
#' with the \emph{P}-value estimated using a parametric bootstrap -- by repeatedly
#' simulating new responses from the fitted model under the null hypothesis.
#' 
#' @param objectNull is the fitted model under the null hypothesis. This can be
#' \emph{any} object that responds to \code{simulate}, \code{update} and \code{anova}.
#' @param object is the fitted model under the alternative hypothesis. This can be
#' \emph{any} object that responds to \code{update}, \code{anova} and \code{model.frame}.
#' It must be larger than \code{objectNull} and its model frame should contain
#' all the terms required to fit \code{objectNull}.
#' @param nBoot the number of simulated datasets to generate for parametric
#' bootstrap testing. Defaults to 999.
#' @param rowRef the row number where the test statistic of interest can be found
#' in the standard \code{anova} output when calling \code{anova(objectNull,object)}.
#' Defaults to \code{2}, because it is on the second row for most common models.
#' @param colRef the column number where the test statistic of interest can be found
#' in the standard \code{anova} output when calling \code{anova(objectNull,object)}.
#' Default choices have been set for common models (\code{colRef=5} for \code{lm} objects,
#' \code{colRef=6} for \code{lmer} objects and \code{colRef=4} otherwise, which is correct 
#' for \code{glm} and \code{gam} objects).
#' @param ... further arguments sent through to \code{anova}.
#' 
#' @details
#' The \code{anova} function typically relies on classical asymptotic results
#' which sometimes don't apply well, e.g. when using penalised likelihood to fit a
#' generalised additive model, or when testing whether a random effect should be
#' included in a mixed model. In such cases we can get a more accurate test by
#' using simulation to estimate the \emph{P}-value -- repeatedly simulating 
#' new sets of simulated responses, refitting the null and alternative models, and 
#' recomputing the test statistic. This process allows us to estimate the 
#' distribution of the test statistic when the null hypothesis is true, hence 
#' draw a conclusion about whether the observed test statistic is large relative
#' to what would be expected under the null hypothesis. The process is often
#' referred to as a \emph{parametric bootstrap} (Davison & Hinkley 1997), which
#' the \emph{PB} in the function name (\code{anovaPB}) is referring to.
#' 
#' This function will take any pair of fitted models, a null (\code{objectNull})
#' and an alternative (\code{object}), and use simulation to estimate the 
#' \emph{P}-value of the test of whether there is evidence against the model in
#' \code{objectNull} in favour of the model in \code{object}. This function works
#' by repeatedly performing an \code{anova} to compare \code{objectNull}  to
#' \code{object}, where the responses in the \code{model.frame} have been replaced
#' with values simulated by applying \code{simulate} to \code{objectNull}. Hence
#' for this function to work, the objects being compared need to respond 
#' to the \code{anova}, \code{simulate} and \code{model.frame} functions.
#' 
#' This function needs to be able to find the test statistic in the \code{anova}
#' output, but it is stored in different places for different types of objects.
#' It is assumed that \code{anova(objectNull,object)} returns a matrix and that the
#' test statistic is stored in the \code{(rowRef,colRef)}th element, where both
#' \code{rowRef} and \code{colRef} have been chosen to be correct for common
#' model types (\code{glm}, \code{gam}, \code{lmer}).
#' 
#' @return A matrix which has the appearance and behaviour of an object returned by 
#' \code{anova(objectNull,object)}, except that the \emph{P}-value has been computed
#' by parametric bootstrap, and the matrix now inherits class \code{anovaPB}.
#' 
#' @author David Warton <david.warton@@unsw.edu.au>
#' @references 
#' Davison A.C. & Hinkley D. V. (1997) Bootstrap methods and their application, Cambridge University Press.
#' 
#' @seealso \code{\link[stats]{anova}}
#' @examples
#' # fit a Poisson regression to random data:
#' y = rpois(50,lambda=1)
#' x = 1:50
#' rpois_glm = glm(y~x,family=poisson())
#' rpois_int = glm(y~1,family=poisson())
#' anovaPB(rpois_int,rpois_glm,n.sim=99)
#' 
#' @importFrom stats anova coef model.frame printCoefmat simulate update

#' @export

anovaPB=function(objectNull, object, nBoot=999, colRef = switch(class(object)[1],"lm"=5,"lmerMod"=6,4), rowRef=2,...)
{
  # check the second model is the larger one... otherwise this will be a prob later
  if(length(unlist(coef(objectNull)))>length(unlist(coef(object))))
    stop(paste("Whoah, the first object fitted a larger model than the second object...","\n", 
  " it should be a smaller 'null' model fit, nested in the second object."))

  # get the observed stat
  statObs=anova(objectNull,object,...)
  
  # set up for the bootstrap: assign a stats vector and get the model frame
  stats = rep(NA,nBoot+1)
  stats[1]=statObs[rowRef,colRef]
  mf = model.frame(object)
  
  # now nBoot times, simulate new response, refit models and get anova again
  for(iBoot in 1:nBoot+1)
  {
    mf[1]=simulate(objectNull)
    objectiNull = update(objectNull, data=mf)
    objecti = update(object, data=mf)
    stats[iBoot]=anova(objectiNull,objecti,...)[rowRef,colRef]
  }  
  
  # now take the original anova table, get rid of unneeded columns, stick on P-value
  statReturn=statObs[,1:colRef] #get rid of the extra columns we don't need
  statReturn$'P-value'=NA
  statReturn$'P-value'[rowRef] = mean(stats>(stats[1]-1.e-8), na.rm=TRUE)
  attr(statReturn,"nBoot")=sum(is.na(stats)==FALSE)-1
  
  # add some bells and whistles to the anova table, useful for printing
  attr(statReturn, "heading") = attr(statObs,"heading")
  class(statReturn)=c("anovaPB",class(statObs))
  return(statReturn)
}

print.anovaPB=function(x, digits = max(getOption("digits") - 3, 3),
                       signif.stars = getOption("show.signif.stars"),
                       dig.tst = max(1, min(5, digits - 1)), ...) 
{
  
  if (!is.logical(signif.stars) || is.na(signif.stars)) {
    warning("option \"show.signif.stars\" is invalid: assuming TRUE")
    signif.stars <- TRUE
  }
  
  if (!is.null(heading <- attr(x, "heading"))) 
  { cat(heading, sep = "\n")
    cat("\n")
  }
  printCoefmat(round(x, digits=dig.tst), digits = digits, signif.stars = signif.stars, P.values = TRUE, has.Pvalue=TRUE, cs.ind = NULL, tst.ind = NCOL(x)-1, zap.ind=NCOL(x)-1, na.print = "", ...)
  cat("\n")
  cat("P-value calculated by simulating", attr(x,"nBoot"), "samples from Model 1.\n")
  
}
