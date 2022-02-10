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
#' @param n.sim the number of simulated datasets to generate for parametric
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
#' @importFrom stats anova coef logLik model.frame printCoefmat simulate update

#' @export

anovaPB=function(objectNull, object, n.sim=999, colRef = switch(class(object)[1],"lm"=5,"lmerMod"=6,4), rowRef=2,...)
{
  # check the second model is the larger one... otherwise this will be a prob later
  if(length(unlist(coef(objectNull)))>length(unlist(coef(object))))
    stop(paste("Whoah, the first object fitted a larger model than the second object...","\n", 
  " it should be a smaller 'null' model fit, nested in the second object."))

  # get the observed stat
  targs <- match.call(expand.dots = FALSE)
  anovaFn=anova
#  statObs=try(anova(targs[[2]],targs[[3]],...)) #using targs[[c(2,3)]] instead of c(objectNull,object) so user names for models are in output 
  statObs=try(anova(objectNull,object,...)) 
  if(inherits(statObs,"try-error"))
  {
    # OK so if that didn't work, let's define a new anova function via logLik
    anovaFn=function(objectNull,object,...)
    {
      llAlt=logLik(object)
      llNull=logLik(objectNull)
      table=data.frame(df=c(attr(llNull,"df"),attr(llAlt,"df")),
                       deviance=-2*c(llNull,llAlt),
                       LRT = c(NA,-2*llNull+2*llAlt))
      return(table)
    }
    # apply this function, add P-values column, change colRef to 3
    statObs=anovaFn(objectNull,object)
    statObs$P=c(NA,NA)
    names(statObs)[4]='Pr(>LRT)'
    colRef=3
    
    # add model details so the printed output looks nice
    modelnamelist <- c(deparse(targs[[2]]), deparse(targs[[3]]))
    Xnames <- list(paste(deparse(formula(objectNull),width.cutoff=500), collapse = "\n"),
                   paste(deparse(formula(object),width.cutoff=500), collapse = "\n"))
    topnote <- paste(modelnamelist, ": ", Xnames, sep = "", collapse = "\n")
    title <- "Analysis of Deviance Table\n"
    rownames(statObs)=modelnamelist
    attr(statObs, "heading") <- c(title, topnote)
  }

  # set up for the bootstrap: assign a stats vector and get the model frame
  stats = rep(NA,n.sim+1)
  stats[1]=statObs[rowRef,colRef]
#  mf = model.frame(object)
  mf <- match.call(call=object$call)
  m <- match(c("formula", "data", "subset", 
               "weights", "na.action", "etastart", 
               "mustart", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  #    mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  modelF <- try( eval(mf, parent.frame()), silent=TRUE )
  
  # if for some reason this didn't work (mgcv::gam objects cause grief) then just call model.frame on object:    
  if(inherits(modelF, "try-error"))
    modelF = model.frame(object)
  
  # if there is an offset, add it, as a separate argument when updating
  offs=try(model.offset(modelF))

  # now n.sim times, simulate new response, refit models and get anova again
  for(iBoot in 1:n.sim+1)
  {
    modelFi      = cbind(simulate(objectNull), modelF[-1])
    if(inherits(offs,"try-error"))
    {
      objectiNull  = update(objectNull, data=modelFi)
      objecti      = update(object, data=modelFi)
    }
    else
    {
      objectiNull = update(object,data=modelFi,offset=offs)
      objecti = update(object,data=modelFi,offset=offs)
    }    
    stats[iBoot] = anovaFn(objectiNull,objecti,...)[rowRef,colRef]
  }  
  
  # now take the original anova table, get rid of unneeded columns, stick on P-value
  statReturn=statObs[,1:colRef] #get rid of the extra columns we don't need
  statReturn$'P-value'=NA
  statReturn$'P-value'[rowRef] = mean(stats>(stats[1]-1.e-8), na.rm=TRUE)
  attr(statReturn,"n.sim")=sum(is.na(stats)==FALSE)-1

  # change the name of the P-value column to whatever it was in the original anova table
  hasP = grep("Pr",colnames(statObs))
  colnames(statReturn)[colRef+1]=colnames(statObs)[hasP[1]]
    
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
  { cat("\n")
    cat(heading, sep = "\n")
  }
  printCoefmat(round(x, digits=dig.tst), digits = digits, signif.stars = signif.stars, P.values = TRUE, has.Pvalue=TRUE, cs.ind = NULL, tst.ind = NCOL(x)-1, zap.ind=NCOL(x)-1, na.print = "", ...)
  cat("\n")
  cat(paste("P-value calculated by simulating ", attr(x,"n.sim"), " samples from ", rownames(x)[1],".\n",sep=""))
  
}
