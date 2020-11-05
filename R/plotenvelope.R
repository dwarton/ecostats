#' Diagnostic Plots for a Fitted Object with Simulation Envelopes
#'
#' Produces diagnostic plots of a fitted model \code{y}, and
#' adds global envelopes constructed by simulating new residuals, to see how
#' departures from expected trends compare to what might be expected if the 
#' fitted model were correct. Global envelopes are constructed using the
#' \code{GET} package (Myllymäki et al 2017) for simultaneous control of error rates over the
#' whole plot, and residuals can be simulated by first simulating new responses from the fitted model and recomputing residuals
#' (can be computationally intensive) or just simulating directly from the (multivariate) normal distribution
#' (faster, but not always a smart move). Options for diagnostic plots to construct are a residual vs fits,
#' a normal quantile plot, or scale-location plot, along the lines of \code{\link{plot.lm}}.
#'
#' @param y is \emph{any} object that responds to \code{residuals}, \code{predict} and 
#' (if \code{sim.method="refit"}) \code{simulate}.
#' @param which a vector specifying the diagnostic plots to construct: \enumerate{
#' \item{residual vs fits, with a smoother}
#' \item{normal quantile plot}
#' \item{scale-location plot, with smoother}
#' }
#' These are the first three options in \code{\link{plot.lm}} and borrow 
#' a little from that code. A global envelope is included with each plot around where we expect to see
#' the data (or the smoother, if requested, for plots 1 and 3) when model assumptions are satisfied.
#' If not fully captured by the global envelope, there is some evidence of assumption violations.
#' @param n.sim the number of simulated sets of residuals to be generated, to which
#'  the observed residuals will be compared. The default is 199 datasets, you should bump it up to
#'  something more like 999 for publication (but note that would take about five times longer to run).
#' @param conf.level the confidence level to use in constructing the envelope.
#' @param type the type of global envelope to construct, see 
#'  \code{\link[GET]{global_envelope_test}} for details. Default \code{"st"} uses 
#'  studentized envelope tests to protect for unequal variance, which has performed well
#'  in simulations for diagnosing normality.
#' @param transform a character vector pointing to a function that should be applied to both
#'  axes of the normal quantile plot. The most common use is to set \code{transform="pnorm"}
#'  for a PP-plot.
#' @param sim.method How should residuals be simulated? The default \code{"refit"} option constructs new responses (via \code{\link{simulate}}),
#' refits the model (via \code{\link{update}}), then recomputes residuals, often known as a \emph{parametric bootstrap}.
#' This is computationally intensive but gives a robust answer. This is the only suitable option if
#' residuals are not normal when assumptions are satisfied (like when using \code{\link[stats]{glm}} with a non-Gaussian family). In fact
#' for an object that has class \code{glm} and non-Gaussian family, \code{"refit"} will override any user choice of \code{sim.method}. 
#' Alternatively, \code{"norm"} simulates from a normal distribution, matches means and variances (and covariances for multivariate responses) to the observed residuals.
#' The \code{"stand.norm"} option sets means to zero and variances to one, which is appropriate when residuals
#' should be standard normal when assumptions are satisfied (as when using \code{residuals.manyglm} from the \code{mvabund}
#' package, for example). These options are faster but more approximate than \code{"refit"}.
#' @param main the plot title (if a plot is produced). A vector of three titles, one for each plot.
#'  If only one value is given that will be used for all plots.
#' @param xlab \code{x} axis label (if a plot is produced). A vector of three labels, one for each plot.
#'  If only one value is given that will be used for all plots.
#' @param ylab \code{y} axis label (if a plot is produced). A vector of three labels, one for each plot.
#'  If only one value is given that will be used for all plots.
#' @param col color of points
#' @param line.col color of the line on diagnostic plots about which points should show no trend if model assumptions are satisfied.
#' Defaults to "olivedrab". Because it's cool. If \code{do.smooth=TRUE} this is the color of the smoother (defaulting to "slateblue4").
#' @param envelope.col color of the global envelope around the expected trend. All data points should always stay within this envelope
#' (and will for a proportion \code{conf.level} of datasets satisfying model assumptions). If a smoother has been used on
#' the residual vs fits or scale-location plot, the envelope is constructed around the smoother, that is, the smoother should always stay
#' within the envelope.
#' @param do.smooth logical defaults to \code{FALSE} (no smoother). If \code{TRUE}, a smoother is drawn on residual vs fits and
#' scale-location plots, and the global envelope is around the \emph{smoother} not the data. No smoother is added to the normal quantile plot.
#' @param plot.it logical. Should the result be plotted? If not, a list of analysis outputs is returned, see \emph{Value}.
#' @param ... further arguments sent through to \code{plot}
#' 
#' @details
#' A challenge when interpreting diagnostic plots is understanding the extent to which
#' deviations from the expected pattern could be due to random noise (sampling variation)
#' rather than actual assumption violations. This function is intended to assess this, 
#' by simulating multiple realizations of residuals (and fitted values) in situations where 
#' assumptions are satisfied, and plotting a global simulation envelope around these at level \code{conf.level}.
#' 
#' This function can take any fitted model, and construct any of three diagnostic plots, as determined by \code{which}:
#' \enumerate{
#' \item Residual vs fits plot (optionally, with a smoother)
#' \item Normal quantile plot
#' \item Scale-Location plot (optionally, with smoother)
#' }
#' and see if the trend is behaving as expected if the model were true. As long as 
#' the fitted model responds to \code{\link{residuals}} and \code{\link{predict}} 
#' (and when \code{sim.method="refit"}, \code{\link{simulate}}) then a simulation envelope
#' will be constructed for each plot.
#' 
#' Simulation envelopes are global, constructed using the \code{\link[GET]{GET-package}}, meaning that
#' (for example) a 95\% global envelope on a quantile plot should contain \emph{all} residuals for 95\% of datasets
#' that satisfy model assumptions. So if \emph{any} data points lie outside the
#' quantile plot's envelope we have evidence that assumptions of the fitted model are not satisfied. 
#' The \code{\link[GET]{GET-package}} was originally constructed to place envelopes around functions, motivated by
#' the problem of diagnostic testing of spatial processes (Myllymäki et al 2017), but it can equally
#' well be applied here, by treating the set of residuals (order according to the x-axis) as point-wise evaluations of a function.
#' For residual vs fits and scale-location plots, global envelopes are constructed for
#' the \emph{smoother}, not for the data, hence we are looking to see if the smoother
#' is wholly contained within the envelope. The smoother is constructed using \code{\link[mgcv]{gam}}.
#' 
#' The method used to simulate data for the global envelopes is controlled by \code{sim.method}.
#' Unless \code{y} has class \code{glm}, the default (\code{sim.method="norm"}) is to
#' simulate new (multivariate) normal residuals repeatedly and use these to assess whether trends 
#' in the observed data depart from what would be expected
#' for independent (multivariate) normal residuals. If residuals are expected to be standard
#' normal, a more refined check is to simulate from the standard normal using (\code{sim.method="stand.norm"}).
#' This might for example be appropriate when diagnosing a \code{manyglm} object (from the \code{mvabund} package),
#' since Dunn-Smyth residuals are approximately standard normal when assumptions are satisfied.  
#' A more rigorous but computationally intensive approach (\code{sim.method="refit"}) is to use a 
#' \emph{parametric bootstrap} approach: simulate new responses from the fitted model, refit 
#' the model and recompute residuals and fitted values. This directly assesses whether trends 
#' in observed trends depart from what would be expected if the fitted model were correct, 
#' without any further assumptions. For complex models or large datasets this would however be super-slow.
#' If \code{y} is a \code{glm} with non-Gaussian
#' family then residuals will not be normal and \code{"refit"} is the only appropriate option, hence default
#' behavior for such a model is to use \code{sim.method="refit"} and return a warning indicating this
#' has been done. 
#' 
#' Note that for Multivariate Linear Models (\code{mlm}), \code{\link{cresiduals}} and \code{\link{cpredict}} 
#' are used to construct residuals and fitted values (respectively) from the \emph{full conditional models}
#' (that is, models constructed by regressing each response against all other responses
#' together with predictors). This is done because full conditionals are diagnostic of joint 
#' distributions, so \emph{any} violation of multivariate normality is expressed as a violation of 
#' linear model assumptions on full conditionals. By default, results for all responses are
#' overlaid on a single plot.
#' 
#' The simulated data and subsequent analysis are also used to obtain a \emph{P}-value 
#' for the test that model assumptions are correct, for each plot. This tests if sample residuals or their smoothers
#' are unusually far from the values expected of them if model assumptions were satisfied. For details see
#' \code{\link[GET]{global_envelope_test}}.
#' 
#' @return up to three diagnostic plots with simulation envelopes are returned, and additionally a list of 
#' three objects used in plotting, for plots 1-3 respectively. Each is a list with five components:\describe{
#' \item{\code{x}}{X-values used for the envelope. In plots 1 and 3, this is 500 equally spaced points covering
#' the range of fitted values. For plot 2, this is sorted normal quantiles corresponding to observed data.}
#' \item{\code{y}}{In plots 1 and 3, this is the values of the smoother corresponding to \code{x}. For plot 2,
#' this is the sorted residuals.}
#' \item{\code{lo}}{The lower bound on the global envelope, for each value of \code{x}.}
#' \item{\code{hi}}{The upper bound on the global envelope, for each value of \code{x}.}
#' \item{\code{p.value}}{A \emph{P}-value for the test that the observed smoother or data is not unusually far
#' from what was expected, computed in \code{\link[GET]{global_envelope_test}}.} 
#' }
#' 
#' @author David Warton <david.warton@@unsw.edu.au>
#' @references 
#' Myllymäki, M., Mrkvička, T., Grabarnik, P., Seijo, H. and Hahn, U. (2017), Global envelope tests for spatial processes. J. R. Stat. Soc. B, 79: 381-404. doi:10.1111/rssb.12172
#' 
#' @seealso \code{\link{cpredict}}, \code{\link{cresiduals}}, \code{\link{qqenvelope}}
#' @examples
#' # fit a Poisson regression to random data:
#' y = rpois(50,lambda=1)
#' x = 1:50
#' rpois_glm = glm(y~x,family=poisson())
#' plotenvelope(rpois_glm,n.sim=99)
#' 
#' # fit a multivariate linear model to the iris dataset:
#' data(iris)
#' Y = with(iris, cbind(Sepal.Length,Sepal.Width,Petal.Length,Petal.Width))
#' iris_mlm=lm(Y~Species,data=iris)
#' # check normality assumption:
#' plotenvelope(iris_mlm,n.sim=99,which=2)
#' \donttest{plotenvelope(iris_mlm, sim.method="refit",which=1:3)
#' ## This line takes a few seconds to run.
#' ## Note violation on the scale/location plot.}
#' 
#' @importFrom grDevices adjustcolor 
#' @importFrom graphics plot points lines par polygon
#' @importFrom stats cor fitted formula lm model.frame model.matrix model.response predict qnorm qqnorm quantile residuals residuals.lm rnorm rstandard runif sd update var

#' @export
plotenvelope = function (y, which = 1:2, sim.method="refit", 
                       n.sim=199, conf.level=0.95, type="st", transform = NULL, 
                       main = c("Residuals vs Fitted Values", "Normal Quantile Plot", "Scale-Location Plot"), xlab = c("Fitted values", "Theoretical Quantiles", "Fitted Values"),
                       ylab = c("Residuals", "Residuals", expression(sqrt("|Residuals|"))), col=par("col"), 
                       line.col=if(do.smooth) "slateblue4" else "olivedrab", envelope.col = adjustcolor(line.col, 0.1), do.smooth=FALSE,
                       plot.it = TRUE, ...) 
{
  ### which plots to construct? And clean up arguments
  show = rep(FALSE, 6)
  show[which] = TRUE
  if(length(main)==1) main = rep(main,3)
  if(length(xlab)==1) xlab = rep(xlab,3)
  if(length(ylab)==1) ylab = rep(ylab,3)

  ### first work out what object we have and get observed residuals and fits ###
  
  # is y data or model fit? Make y data, object model fit
  if(is.numeric(y))
    stop("y must be an object obtained by fitting a model to data, it cannot be the dataset itself")
  else
  {
    object = y
    y      = model.response(model.frame(object))
  }
  
  # define residual function to be cresiduals for mlm, rstandard for lm, otherwise residuals
  resFunction = 
  if(inherits(object,"mlm")) 
    cresiduals
  else
  {
    if(all(class(object)=="lm")) 
      rstandard
    else 
      residuals
  } 

  # define predict function to be cpredict for mlm, otherwise predict
  predFunction = 
  if(inherits(object,"mlm")) 
    cpredict
  else
  {
    predict
  } 
  y = resFunction(object)
  x = predFunction(object)

  # if it's a non-Gaussian GLM, change sim.method to "refit" with warning if required
  if(inherits(object,"glm")) 
  {
    if(sim.method!="refit" & object$family$family!="gaussian")
    {
      sim.method="refit"
      warning("A GLM for non-Gaussian data has been detected. Residuals will not be Gaussian when the model is correct, so sim.method has been changed to refit")
    }
  }
  
  # check if it is multivariate data, if so set a flag for later and vectorise res
  if(length(dim(y))>1)
  {
    n.resp = dim(y)[2]
    n.rows = dim(y)[1]
    is.mva = TRUE
    mu     = switch(sim.method, "stand.norm" = rep(0, n.resp), colMeans(y) )
    Sigma  = switch(sim.method, "stand.norm" = cor(y), var(y) )
    if(col==par("col"))
      col  = rep(1:n.resp,each=n.rows)
  }
  else
  {
    is.mva = FALSE
    n.resp = 1
    mu     = switch(sim.method,"stand.norm" = 0, mean(y))
    sigma  = switch(sim.method,"stand.norm" = 1, sd(y))
  }
  
  ### now simulate new residuals and get fitted values too ###
  n.obs=length(y)
  
  # simulate to get limits around these
  if(sim.method=="refit")
  {
    modelF = model.frame(object)
    yNew   = simulate(object,n.sim)
    if(is.mva) # for multivariate data, vectorise res for each sim dataset
      yNew = apply(yNew,3,cbind)
    resids = fits = matrix(NA,nrow(yNew),ncol(yNew))
#    logLiks=rep(NA,n.sim)
    for(i.sim in 1:n.sim)
    {
      modelF[1]      = matrix(yNew[,i.sim],ncol=n.resp,dimnames=dimnames(y))
      newFit         = try(update(object,data=modelF))
      resids[,i.sim] = resFunction(newFit)
      fits[,i.sim]   = predFunction(newFit)
#      logLiks[i.sim] = logLik(newFit)
      # if no variation in preds, add v small random noise to avoid error later
      if(var(fits[,i.sim])<1.e-8) fits[,i.sim] = fits[,i.sim] + 1.e-6*rnorm(n.obs)
    }
  }
  else  
  {
    if(is.mva) # for multivariate data, generate multivariate normal and resize to long format
    {
      rnormMat    = mvtnorm::rmvnorm(n.rows*n.sim, mean=mu, sigma=Sigma)
      resids      = array(rnormMat,c(n.rows,n.sim,n.resp)) #stack as an array, which defaults to sims second, response last
      rm(rnormMat)
      resids      = aperm(resids,c(1,3,2)) #switch so response is second, sims last
      dim(resids) = c(n.rows*n.resp,n.sim) #make a matrix with sims in columns and multivariate data vectorised ("long format")
    }
    else
      resids = matrix(rnorm(n.obs*n.sim,mean=mu,sd=sigma),ncol=n.sim)
    fits = matrix(rep(x,n.sim),ncol=n.sim)
  }
  out=vector("list",3)
  
  ### plot 1 - residuals vs fits plot ###
  
  if (show[1L])
  {
    # get observed smoother
    if(do.smooth)
    {
      nPred = 500
      xPred = seq(min(x),max(x),length=nPred)
      # get smoother for observed data
      gam.yx=mgcv::gam(c(y)~s(c(x)))
      resObs=as.vector(predict(gam.yx,newdata=data.frame(x=xPred)))

      # get smoothers for simulated data
      residFn=matrix(NA,nPred,n.sim)
      for(i.sim in 1:n.sim)
      {
        xiSim=fits[,i.sim]
        gam.yxi=mgcv::gam(resids[,i.sim]~s(xiSim))
        residFn[,i.sim]=predict(gam.yxi,newdata=data.frame(xiSim=xPred))
      }
    }
    else
    {
      # get observed smoother
      nPred=n.obs
      xSort = sort(x,index.return=TRUE)
      xPred=xSort$x
      resObs = y[xSort$ix]
      
      residFn=resids
      for(i.sim in 1:n.sim)
      {
        fitSort = sort(fits[,i.sim],index.return=TRUE)
        residFn[,i.sim]=resids[,i.sim][fitSort$ix]
      }
      
    }
    
    #use the Global Envelope Test package to get global envelope
    datCurves = GET::create_curve_set(list(obs=resObs, sim_m=residFn))
    cr = GET::global_envelope_test(datCurves, type=type, alpha=1-conf.level)
    
    if(plot.it==TRUE)      #do a res vs fits plot and add sim envelope
    {
      # make axes and scales as appropriate
      if(do.smooth) #for smoother, keep ylim to data only
      {
        plot(x,y, main=main[1], 
             xlab=xlab[1], ylab=ylab[1], type="n", ...)
      }
      else #otherwise ylim should cover envelope
      {
        plot(c(x,xPred,xPred), c(y,cr$lo,cr$hi), main=main[1], 
             xlab=xlab[1], ylab=ylab[1], type="n", ...)
      }
      # add envelope and line
      polygon(xPred[c(1:nPred,nPred:1)], c(cr$lo,cr$hi[nPred:1]), 
              col=envelope.col, border=NA)
      # add smooth, if applicable
      if(do.smooth)
        lines(xPred,resObs,col=line.col)
      else
        lines(range(x),c(0,0),col=line.col,...)
      # add data
      points(x, y, col=col, ...)
      
    }

    # return a list with smoother predictions and limits
    out[[1]]=list(x = xPred, y = resObs, lo=cr$lo, hi=cr$hi, p.value=attr(cr,"p"))
  }
    
  ### plot 2 - normal quantile plot ###
  
  if (show[2L]) 
  {
      qq.x=qqnorm(as.vector(y),plot.it=FALSE)
      qqSort=apply(resids,2,sort)
    
      # if required, transform data
      if (is.null(transform)==FALSE)
      {
        qq.x$x = do.call(transform,list(qq.x$x))
        qq.x$y = do.call(transform,list(qq.x$y))
        qqSort = do.call(transform,list(qqSort))
        y      = do.call(transform,list(y))
      }
    
      #use the Global Envelope Test package to get global envelope
      ySort = sort(y, index.return=TRUE)
      datCurves = GET::create_curve_set(list(obs=ySort$x, sim_m=qqSort))
      cr = GET::global_envelope_test(datCurves, type=type, alpha=1-conf.level)
    
      if(plot.it==TRUE)      #do a qq plot and add sim envelope
      {
        # make qqplot
        plot(qq.x$x, qq.x$y, col=col, main=main[2], 
             xlab=xlab[2], ylab=ylab[2], ...)
        #add a qq line as is done in qqline, but limited to range of data
        probs=c(0.25,0.75)
        qs=quantile(qq.x$y,probs)
        xs=qnorm(probs)
        slope=diff(qs)/diff(xs)
        int=qs[1]-slope*xs[1]
        lines(range(qq.x$x),int+slope*range(qq.x$x),col=line.col)
        # add envelope
        polygon(sort(qq.x$x)[c(1:n.obs,n.obs:1)], c(cr$lo,cr$hi[n.obs:1]), 
                col=envelope.col, border=NA)
      }
    
      # return a list with qq values and limits, ordered same way as input data
      qqLo = qqHi = rep(NA,n.obs)
      qqLo[ySort$ix] = cr$lo
      qqHi[ySort$ix] = cr$hi
      out[[2]]=list(x = qq.x$x, y = qq.x$y, lo=qqLo, hi=qqHi, p.value=attr(cr,"p"))
  }

  ### plot 3 - scale-location plot ###

  if (show[3L])
  {
    yAbs    = sqrt(abs(y))
    residsAbs = sqrt(abs(resids))
    # get observed smoother
    if(do.smooth)
    {
      nPred = 500
      xPred = seq(min(x),max(x),length=nPred)
      # get smoother for observed data
      gam.yx=mgcv::gam(c(yAbs)~s(c(x)))
      resObs=as.vector(predict(gam.yx,newdata=data.frame(x=xPred)))
      
      # get smoothers for simulated data
      residFn=matrix(NA,nPred,n.sim)
      for(i.sim in 1:n.sim)
      {
        xiSim          = fits[,i.sim]
        gam.yxi        = mgcv::gam(residsAbs[,i.sim]~s(xiSim))
        residFn[,i.sim] = predict(gam.yxi,newdata=data.frame(xiSim=xPred))
      }
    }
    else
    {
      # get observed smoother
      nPred=n.obs
      xSort = sort(x,index.return=TRUE)
      xPred=xSort$x
      resObs = yAbs[xSort$ix]
      
      residFn=residsAbs
      for(i.sim in 1:n.sim)
      {
        fitSort = sort(fits[,i.sim],index.return=TRUE)
        residFn[,i.sim]=residsAbs[,i.sim][fitSort$ix]
      }
    }
    
    #use the Global Envelope Test package to get global envelope
    datCurves = GET::create_curve_set(list(obs=resObs, sim_m=residFn))
    cr = GET::global_envelope_test(datCurves, type=type, alpha=1-conf.level)
    cr$lo=pmax(cr$lo,0) #ensure non-negative
    
    if(plot.it==TRUE)      #do a res vs fits plot and add sim envelope
    {
      # make axes and scales as appropriate
      if(do.smooth) #for smoother, keep ylim to data only
      {
        plot(x,yAbs, main=main[3], 
             xlab=xlab[3], ylab=ylab[3], type="n", ...)
      }
      else #otherwise ylim should cover envelope
      {
        plot(c(x,xPred,xPred), c(yAbs,cr$lo,cr$hi), main=main[3], 
             xlab=xlab[3], ylab=ylab[3], type="n", ...)
      }
      # add envelope and line
      polygon(xPred[c(1:nPred,nPred:1)], c(cr$lo,cr$hi[nPred:1]), 
              col=envelope.col, border=NA)
#      lines(range(xObs),median(yAbs)*c(1,1),col=line.col,...)
      # add smooth, if applicable
      if(do.smooth)
        lines(xPred,resObs,col=line.col,...)
      # add data
      points(x, yAbs, col=col, ...)
      
    }
    
    # return a list with smoother predictions and limits
    out[[3]]=list(x = xPred, y = resObs, lo=cr$lo, hi=cr$hi, p.value=attr(cr,"p"))
  }
  
  
  # an envelope around points on res vs fits plot
  if (show[4L])
  {
    # get observed smoother
    xSort = sort(x,index.return=TRUE)
    ySort = y[xSort$ix]

    residSort=resids
    for(i.sim in 1:n.sim)
    {
      fitSort = sort(fits[,i.sim],index.return=TRUE)
      residSort[,i.sim]=resids[,i.sim][fitSort$ix]
    }
      
    #use the Global Envelope Test package to get global envelope
    datCurves = GET::create_curve_set(list(obs=ySort, sim_m=residSort))
    cr = GET::global_envelope_test(datCurves, type=type, alpha=1-conf.level)
      
    if(plot.it==TRUE)      #do a qq plot and add sim envelope
    {
      # make axes and scales as appropriate
      plot(rep(xSort$x,3), c(ySort,cr$lo,cr$hi), main=main[1], 
           xlab=xlab[1], ylab=ylab[1], type="n", ...)
      # add data
      points(xSort$x, ySort, col=col, ...)
      polygon(xSort$x[c(1:n.obs,n.obs:1)], c(cr$lo,cr$hi[n.obs:1]), 
              col=envelope.col, border=NA)
    }
      
    # return a list with qq values and limits, ordered same way as input data
    out[[4]]=list(x = xSort, y = ySort, lo=cr$lo, hi=cr$hi, p.value=attr(cr,"p"))
  }
    
  # an envelope around points on scale-location plot
  if (show[6L])
  {
    # get observed smoother
    xSort = sort(x,index.return=TRUE)
    ySort = sqrt(abs(y))[xSort$ix]

    residsAbs = sqrt(abs(resids))
    residSort = residsAbs
    for(i.sim in 1:n.sim)
    {
      fitSort = sort(fits[,i.sim],index.return=TRUE)
      residSort[,i.sim]=residsAbs[,i.sim][fitSort$ix]
    }
    
    #use the Global Envelope Test package to get global envelope
    datCurves = GET::create_curve_set(list(obs=ySort, sim_m=residSort))
    cr = GET::global_envelope_test(datCurves, type=type, alpha=1-conf.level)
    
    if(plot.it==TRUE)      #do a qq plot and add sim envelope
    {
      # make axes and scales as appropriate
      plot(rep(xSort$x,3), c(ySort,cr$lo,cr$hi), main=main[3], 
           xlab=xlab[3], ylab=ylab[3], type="n", ...)
      # add data
      points(xSort$x, ySort, col=col, ...)
      polygon(xSort$x[c(1:n.obs,n.obs:1)], c(cr$lo,cr$hi[n.obs:1]), 
              col=envelope.col, border=NA)
    }
    
    # return a list with qq values and limits, ordered same way as input data
    out[[6]]=list(x = xSort, y = ySort, lo=cr$lo, hi=cr$hi, p.value=attr(cr,"p"))
  }
#  if(sim.method=="refit")
#    out[[7]]=list(logLik=logLik(object), p.value=sum(logLiks<=logLik(object))/(n.sim+1))
  
  invisible(out)
}
