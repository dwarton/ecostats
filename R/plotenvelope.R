#' Diagnostic Plots for a Fitted Object with Simulation Envelopes
#'
#' Produces diagnostic plots of a \emph{fitted model} \code{y}, and
#' adds global envelopes constructed by simulating new residuals to
#' provide guidance as to whether departures from expected trends are
#' large relative to sampling variation if the fitted model were
#' correct. Global envelopes are constructed using the
#' \code{GET} package for simultaneous control of error rates over the
#' whole plot, and residuals can wither be simulated from the true model
#' (computationally intensive) or from the standard normal (not always
#' a smart move). Diagnostic plots presented are residual vs fits,
#' a normal quantile plot, and a scale-location plot.
#'
#' @param y is \emph{any} object that responds to the \code{residuals} and \code{simulate}
#' functions. The function was designed with models in mind whose residuals would be 
#' approximately normally distributed if the model were correct.
#' @param which which diagnostic plots to construct: 1=residual vs fits with a smoother,
#' 2=normal quantile plot, 3=scale-location plot with smoother. These are the first three
#' options in \code{\link{plot.lm}} and borrow heavily from that code.
#' @param n.sim the number of simulated sets of residuals to be generated, to which
#'  the observed residuals will be compared. The default is 199 datasets, you should bump it up to
#'  something more like 999 for publication (but note that would take about five times longer to run).
#' @param conf.level the confidence level to use in constructing the envelope.
#' @param type the type of global envelope to construct, see 
#'  \code{\link[GET]{global_envelope_test}} for details. Defaults \code{"st"} uses 
#'  studentized envelope tests to protect for unequal variance, which has performed well
#'  in \code{qqnorm} simulations.
#' @param transform a character vector pointing to a function that should be applied to both
#'  axes of the normal quantile plot. The most common use is to set \code{transform="pnorm"}
#'  for a PP-plot.
#' @param sim.method How should residuals be simulated? The default \code{"norm"} simulates from a normal
#' distribution, matching means and variances (and covariances for multivariate responses) to the observed residuals.
#' The \code{"stand.norm"} option sets means to zero and variances to one, which is appropriate when residuals
#' should be standard normal when assumptions are satisfied (as when using \code{\link[mvabund]{residuals.manyglm}}, for example).. If response is multivariate, this will use \code{\link[mvtnorm]{rmvnorm}} and
#' The \code{"refit"} option constructs new responses (via \code{\link{simulate}}),
#' refits the model (via \code{\link{update}}), then recomputes residuals, often known as a \emph{parametric bootstrap}.
#' This is computationally intensive but gives a more robust answer. This is the only suitable option if
#' residuals are not normal when assumptions are satisfied (like when using \code{\link[stats]{glm}}). 
#' @param main the plot title (if a plot is produced). A vector of three titles, one for each plot.
#'  If only one value is given that will be used for all plots.
#' @param xlab \code{x} axis label (if a plot is produced). A vector of three labels, one for each plot.
#'  If only one value is given that will be used for all plots.
#' @param ylab \code{y} axis label (if a plot is produced). A vector of three labels, one for each plot.
#'  If only one value is given that will be used for all plots.
#' @param col color of points
#' @param smooth.col color of smoother in residual vs fits and scale-location plots
#' @param line.col color of line on each plot, following the expected trend when assumptions
#' are satisfied. Defaults to "darkolivegreen". Because it's cool.
#' @param envelope.col color of the global envelope around the expected trend. The smoother
#' should always stay within this (because it will do for a proportion \code{conf.level} of 
#' datasets satisfying model assumptions). Color defaults to a transparent version of \code{line.col}.
#' @param plot.it logical. Should the result be plotted? If not, a list of analysis outputs
#' is returned, see \emph{Value}.
#' @param ... further arguments sent through to \code{plot}
#' 
#' @details
#' A challenge when interpreting diagnostic plots is understanding the extent to which
#' deviations from the expected pattern could be due to random noise (sampling variation)
#' rather than actual assumption violations. This function is intended to assess this, 
#' by simulating multiple realizations of residuals (and fitted values) in situations where 
#' assumptions are satisfied, and plotting a simulation envelope around these at level \code{conf.level}.
#' 
#' This function can take data (univariate or multivariate) and check for (multivariate)
#' normality, or it can take a fitted model and use qq plots to interrogate residuals and
#' see if they are behaving as we would expect them to if the model were true. As long as 
#' the fitted model accepts \code{\link{simulate}}, \code{\link{residuals}} and 
#' \code{\link{predict}} it should be fine.
#' 
#' The \code{which} argument determines which plots are constructed, out of a residual vs fits
#' plot, a normal quantile plot, and a scale-location plot. These are the first three options 
#' in \code{link{plot.lm}}, code was borrowed from there in constructing these plots.
#' 
#' Envelopes are generated using a method controlled by \code{sim.method}.
#' The default (\code{sim.method="norm"}) is to simulate new (multivariate) normal residuals repeatedly and hence these
#' can be used to assess whether trends in the observed data depart from what would be expected
#' for independent (multivariate) normal residuals. If residuals are expected to be standard
#' normal, a more refined check is to simulate from the standard normal using (\code{sim.method="stand.norm"}).
#' This would for example be appropriate when diagnosing a \code{\link[mvabund]{manyglm}} object, since
#' Dun-Smyth residuals are approximately standard normal when assumptions are satisfied.  
#' A more rigorous but computationally intensive approach (\code{sim.method="refit"}) is to use a 
#' \emph{parametric bootstrap} approach: simulate new responses from the fitted model, refit 
#' the model and recompute residuals and fitted values. This directly assesses whether trends 
#' in observed trends depart from what would be expected if the fitted model were correct, 
#' without any further assumptions. For complex models or large datasets this would however be super-slow.
#' 
#' Note that for multivariate normal data, \code{\link{cresiduals}} and \code{\link{cpredict}} 
#' are used to construct residuals and fitted values (respectively) from the \emph{full conditional models}
#' (that is, models constructed by regressing each response against all other responses
#' together with predictors). This is done because full conditionals are diagnostic of joint 
#' distributions, so \emph{any} violation of multivariate normality is expressed as a violation of 
#' linear model assumptions on full conditionals. By default, results for all responses are
#' overlaid on a single plot, but separate plots for each response can be constructed using
#' \code{overlay=FALSE}.
#' 
#' Simulation envelopes are global, constructed using the \code{\link[GET]{GET-package}}, meaning that
#' (for example) a 95% global envelope on a QQ plot should contain \emph{all} residuals for 95% of datasets
#' for which assumptions are actually satisfied. So if \emph{any} data points lie outside the
#' envelope we have evidence that assumptions of the fitted model are not satisfied. 
#' For residual vs fits and scale-location plots, global envelopes are constructed for
#' the \emph{smoother}, not for the data, hence we are looking to see if the smoother
#' is wholly contained within the envelope. The smoother is constructed using \code{link[mgcv]{gam}}.
#' 
#' The simulated data and subsequent analysis are also used to obtain a \emph{P}-value 
#' for the test that model assumptions are correct, for each plot. This test uses a 
#' "parametric bootstrap" test (for \code{do.fit=TRUE}), and tests if sample residuals or their smoothers
#' are unusually far from the values expected of them if model assumptions were satisfied. For details see
#' \code{\link[GET]{global_envelope_test}}.
#' 
#' @return up to three diagnostic plots with simulation envelopes are returned, and additionally a list of 
#' three objects used in plotting.
#' 
#' @author David Warton <david.warton@@unsw.edu.au>
#' 
#' @seealso \code{\link{cpredict}}, \code{\link{cresiduals}}, \code{\link{qqenvelope}}
#' @examples
#' # fit a multivariate linear model to the iris dataset:
#' data(iris)
#' Y = with(iris, cbind(Sepal.Length,Sepal.Width,Petal.Length,Petal.Width))
#' iris.mlm=lm(Y~Species,data=iris)
#' # check normality assumption:
#' plotenvelope(iris.mlm,n.sim=199)
#' 
#' @importFrom grDevices adjustcolor 
#' @importFrom graphics lines par polygon
#' @importFrom stats cor fitted formula lm model.frame model.matrix model.response predict qnorm qqnorm quantile residuals residuals.lm rnorm rstandard runif update var
#' @export
plotenvelope = function (y, which = 1:3, main = c("Residuals vs Fitted Values", "Normal Quantile Plot", "Scale-Location Plot"), xlab = c("Fitted values", "Theoretical Quantiles", "Fitted Values"),
                       ylab = c("Residuals", "Residuals", expression(sqrt("|Residuals|"))), 
                       n.sim=199, conf.level=0.95, type="st", transform = NULL, sim.method="norm", col=par("col"), smooth.col="black", line.col="darkolivegreen3", envelope.col = adjustcolor(line.col, 0.1),
                       plot.it = TRUE, ...) 
{
  ### which plots to construct? And clean up arguments
  show = rep(FALSE, 3)
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
  
  resFunction = 
  if(inherits(object,"mlm")) 
    cresiduals
  else
  {
    if(inherits(object,"lm"))
      rstandard
    else 
      residuals
  } 
  predFunction = 
  if(inherits(object,"mlm")) 
    cpredict
  else
  {
    predict
  } 
  y = resFunction(object)
  x = predFunction(object)

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
    sigma  = switch(sim.method,"stand.norm" = 1, var(y))
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
    for(i.sim in 1:n.sim)
    {
      modelF[1]      = matrix(yNew[,i.sim],ncol=n.resp,dimnames=dimnames(y))
      newFit         = try(update(object,data=modelF))
      resids[,i.sim] = resFunction(newFit)
      fits[,i.sim]   = predFunction(newFit)
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
    nPred = 500
    xPred = seq(min(x),max(x),length=nPred)
    gam.yx=mgcv::gam(c(y)~c(x))
    resPred=predict(gam.yx,newdata=data.frame(x=xPred))

    # get smoothers for simulated data
    yiPred=matrix(NA,nPred,n.sim)
    for(i.sim in 1:n.sim)
    {
      xiSim=fits[,i.sim]
      gam.yxi=mgcv::gam(resids[,i.sim]~xiSim)
      yiPred[,i.sim]=predict(gam.yxi,newdata=data.frame(xiSim=xPred))
    }
    
    #use the Global Envelope Test package to get global envelope
    ySort = sort(y, index.return=TRUE)
    datCurves = GET::create_curve_set(list(obs=as.vector(resPred), sim_m=yiPred))
    cr = GET::global_envelope_test(datCurves, type=type, alpha=1-conf.level)

    if(plot.it==TRUE)      #do a qq plot and add sim envelope
    {
      # make qqplot
      plot(x, y, col=col, main=main[1], 
           xlab=xlab[1], ylab=ylab[1], ...)
      lines(range(x),c(0,0),col=line.col)
      # add smooth
      lines(xPred,resPred,col=smooth.col)
      
      # add envelope
      polygon(xPred[c(1:nPred,nPred:1)], c(cr$lo,cr$hi[nPred:1]), 
              col=envelope.col, border=NA)
    }
    
    # return a list with smoother predictions and limits
    out[[1]]=list(x = xPred, y = resPred, lo=cr$lo, hi=cr$hi, p.value=attr(cr,"p"))
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
    # get observed smoother
    nPred   = 500
    xPred   = seq(min(x),max(x),length=nPred)
    yAbs    = sqrt(abs(y))
    gam.yx  = mgcv::gam(c(yAbs)~c(x))
    resPred = predict(gam.yx,newdata=data.frame(x=xPred))
    
    # get smoothers for simulated data
    yiPred=matrix(NA,nPred,n.sim)
    residsAbs = sqrt(abs(resids))
    for(i.sim in 1:n.sim)
    {
      xiSim          = fits[,i.sim]
      gam.yxi        = mgcv::gam(residsAbs[,i.sim]~xiSim)
      yiPred[,i.sim] = predict(gam.yxi,newdata=data.frame(xiSim=xPred))
    }
    
    #use the Global Envelope Test package to get global envelope
    ySort     = sort(yAbs, index.return=TRUE)
    datCurves = GET::create_curve_set(list(obs=as.vector(resPred), sim_m=yiPred))
    cr        = GET::global_envelope_test(datCurves, type=type, alpha=1-conf.level)
    
    if(plot.it==TRUE)      #do a qq plot and add sim envelope
    {
      # make qqplot
      plot(x, yAbs, col=col, main=main[3], 
           xlab=xlab[3], ylab=ylab[3], ...)
      lines(range(x),mean(yAbs)*c(1,1),col=line.col)
      # add smooth
      lines(xPred,resPred,col=smooth.col)
      
      # add envelope
      polygon(xPred[c(1:nPred,nPred:1)], c(cr$lo,cr$hi[nPred:1]), 
              col=envelope.col, border=NA)
    }
    
    # return a list with smoother predictions and limits
    out[[3]]=list(x = xPred, y = resPred, lo=cr$lo, hi=cr$hi, p.value=attr(cr,"p"))
  }
}
