#' Plot Diagnostics for a Multivariate Linear Model Object
#'
#' @description
#' A plot function for \code{mlm} objects, constructing diagnostic plots as for \code{lm}.
#' Iff multivariate linear model assumptions are satisfied, each response satisfies a linear
#' model as a function of all other response variables as well as other predictors. Hence
#' we generate residual diagnostics for each of these linear models (and optionally overlay 
#' them on a single plot).
#' 
#' Borrows heavily from \code{\link{plot.lm}}, with just a few cosmetic changes.
#'
#' @param x a \code{mlm} object, typically the result of calling \code{lm} with a matrix response.
#' @param which as for \code{\link{plot.lm}}. If a subset of the plots is required, 
#'    specify a subset of the numbers 1:6, see caption below (and the ‘Details’)
#'    for the different kinds.
#' @param overlay logical \code{TRUE} (default) means responses are overlaid on a single plot,
#'    after standardizing to ensure they are all on a comparable scale. Otherwise separate plots
#'    are generated for each response variable.
#' @param caption captions to appear above the plots, of length 6, corresponding to the 6 potential plots.
#' @param panel panel function. Default is \code{\link{gamenvelope}} to plot smoothers with
#'    confidence bands.
#' @param col color of points. When responses are overlaid on a single plot, each response
#'    has a different color.
#' @param sub.caption common title across figures. If \code{NULL}, as by default, a possibly
#'     abbreviated version of \code{deparse(x$call)} is used.
#' @param main title to each plot -- in addition to \code{caption}.
#' @param ask logical; if TRUE, the user is asked before each plot, see \code{par(ask=.)}.
#' @param ... other parameters to be passed through to plotting functions.
#' @param id.n number of points to be labelled in each plot, starting with the most extreme.
#' @param labels.id	vector of labels, from which the labels for extreme points will be chosen.
#'          NULL uses observation numbers.
#' @param cex.id magnification of point labels.
#' @param qqline logical indicating if a qqline() should be added to the normal Q-Q plot.
#' @param cook.levels	levels of Cook's distance at which to draw contours.
#' @param add.smooth logical indicating if a smoother should be added to most plots; see also panel above.
#' @param label.pos	positioning of labels, for the left half and right half of the graph respectively, for plots 1-3.
#' @param cex.caption	controls the size of caption.
#' @param cex.oma.main controls the size of the sub.caption only if that is above the figures when there is more than one.
#'
#' @details
#' This function borrows heavily from \code{\link{plot.lm}} so see that file for usage details.
#' The essential difference is that residuals and fitted values are calculated from the
#' `full conditional` model, that is, as a linear function of all other response variables as 
#' well as a function of predictors. The full conditionals of a model are diagnostic of the
#' joint distribution, and in the special case of a multivariate linear model, this means that
#' the model is correct if each response satisfies a linear model as a function of all other
#' responses and other predictors. Most of the work is done in \code{\link{predict.mlm}} and
#' \code{\link{residuals.mlm}} functions.
#' 
#' @return a graph or set of graphs
#' 
#' @author David Warton <david.warton@@unsw.edu.au> and Christopher Chung
#' 
#' @seealso \code{\link{lm}}, \code{\link{residuals.mlm}}, \code{\link{predict.mlm}}, \code{\link{plot.lm}}
#' @examples
#' data(iris)
#' # fit a mlm:
#' iris.mlm=lm(cbind(Sepal.Length,Sepal.Width,Petal.Length,Petal.Width)~Species,data=iris)
#' # construct residual vs fits and QQ plot
#' plot(iris.mlm, which=1:2)
#'
#' @importFrom grDevices adjustcolor as.graphicsAnnot dev.flush dev.hold dev.interactive devAskNewPage extendrange rgb
#' @importFrom graphics abline axis frame legend lines matplot matpoints mtext par plot points polygon strheight strwidth text title
#' @importFrom stats cooks.distance deviance df.residual fitted influence lm model.frame model.matrix model.response predict qnorm qqnorm residuals residuals.lm rnorm rstandard rstudent runif terms update var weights
#' @export
plot.mlm=function(x, which = c(1L:3L, 5L), overlay = TRUE, caption = list("Residuals vs Fitted", 
                                                   "Normal Q-Q", "Scale-Location", "Cook's distance", "Residuals vs Leverage", 
                                                   expression("Cook's dist vs Leverage  " * h[ii]/(1 - h[ii]))), 
                  panel = if (add.smooth == T) gamenvelope else points,
                  col=if(overlay) 1:dim(x$model[[1]])[2] else 1,
                  sub.caption=NULL, main="",
                  ask = prod(par("mfcol")) < length(which) && dev.interactive(), 
                  ..., id.n = 3, labels.id = names(residuals(x)), cex.id = 0.75, 
                  qqline = TRUE, cook.levels = c(0.5, 1), add.smooth=TRUE,
                  label.pos = c(4, 2), cex.caption = 1, cex.oma.main = 1.25
                  )
{
  plottedObject = if(overlay == TRUE) 
    plot.mlm.default(x, which = which, caption=caption, 
                     panel = panel, sub.caption = sub.caption, 
                     main = main, ask = ask, 
                     ..., id.n = id.n, labels.id = labels.id, cex.id = cex.id, 
                     qqline = qqline, cook.levels = cook.levels, add.smooth=add.smooth,
                     label.pos = label.pos, cex.caption = cex.caption, cex.oma.main = cex.oma.main)
  else
    plotmlm_separate(x, which = which, caption=caption, 
               panel = panel, sub.caption = sub.caption, 
               main = main, ask = ask, 
               ..., id.n = id.n, labels.id = labels.id, cex.id = cex.id, 
               qqline = qqline, cook.levels = cook.levels, add.smooth=add.smooth,
               label.pos = label.pos, cex.caption = cex.caption, cex.oma.main = cex.oma.main)
} 


plot.mlm.default = function(x, which = c(1L:3L, 5L), caption = list("Residuals vs Fitted", 
                                                                            "Normal Q-Q", "Scale-Location", "Cook's distance", "Residuals vs Leverage", 
                                                                            expression("Cook's dist vs Leverage  " * h[ii]/(1 - h[ii]))), 
                    panel = if (add.smooth == T) gamenvelope else points,
                    col=1:dim(x$model[[1]])[2],
                    sub.caption=NULL, main="",
                    ask = prod(par("mfcol")) < length(which) && dev.interactive(), 
                    ..., id.n = 3, labels.id = names(residuals(x)), cex.id = 0.75, 
                    qqline = TRUE, cook.levels = c(0.5, 1), add.smooth=TRUE,
                    label.pos = c(4, 2), cex.caption = 1, cex.oma.main = 1.25
  )
  {
  
  dropInf <- function(x, h) {
    if (any(isInf <- h >= 1)) {
      warning(gettextf("not plotting observations with leverage one:\n  %s", 
                       paste(which(isInf), collapse = ", ")), call. = FALSE, 
              domain = NA)
      x[isInf] <- NaN
    }
    x
  }
  if (!inherits(x, "lm")) 
    stop("use only with \"lm\" objects")
  if (!is.numeric(which) || any(which < 1) || any(which > 6)) 
    stop("'which' must be in 1:6")
  isGlm <- inherits(x, "glm")
  show <- rep(FALSE, 6)
  show[which] <- TRUE
  r <- residuals(x,conditional=TRUE,standardize=TRUE) #changed to ensure conditional resids
  yh <- predict(x,conditional=TRUE,standardize=TRUE)  # same # reverse it in our functions
  w <- weights(x)
  if (!is.null(w)) {
    wind <- w != 0
    r <- r[wind]
    yh <- yh[wind]
    w <- w[wind]
    labels.id <- labels.id[wind]
  }
  n <- length(r)
  if (any(show[2L:6L])) {
    s <- if (inherits(x, "rlm")) 
      x$s
    else if (isGlm)  
      sqrt(summary(x)$dispersion)
    else sqrt(deviance(x)/df.residual(x))
    hii <- (infl <- influence(x, do.coef = FALSE))$hat
    if (any(show[4L:6L])) {
      cook <- if (isGlm) 
        cooks.distance(x, infl = infl)
      else cooks.distance(x, infl = infl, sd = s, res = r, 
                          hat = hii)
    }
  }
  if (any(show[2L:3L])) {
    ylab23 <- if (isGlm) 
      "Std. deviance resid."
    else "Standardized residuals"
    r.w <- if (is.null(w)) 
      r
    else sqrt(w) * r
#DW change, 12/7/19: I think r was already standardised so dont divide by s again
    #    rs <- dropInf(r.w/(s * sqrt(1 - hii)), hii)
    rs <- dropInf(r.w/sqrt(1 - hii), hii)
# end DW change
  }
  if (any(show[5L:6L])) {
    r.hat <- range(hii, na.rm = TRUE)
    isConst.hat <- all(r.hat == 0) || diff(r.hat) < 1e-10 * 
      mean(hii, na.rm = TRUE)
  }
  if (any(show[c(1L, 3L)])) 
    l.fit <- if (isGlm) 
      "Predicted values"
  else "Fitted values (standardized)" #sandardise added to label since standardised across predictors
  if (is.null(id.n)) 
    id.n <- 0
  else {
    id.n <- as.integer(id.n)
    if (id.n < 0L || id.n > n) 
      stop(gettextf("'id.n' must be in {1,..,%d}", n), 
           domain = NA)
  }
  if (id.n > 0L) {
    if (is.null(labels.id)) 
      labels.id <- paste(1L:n)
    iid <- 1L:id.n
    show.r <- sort.list(abs(r), decreasing = TRUE)[iid]
    if (any(show[2L:3L])) 
      show.rs <- sort.list(abs(rs), decreasing = TRUE)[iid]
    text.id <- function(x, y, ind, adj.x = TRUE) {
      labpos <- if (adj.x) 
        label.pos[1 + as.numeric(x > mean(range(x)))]
      else 3
      text(x, y, labels.id[ind], cex = cex.id, xpd = TRUE, 
           pos = labpos, offset = 0.25)
    }
  }
  getCaption <- function(k) if (length(caption) < k) 
    NA_character_
  else as.graphicsAnnot(caption[[k]])
  if (is.null(sub.caption)) {
    cal <- x$call
    if (!is.na(m.f <- match("formula", names(cal)))) {
      cal <- cal[c(1, m.f)]
      names(cal)[2L] <- ""
    }
    cc <- deparse(cal, 80)
    nc <- nchar(cc[1L], "c")
    abbr <- length(cc) > 1 || nc > 75
    sub.caption <- if (abbr) 
      paste(substr(cc[1L], 1L, min(75L, nc)), "...")
    else cc[1L]
  }
  one.fig <- prod(par("mfcol")) == 1
  if (ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }
  if (show[1L]) {
    ylim <- range(r, na.rm = TRUE)
    if (id.n > 0) 
      ylim <- extendrange(r = ylim, f = 0.08)
    dev.hold()
    plot(yh, r, xlab = l.fit, ylab = "Standardized residuals", main = main, #ylab= "standardized residuals"
         ylim = ylim, type = "n", ...)
    panel(yh, r, ...)
    if (one.fig) 
      title(sub = sub.caption, ...)
    mtext(getCaption(1), 3, 0.25, cex = cex.caption)
    if (id.n > 0) {
      y.id <- r[show.r]
      y.id[y.id < 0] <- y.id[y.id < 0] - strheight(" ")/3
      text.id(yh[show.r], y.id, show.r)
    }
    abline(h = 0, lty = 3, col = "gray")
    dev.flush()
  }
  if (show[2L]) {
    ylim <- range(rs, na.rm = TRUE)
    ylim[2L] <- ylim[2L] + diff(ylim) * 0.075
    dev.hold()
    qq <- qqenvelope(rs, main = main, ylab = ylab23, ylim = ylim, #changed from qqnorm ->qqnormEnvelope. But this seems to miss only in the overlay=T version?
                         ...)
    if (qqline) 
      qqline(rs, lty = 3, col = "gray50")
    if (one.fig) 
      title(sub = sub.caption, ...)
    mtext(getCaption(2), 3, 0.25, cex = cex.caption)
    if (id.n > 0) 
      text.id(qq$x[show.rs], qq$y[show.rs], show.rs)
    dev.flush()
  }
  if (show[3L]) {
    sqrtabsr <- sqrt(abs(rs))
    ylim <- c(0, max(sqrtabsr, na.rm = TRUE))
    yl <- as.expression(substitute(sqrt(abs(YL)), list(YL = as.name(ylab23))))
    yhn0 <- if (is.null(w)) 
      yh
    else yh[w != 0]
    dev.hold()
    plot(yhn0, sqrtabsr, xlab = l.fit, ylab = yl, main = main, 
         ylim = ylim, type = "n", ...)
    panel(yhn0, sqrtabsr, ...)
    if (one.fig) 
      title(sub = sub.caption, ...)
    mtext(getCaption(3), 3, 0.25, cex = cex.caption)
    if (id.n > 0) 
      text.id(yhn0[show.rs], sqrtabsr[show.rs], show.rs)
    dev.flush()
  }
  if (show[4L]) {
    if (id.n > 0) {
      show.r <- order(-cook)[iid]
      ymx <- cook[show.r[1L]] * 1.075
    }
    else ymx <- max(cook, na.rm = TRUE)
    dev.hold()
    plot(cook, type = "h", ylim = c(0, ymx), main = main, 
         xlab = "Obs. number", ylab = "Cook's distance", ...)
    if (one.fig) 
      title(sub = sub.caption, ...)
    mtext(getCaption(4), 3, 0.25, cex = cex.caption)
    if (id.n > 0) 
      text.id(show.r, cook[show.r], show.r, adj.x = FALSE)
    dev.flush()
  }
  if (show[5L]) { ### maybe try deleting this
    ylab5 <- if (isGlm) 
      "Std. Pearson resid."
    else "Standardized residuals"
    r.w <- r # changed to use residuals as previously computed 
    if (!is.null(w)) 
      r.w <- r.w[wind]
    rsp <- dropInf(r.w/(s * sqrt(1 - hii)), hii)
    ylim <- range(rsp, na.rm = TRUE)
    if (id.n > 0) {
      ylim <- extendrange(r = ylim, f = 0.08)
      show.rsp <- order(-cook)[iid]
    }
    do.plot <- TRUE
    if (isConst.hat) {
      if (missing(caption)) 
        caption[[5L]] <- "Constant Leverage:\n Residuals vs Factor Levels"
      aterms <- attributes(terms(x))
      dcl <- aterms$dataClasses[-aterms$response]
      facvars <- names(dcl)[dcl %in% c("factor", "ordered")]
      mf <- model.frame(x)[facvars]
      if (ncol(mf) > 0) {
        dm <- data.matrix(mf)
        nf <- length(nlev <- unlist(unname(lapply(x$xlevels, 
                                                  length))))
        ff <- if (nf == 1) 
          1
        else rev(cumprod(c(1, nlev[nf:2])))
        facval <- (dm - 1) %*% ff
        xx <- facval
        dev.hold()
        matplot(facval, rsp, xlim = c(-1/2, sum((nlev - 
                                                   1) * ff) + 1/2), ylim = ylim, xaxt = "n", main = main, 
                xlab = "Factor Level Combinations", ylab = ylab5, 
                type = "n", ...) #changed to matplot so it can handle matrix response and vector predictor
        axis(1, at = ff[1L] * (1L:nlev[1L] - 1/2) - 1/2, 
             labels = x$xlevels[[1L]])
        mtext(paste(facvars[1L], ":"), side = 1, line = 0.25, 
              adj = -0.05)
        abline(v = ff[1L] * (0:nlev[1L]) - 1/2, col = "gray", 
               lty = "F4")
        panel(facval, rsp, ...)
        abline(h = 0, lty = 3, col = "gray")
        dev.flush()
      }
      else {
        message(gettextf("hat values (leverages) are all = %s\n and there are no factor predictors; no plot no. 5", 
                         format(mean(r.hat))), domain = NA)
        frame()
        do.plot <- FALSE
      }
    }
    else {
      xx <- hii
      xx[xx >= 1] <- NA
      dev.hold()
      matplot(xx, rsp, xlim = c(0, max(xx, na.rm = TRUE)), 
              ylim = ylim, main = main, xlab = "Leverage", 
              ylab = ylab5, type = "n", ...) #changed to matplot so it cn handle matrix response, vector predictor
      matpoints(xx, rsp, pch=1,...) #as above... but better to change panel.smooth definition
      abline(h = 0, v = 0, lty = 3, col = "gray")
      if (one.fig) 
        title(sub = sub.caption, ...)
      if (length(cook.levels)) {
        p <- x$rank
        usr <- par("usr")
        hh <- seq.int(min(r.hat[1L], r.hat[2L]/100), 
                      usr[2L], length.out = 101)
        for (crit in cook.levels) {
          cl.h <- sqrt(crit * p * (1 - hh)/hh)
          lines(hh, cl.h, lty = 2, col = 2)
          lines(hh, -cl.h, lty = 2, col = 2)
        }
        legend("bottomleft", legend = "Cook's distance", 
               lty = 2, col = 2, bty = "n")
        xmax <- min(0.99, usr[2L])
        ymult <- sqrt(p * (1 - xmax)/xmax)
        aty <- sqrt(cook.levels) * ymult
        axis(4, at = c(-rev(aty), aty), labels = paste(c(rev(cook.levels), 
                                                         cook.levels)), mgp = c(0.25, 0.25, 0), las = 2, 
             tck = 0, cex.axis = cex.id, col.axis = 2)
      }
      dev.flush()
    }
    if (do.plot) {
      mtext(getCaption(5), 3, 0.25, cex = cex.caption)
      #      if (id.n > 0) { #commenting out text labels because of problems aligning with obs in matrix response (resolvable but not a priority for us)
      #        y.id <- rsp[show.rsp]
      #        y.id[y.id < 0] <- y.id[y.id < 0] - strheight(" ")/3
      #        text.id(xx[show.rsp], y.id, show.rsp)
      #      }
    }
  }
  if (show[6L]) {
    g <- dropInf(hii/(1 - hii), hii)
    ymx <- max(cook, na.rm = TRUE) * 1.025
    dev.hold()
    plot(g, cook, xlim = c(0, max(g, na.rm = TRUE)), ylim = c(0, 
                                                              ymx), main = main, ylab = "Cook's distance", xlab = expression("Leverage  " * 
                                                                                                                               h[ii]), xaxt = "n", type = "n", ...)
    panel(g, cook, ...)
    athat <- pretty(hii)
    axis(1, at = athat/(1 - athat), labels = paste(athat))
    if (one.fig) 
      title(sub = sub.caption, ...)
    p <- x$rank
    bval <- pretty(sqrt(p * cook/g), 5)
    usr <- par("usr")
    xmax <- usr[2L]
    ymax <- usr[4L]
    for (i in seq_along(bval)) {
      bi2 <- bval[i]^2
      if (p * ymax > bi2 * xmax) {
        xi <- xmax + strwidth(" ")/3
        yi <- bi2 * xi/p
        abline(0, bi2, lty = 2)
        text(xi, yi, paste(bval[i]), adj = 0, xpd = TRUE)
      }
      else {
        yi <- ymax - 1.5 * strheight(" ")
        xi <- p * yi/bi2
        lines(c(0, xi), c(0, yi), lty = 2)
        text(xi, ymax - 0.8 * strheight(" "), paste(bval[i]), 
             adj = 0.5, xpd = TRUE)
      }
    }
    mtext(getCaption(6), 3, 0.25, cex = cex.caption)
    if (id.n > 0) { 
      show.r <- order(-cook)[iid]
      text.id(g[show.r], cook[show.r], show.r)
    }
    dev.flush()
  }
  if (!one.fig && par("oma")[3L] >= 1) 
    mtext(sub.caption, outer = TRUE, cex = cex.oma.main)
  invisible()
}
 
# function by DW and Chris Chung to regress each response against predictors+X then plot separately        
plotmlm_separate = function(x, which = c(1L:3L, 5L), caption = list("Residuals vs Fitted", 
                                                                    "Normal Q-Q", "Scale-Location", "Cook's distance", "Residuals vs Leverage", 
                                                                    expression("Cook's dist vs Leverage  " * h[ii]/(1 - h[ii]))), 
                            panel = if (add.smooth == T) gamenvelope else points,
                            col=1,
                            sub.caption=NULL, main="",
                            ask = prod(par("mfcol")) < length(which) && dev.interactive(), 
                            ..., id.n = 3, labels.id = rownames(residuals(x)), cex.id = 0.75, 
                            qqline = TRUE, cook.levels = c(0.5, 1), add.smooth=TRUE,
                            label.pos = c(4, 2), cex.caption = 1, cex.oma.main = 1.25)
{
  mf = model.frame(x) 
  Ys = model.response(mf)
  Ys = as.matrix(Ys) # think this fixes dim = NULL problemsfit
  X = model.matrix(x) # gets the design matrix
  nYs = dim(Ys)[2]
  for (iVar in 1:nYs){
    y = Ys[ ,iVar]
    XandY = data.frame(X, Ys[, -iVar])
    fit = lm(y~., data = XandY)
    captionI=paste0(caption,": ",colnames(Ys)[iVar])
    plot.mlm.default(fit, which = which, caption=captionI, 
                     panel = panel, sub.caption = sub.caption, 
                     main = "", ask = ask, 
                     ..., id.n = id.n, labels.id = labels.id, cex.id = cex.id, 
                     qqline = qqline, cook.levels = cook.levels, add.smooth=add.smooth,
                     label.pos = label.pos, cex.caption = cex.caption, cex.oma.main = cex.oma.main)
  }
}