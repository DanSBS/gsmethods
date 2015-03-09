\name{smtlasso}
\alias{smtlasso}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{This is a title
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
smtlasso(x, y, grp, lambda, rho, ngam, max.iter, eps, method = "smtl", intercept = TRUE, parallel = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{grp}{
%%     ~~Describe \code{grp} here~~
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
}
  \item{rho}{
%%     ~~Describe \code{rho} here~~
}
  \item{ngam}{
%%     ~~Describe \code{ngam} here~~
}
  \item{max.iter}{
%%     ~~Describe \code{max.iter} here~~
}
  \item{eps}{
%%     ~~Describe \code{eps} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{intercept}{
%%     ~~Describe \code{intercept} here~~
}
  \item{parallel}{
%%     ~~Describe \code{parallel} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, y, grp, lambda, rho, ngam, max.iter, eps, method = "smtl", 
    intercept = TRUE, parallel = TRUE) 
{
    p <- ncol(x)
    lambda <- rep(lambda, p)
    if (intercept) 
        lambda[1] <- 0
    sx <- split(as.data.frame(x), grp)
    k <- length(sx)
    sy <- split(y, grp)
    Clist <- lapply(names(sx), function(g) {
        sxg <- as.matrix(sx[[g]])
        sxy <- sy[[g]]
        t(sxg) \%*\% sxy/nrow(sxg)
    })
    C <- do.call("cbind", Clist)
    if (method == "smtl") {
        lambda <- lambda/rowSums(abs(C))
        rho <- rep(rho, p)/rowSums(abs(C))
        if (intercept) 
            rho[1] <- 0
    }
    if (method == "sgl") 
        lambda <- lambda/sqrt(rowSums(C^2))
    Dlist <- lapply(names(sx), function(g) {
        gmat <- as.matrix(sx[[g]])
        t(gmat) \%*\% gmat/nrow(gmat)
    })
    D <- abind(Dlist, along = 3)
    if (ngam <= 1) {
        ngam <- 2
        g_zero <- TRUE
    }
    else {
        g_zero <- TRUE
    }
    gamma <- apply(C, 2, function(u) {
        maxu <- max(abs(u))
        minu <- min(abs(u))
        if (minu == 0) 
            minu <- 1e-10
        exp(seq(log(0.975 * maxu), log(1e-10), , ngam))
    })
    if (g_zero) 
        gamma <- gamma * 0
    B <- array(0, dim = c(p, k, ngam))
    if (method == "smtl") 
        Best <- array(smtl(C, D, B, max.iter, rho, lambda, gamma, 
            ngam, eps), dim = c(p, k, ngam))
    if (method == "sgl") 
        Best <- array(sgl(C, D, B, max.iter, 0, lambda, gamma, 
            ngam, eps), dim = c(p, k, ngam))
    indx <- 1:ncol(C)
    if (parallel) {
        ncores <- detectCores()
        cl <- makeCluster(ncores)
        clusterExport(cl, c("Best", "indx", "sx", "sy", "g_zero"), 
            envir = environment())
        Bopt <- parLapply(cl, indx, fun = function(u) {
            candB <- Best[, u, ]
            if (!is.null(dim(candB))) {
                if (!g_zero) 
                  candB <- (1 + rho) * candB
                predi <- as.matrix(sx[[u]]) \%*\% candB
                if (!is.null(dim(candB))) {
                  npars <- colSums(candB > 0)
                  rss <- colSums((predi - sy[[u]])^2)/(nrow(sx[[u]]) - 
                    npars)
                  score <- nrow(sx[[u]]) * log(rss) + 2 * npars
                  o <- order(score, decreasing = FALSE)[1]
                  candB[, o]
                }
                else candB
            }
            else {
                candB
            }
        })
    }
    else Bopt <- lapply(indx, FUN = function(u) {
        candB <- Best[, u, ]
        if (!is.null(dim(candB))) {
            if (!g_zero) 
                candB <- (1 + rho) * candB
            predi <- as.matrix(sx[[u]]) \%*\% candB
            if (!is.null(dim(candB))) {
                npars <- colSums(candB > 0)
                rss <- colSums((predi - sy[[u]])^2)/(nrow(sx[[u]]) - 
                  npars)
                score <- nrow(sx[[u]]) * log(rss) + 2 * npars
                o <- order(score, decreasing = FALSE)[1]
                candB[, o]
            }
            else candB
        }
        else {
            candB
        }
    })
    Bbest <- t(do.call("rbind", Bopt))
    if (parallel) 
        stopCluster(cl)
    return(Bbest)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
