#k <- 4
#lamg <- 1e-1
#rho <- .1
#max.iter <- 10
#eps <- 1e-2
#pos <- 1
#nls <- 20
#nls_eps <- 1e-5
#thresh <- 0.975
#save(x,y,file='data/gsplass.RData')
#load('data/gsplass.RData')
require(Rcpp)
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")
sourceCpp('C:/eclipse/workspace/Tools/smtl/src/gsplasso_v4.cpp')


#gspl <-
ssgl <- function(x, y, ys, k, lamg, rho, nls, max.iter = 1000,
		eps = 1e-5, thresh = 0.975,
		parallel = FALSE, pos = 0,  
		nls_eps = NULL, ridge_df = FALSE,
		B_start = NULL, weight_ones = FALSE,
		lamsw = NULL){
	
	## Here y is a hyperspectral data cube
	## Getting dimensions
	dims <- dim(y)
	nr <- dims[1]
	nc <- dims[2]
	ns <- dims[3]
	
	## Number of variables (endmembers)
	p <- ncol(x)
	
	## converting y to a matrix
	ymat <- matrix(c(y), nr * nc, ns)
	
	
	## Check to see if element-wise sparsity has been selected
	if(nls <= 1){
		nls <- 2
		nls_zero <- TRUE
	}
	else{
		nls_zero <- FALSE
	}
	
	## This will be used to get the range for the element-wise
	## sparse solutions path
	C <- t(x) %*% t(ymat) / ns
	lams <- apply(C, 2, function(u) {
				maxu <- max(abs(u))
				minu <- min(abs(u)) 
				exp(seq(log(maxu), log(1e-6 * maxu),, nls))
			})
#	browser()
	## Get storage array for coefficient estimates
	B <- array(0, dim = c(p, nr * nc, nls))
	
	if(!is.null(B_start))
		B[,,1] <- B_start
	
	## Checks to see if/where element-wise sparse solution path
	## should be cut off (i.e. differences between regularization
	## parameters is effectively 0).
	if(!is.null(nls_eps) & !nls_zero){
		nls <- ceiling(mean(rowSums(abs(t(apply(lams, 2, diff))) > nls_eps)))
		lams <- lams[1:nls, ]
		B <- B[,,1:nls]
	}
	
	## If no element-wise sparsity just return 0
	if(nls_zero)
		lams <- lams * 0
	
	## Spatial weight matrix
	w <- 1/outer((-k:k)^2, (-k:k)^2, FUN = function(X ,Y) X + Y)
	w[!is.finite(w)] <- 0
#	browser()
	## Check to see if spatial weights should be used or just
	## spectral
	if(weight_ones){
		w <- w^0
		w[(k+1),(k+1)] <- 0
	}
	
	
	## Group sparse parameter
	if(length(lamg) < p)
		lamg <- rep(lamg, p)
	
	
#	browser()
	## Fit gsplasso
#	browser()
	system.time({Best <- array(gsplasso(x, y, B, k, pos, 
		 			lams, lamg, rho, w, nls, 
					nr, nc, ns, eps, max.iter,
					thresh, ys), 
			dim = c(p, nr * nc, nls))})
#	browser()
#	tmp <- gsplasso(x, y, B, k, pos, 
#			lams, lamg, rho, w, nls, 
#			nr, nc, ns, eps, max.iter,
#			thresh, isadapt, abs(lamsw), ys)
	
	## Grabbing index to run estimated optimal sparse solution path
	indx <- 1:ncol(C)
	
	if(parallel){
		ncores <- detectCores()
		cl <- makeCluster(ncores)
		
		if(singled == 0)
			clusterExport(cl, c('Best','indx', 'sx', 'sy', 'g_zero', 'ginv'),
					envir = environment())
		else
			clusterExport(cl, c('Best','indx', 'x', 'y', 'g_zero', 'singled', 'ginv'),
					envir = environment())
		
		Bopt <- parLapply(cl, indx, fun = function(u){
					
					candB <- Best[,u,]
					
					if(!is.null(dim(candB))){
						
						if(!g_zero)
							candB <- (1+rho[u])*candB
						
						if(singled == 0){
							xm <- sx[[u]]
							ym <- sy[[u]]
							predi <- as.matrix(sx[[u]]) %*% candB
						}
						else{
							predi <- as.matrix(x) %*% candB
							xm <- x
							ym <- y[,u]
						}
						
						if(!is.null(dim(candB))){
							c0 <- abs(candB) > 0
							if(!all(rho == 0) & ridge_df)
								npars <- apply(c0, 2, function(v) {
											xm0 <- xm[, v]
											rho0 <- rho[u]
											if(length(xm0) == 0)
												0
											else{
												if(is.null(dim(xm0))){
													xm0 <- t(rbind(xm0))
												}
												
												xtx0 <- ginv(t(xm0) %*% xm0 + rho0 * diag(ncol(xm0)))
												sum(apply(xm0, 1, function(w) t(w) %*% xtx0 %*% w))
											}
										})
							
							else
								npars <- colSums(abs(candB) > 0)
							rss <- colSums((predi - ym)^2)/(nrow(xm) - npars)
							
							## Using AIC or BIC here, haven't settled on one or the
							## other yet. The standard LASSO algorithm (or the one
							## implemented in R uses Cp, will have to think about whic
							## is most appropriate here).
							score <- nrow(xm) * log(rss) + npars * log(nrow(xm))
#							score <- nrow(xm) * log(rss) + 2 * npars
#							score <- (rss * (nrow(xm) - npars))/rss[ngam] - nrow(xm) + 2*npars						
							score[is.nan(score)] <- Inf
							o <- order(score, decreasing = FALSE)[1]
							candB[,o]
						}
						else
							candB
						
					}
					else{
						candB
					}
				})
	}
	else{
		Bopt <- lapply(indx, FUN = function(u){
					
					candB <- Best[,u,]
					
					if(!is.null(dim(candB))){
						predi <- as.matrix(x) %*% candB
						xm <- x
						ym <- ymat[u,]
						
						if(!is.null(dim(candB))){
							c0 <- abs(candB) > 0
							
							npars <- colSums(abs(candB) > 0)
							rss <- colSums((predi - ym)^2)/(nrow(xm) - npars)
							## Using AIC or BIC here, haven't settled on one or the
							## other yet. The standard LASSO algorithm (or the one
							## implemented in o uses Cp, will have to think about whic
							## is most appropriate here).
							## BIC
#							if(u == 4857)
#								browser()
#							score <- nrow(xm) * log(rss) + npars * log(nrow(xm))
							## AIC
#							score <- nrow(xm) * log(rss) + 2 * npars
							## Cp
							score <- (rss * (nrow(xm) - npars))/rss[nls] - nrow(xm) + 2 *  npars
#							score <- (rss * (nrow(xm) - npars))/rss[nls] + 2 * npars
							score[is.nan(score)] <- Inf
							o <- order(score, decreasing = FALSE)[1]
							candB[,o]
						}
						else
							candB
						
					}
					else{
						candB
					}
#					}
#					else{
#						## If the criteria is not met return a vector of zeros.
#						rep(0, ncol(beta))
#					}
					
				})
	}
	## Collect the final results
	Bbest <- t(do.call("rbind", Bopt))
	
	if(parallel)
		stopCluster(cl)
#	browser()
	return(Bbest)
	
	
}
#
#res_sm <- rbind()
#for(i in 1:14){
##	res_sgl <- rbind(res_sgl, matrix(res_sgen1[i,], 20, ry))
##	res_sgl <- rbind(res_sgl, matrix(res_sgl1[i,], 20, ry))
#	res_sm <- rbind(res_sm, matrix(Bbest[i,], 20, ry))
#}
#
#par(mfrow = c(1, 1))
#par(mar = c(2.6, 2.1, 3.1, 1.1))
#image.plot(1:nrow(res_sm), 1:ncol(res_sm), res_sm,
#		xaxt = 'n', yaxt = 'n', main = 'Results for SGEN',
#		zlim = range(c(res_sm)),
#		cex.main = 2,
#		cex.axis = 2)
#


