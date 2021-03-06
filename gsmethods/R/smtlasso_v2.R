require(Rcpp)
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")
sourceCpp('C:/Users/samarov/git/gsmethods/gsmethods/src/sgl_v2.cpp')
sourceCpp('C:/Users/samarov/git/gsmethods/gsmethods/src/smtl_v2.cpp')

smtlasso <- function(x, y, grp, lambda, rho, ngam, max.iter,
		eps, method = 'smtl', intercept = TRUE,
		parallel = TRUE, pos = 0, singled = 0, 
		scale_lambda = FALSE, scale_rho = FALSE, 
		gam_eps = NULL, ridge_df = FALSE,
		B_start = NULL, clus = TRUE, rho2 = 0,
		chk1 = FALSE){
	
	p <- ncol(x)
	
	if(length(lambda) < p)
		lambda <- rep(lambda, p)
	if(intercept)
		lambda[1] <- 0
	
#	browser()
	if(singled == 0){
		sx <- split(as.data.frame(x), grp)
		
		## The number of tasks
		k <- length(sx)
		sy <- split(y, grp)
		
		
		Clist <- lapply(names(sx), function(g){
					sxg <- as.matrix(sx[[g]])
					sxy <- sy[[g]]
					t(sxg) %*% sxy / nrow(sxg)
				})
		C <- do.call('cbind', Clist)
#		browser()
		Dlist <- lapply(names(sx), function(g){
					gmat <- as.matrix(sx[[g]])
					t(gmat) %*% gmat / nrow(gmat)
				})
		D <- abind(Dlist, along = 3)
		dslices <- dim(D)[3]
#		browser()
	}
	else{	
		k <- ncol(y)
		C <- t(x) %*% y / nrow(x)
		D <- array(t(x) %*% x / nrow(x) + diag(ncol(x)) * rho2, 
				dim = c(p, p, 1))
		dslices <- 1
	}
	if(method == 'smtl'){
#		lambda <- lambda/sqrt(rowSums(C^2))
#		rho <-  rep(rho, p)/sqrt(rowSums(C^2))
		if(scale_lambda)
			lambda <- lambda/rowSums(abs(C))
#		lambda <- lambda#/apply(abs(C),1,max)
#		lambda <- lambda/sqrt(rowSums(C^2))
		if(length(rho) < k){
			if(scale_rho)
				rho <- rep(rho, k)/sqrt(colSums(C^2))
			else
				rho <- rep(rho, k)
		}
		
		
#		if(intercept)
#			rho[1] <- 0
		
		rho[is.nan(rho)] <- 0
#		rho q`<- rep(rho, p)/sum(rowSums(abs(C)))
	}
	if(method == 'sgl'){
		if(scale_lambda)
			lambda <- lambda/sqrt(rowSums(C^2))	
		rho <- rep(0, k)
	}
	
	
	lambda[is.nan(lambda)] <- 0
	## Being lazy here
	if(ngam <= 1){
		ngam <- 2
		g_zero <- TRUE
	}
	else{
		g_zero <- FALSE
	}
	
	gamma <- apply(C, 2, function(u) {
				maxu <- max(abs(u))
				minu <- min(abs(u))
#				if(minu == 0)
#					minu <- 1e-10
#				exp(seq(log(1.05*maxu), log(1e-10),, ngam))
				exp(seq(log(maxu), log(1e-6 * maxu),, ngam))
#				exp(seq(log(maxu), log(1e-6),, ngam))
#				seq(1.05*maxu, 1e-10,, ngam)
#				exp(seq(log(maxu), log(min(minu, 1e-5)),, ngam))
#				seq(maxu, min(minu, 1e-5),, ngam)
			})
#	browser()
	
	## Check to see if beta has been provided
	
	B <- array(0, dim = c(p, k, ngam))
	if(!is.null(B_start)){
		B[,,1] <- B_start
	}
	
#	browser()
	if(!is.null(gam_eps) & !g_zero){
#		browser()
		ngam <- ceiling(mean(rowSums(abs(t(apply(gamma, 2, diff))) > gam_eps)))
		gamma <- gamma[1:ngam, ]
		B <- B[,,1:ngam]
#		ngam <- gam_thr
	}
	
	if(g_zero)
		gamma <- gamma * 0
	
	
#	browser()
	if(method == 'smtl')
		Best <- array(smtl(C, D, B, max.iter, rho, 
						lambda, gamma, ngam, eps, singled,
						pos, dslices, t(y), x), 
				dim = c(p, k, ngam))
	if(method == 'sgl'){
#		browser()
		Best <- array(sgl(C, D, B, max.iter, 
						lambda, gamma, ngam, eps, singled,
						pos, dslices, t(y), x), 
				dim = c(p, k, ngam))
	}
	
#	browser()
	indx <- 1:ncol(C)
	
#	browser()
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
					if(chk1)
						candB <- candB[,colSums(candB) <= 1]
#					candB <- t(t(candB)/colSums(candB))
#					err <- try({chk <- colSums(candAlph) <= 7},TRUE)
					
#					if(!all(!chk)){
#						candAlph <- candAlph[,chk]
					if(!is.null(dim(candB))){
						
						if(!g_zero)
							candB <- candB
#						rho
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
						
#					candAlpha <- candAlpha[,!chk]
						if(!is.null(dim(candB))){
							c0 <- abs(candB) > 0
#							browser()
#							if(!all(rho == 0) & ridge_df)
#								npars <- apply(c0, 2, function(v) {
#											xm0 <- xm[, v]
#											rho0 <- rho[u]
#											if(length(xm0) == 0)
#												0
#											else{
##										browser()
#												if(is.null(dim(xm0))){
#													xm0 <- t(rbind(xm0))
##													rho0 <- rbind(rho0)
#												}
#												
#												xtx0 <- ginv(t(xm0) %*% xm0 + rho0 * diag(ncol(xm0)))
#												sum(apply(xm0, 1, function(w) t(w) %*% xtx0 %*% w))
#											}
#										})
							
#								npars <- apply(c0, 2, function(v) {
#											xm0 <- xm[, v]
#											rho0 <- rho[v]
#											if(length(xm0) == 0)
#												0
#											else{
							##										browser()
#												if(is.null(dim(xm0))){
#													xm0 <- t(rbind(xm0))
#													rho0 <- rbind(rho0)
#												}
#												
#												xtx0 <- ginv(t(xm0) %*% xm0 + diag(rho0))
#												sum(apply(xm0, 1, function(w) t(w) %*% xtx0 %*% w))
#											}
#										})
#							else
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
#					}
#					else{
#						## If the criteria is not met return a vector of zeros.
#						rep(0, ncol(beta))
#					}
					
				})
	}
	else{
		if(clus){
			
			score <- unlist(lapply(1:ngam, function(i){
								grp <- apply(t(Best[,,i]), 1, function(u){
											order(u, decreasing = TRUE)[1]
										})
								
								SSW <- Reduce('+',lapply(split(as.data.frame(t(y)), grp), function(u){
													matu <- scale(as.matrix(u), scale = FALSE)
													sum(colSums(matu^2))
												}))
								my <- scale(t(y), scale = FALSE)
								SST <- sum(colSums(my^2))
								SSB <- SST - SSW
								log(SSB/SSW)
							}))
			o <- order(score, decreasing = TRUE)[1]
			Bopt <- Best[,,o]
		}
		else{
			Bopt <- lapply(indx, FUN = function(u){
#						browser()
						candB <- Best[,u,]
#						browser()
#						candB <- candB[,colSums(candB) <= 1]
						candB <- t(t(candB)/colSums(candB))
#					err <- try({chk <- colSums(candAlph) <= 7},TRUE)
#						apply(candB)
#					if(!all(!chk)){
#						candAlph <- candAlph[,chk]
						if(!is.null(dim(candB))){
							
#						if(!g_zero)
#							candB <- (1+rho[u])*candB
							
							if(singled == 0){
								xm <- as.matrix(sx[[u]])
								ym <- c(sy[[u]])
								predi <- as.matrix(sx[[u]]) %*% candB
							}
							else{
								predi <- as.matrix(x) %*% candB
								xm <- x
								ym <- y[,u]
							}
							
							
#					candAlpha <- candAlpha[,!chk]
							if(!is.null(dim(candB))){
								c0 <- abs(candB) > 0
								if(!all(rho == 0) & ridge_df)
									npars <- apply(c0, 2, function(v) {
												xm0 <- xm[, v]
												rho0 <- rho[u]
												if(length(xm0) == 0)
													0
												else{
#										browser()
													if(is.null(dim(xm0))){
														xm0 <- t(rbind(xm0))
#													rho0 <- rbind(rho0)
													}
													
													xtx0 <- ginv(t(xm0) %*% xm0 + rho0 * diag(ncol(xm0)))
													sum(apply(xm0, 1, function(w) t(w) %*% xtx0 %*% w))
												}
											})
								else
									npars <- colSums(abs(candB) > 0)
								rss <- colSums((predi - ym)^2)/(nrow(xm) - npars)
#							if(any(is.nan(log(rss))))
#								browser()
								## Using AIC or BIC here, haven't settled on one or the
								## other yet. The standard LASSO algorithm (or the one
								## implemented in R uses Cp, will have to think about whic
								## is most appropriate here).
								## BIC
#							score <- nrow(xm) * log(rss) + npars * log(nrow(xm))
								## AIC
#								browser()
								score <- nrow(xm) * log(rss) + 2 * npars
								## Cp
								##score <- (rss * (nrow(xm) - npars))/rss[ngam] - nrow(xm) + 2*npars
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
	}
	gc()
#	browser()
	## Collect the final results
	if(!clus)
		Bbest <- t(do.call("rbind", Bopt))
	else
		Bbest <- Bopt
	
	if(parallel)
		stopCluster(cl)
	return(Bbest)
	
	
}


rss <- function(x,y,b){
#	sum((y - t(b) %*% t(x))^2)
	sum(rowSums((y - t(b) %*% t(x))^2))
}

bicG <- function(x, y, b, thresh = 0.001){
#	DF <- sum(rowMeans(b != 0))
#	b[abs(b) <= thresh] <- 0
	nrow(x) * log(rss(x,y,b)/(nrow(x) - sum(rowSums(b) != 0))) +
#	nrow(y) * log(rss(x,y,b)/(nrow(y) - DF)) +
#			DF * log(nrow(y))
			sum(rowSums(b) != 0) * log(nrow(x))	
}

smtlassoBic <- function(X, Y, grp,
		rho = 0, ngam = 100,
		max.iter = 100,
		eps = 1e-3,
		method = 'smtl',
		intercept = FALSE,
		parallel = TRUE,
		scale_lambda = FALSE,
		pos = 1, singled = 1, gam_eps = 1e-10,
		ridge_df = FALSE,
		clus = FALSE,
		lambda_min = NULL,
		nlam = 10, rho2 = 0,
		chk1 = FALSE,
		minFrac = 1e-4,
		lambda_max = NULL){
	B <- t(Y) %*% X / nrow(X)
	
	if(is.null(lambda_max))
		lambda_max <- max(colSums(abs(B)))
	if(is.null(lambda_min))
		lambda_min <- lambda_max * minFrac
	
	lambda_seq <- exp(seq(log(lambda_min), log(lambda_max),, nlam))
#	lambda_seq <- lambda_seq[floor(seq(1,50,,10))]
#	nlam <- 10
#	lambda_seq <- lambda_seq[1:40][floor(seq(1,40,,20))]
#	nlam <- 20
#	lambda_seq <- lambda_seq[31:50]
	res_list <- vector('list', nlam)
	for(i in 1:nlam){
		print(paste('Iteration', i))
		res_list[[i]] <- smtlasso(X, Y, grp = grp, 
				lambda = lambda_seq[i],
				rho = rho,
				ngam = ngam, 
				max.iter = max.iter,
				eps = eps,
				method = method, intercept = intercept, 
				parallel = parallel,
				scale_lambda = scale_lambda,
				pos = pos, singled = singled, 
				gam_eps = gam_eps, 
				ridge_df = ridge_df, clus = clus,
				rho2 = rho2,
				chk1 = chk1
				)
		gc()
	}
	
#	browser()
	BIC <- unlist(lapply(res_list, function(b){
						bicG(X, t(Y), b)
					}))
	
	return(list(res_list = res_list, BIC = BIC))
}