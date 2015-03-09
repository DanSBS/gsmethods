smtlasso <-
		function(x, y, grp, lambda, rho, ngam, max.iter,
				eps, method = 'smtl', intercept = TRUE,
				parallel = TRUE, pos = 0, singled = 0, 
				scale_lambda = FALSE, scale_rho = FALSE, 
				gam_eps = NULL, ridge_df = FALSE,
				B_start = NULL){
	
	p <- ncol(x)
	
	if(length(lambda) < p)
		lambda <- rep(lambda, p)
	if(intercept)
		lambda[1] <- 0
	
	
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
		D <- array(t(x) %*% x / nrow(x), dim = c(p, p, 1))
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
						pos, dslices), 
				dim = c(p, k, ngam))
	if(method == 'sgl'){
#		browser()
		Best <- array(sgl(C, D, B, max.iter, 
						lambda, gamma, ngam, eps, singled,
						pos, dslices), 
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
					
#					err <- try({chk <- colSums(candAlph) <= 7},TRUE)
					
#					if(!all(!chk)){
#						candAlph <- candAlph[,chk]
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
						
#					candAlpha <- candAlpha[,!chk]
						if(!is.null(dim(candB))){
							c0 <- abs(candB) > 0
#							browser()
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
							else
								npars <- colSums(abs(candB) > 0)
							rss <- colSums((predi - ym)^2)/(nrow(xm) - npars)
							
							## Using AIC or BIC here, haven't settled on one or the
							## other yet. The standard LASSO algorithm (or the one
							## implemented in R uses Cp, will have to think about whic
							## is most appropriate here).
#							score <- nrow(xm) * log(rss) + npars * log(nrow(xm))
#							score <- nrow(xm) * log(rss) + 2 * npars
							score <- (rss * (nrow(xm) - npars))/rss[ngam] - nrow(xm) + 2*npars						
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
#		browser()
		Bopt <- lapply(indx, FUN = function(u){
					
					candB <- Best[,u,]
					
#					err <- try({chk <- colSums(candAlph) <= 7},TRUE)
					
#					if(!all(!chk)){
#						candAlph <- candAlph[,chk]
					if(!is.null(dim(candB))){
						
#						if(!g_zero)
						candB <- (1+rho[u])*candB
						
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
#							score <- nrow(xm) * log(rss) + 2 * npars
							## Cp
							score <- (rss * (nrow(xm) - npars))/rss[ngam] - nrow(xm) + 2*npars
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
	return(Bbest)
	
	
}


sgen <- function(x, y, grp, rho, ngam, nlam, max.iter,
		eps, method = 'smtl', intercept = TRUE,
		parallel = TRUE, pos = 0, singled = 0, 
		scale_lambda = FALSE, scale_rho = FALSE, 
		gam_eps = NULL, ridge_df = FALSE, w = 1){
	
	p <- ncol(x)
	
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
#		Dlist <- lapply(names(sx), function(g){
#					gmat <- as.matrix(sx[[g]])
#					t(gmat) %*% gmat / nrow(gmat)
#				})
#		D <- abind(Dlist, along = 3)
		dslices <- dim(D)[3]
	}
	else{
		k <- ncol(y)
		C <- t(x) %*% y / nrow(x)
#		D <- array(t(x) %*% x / nrow(x), dim = c(p, p, 1))
		dslices <- 1
	}
	
	## Get range of values to search across for
	## group sparsity
	lam <- rowSums(abs(C))
	max_lam <- 1e-1 * max(lam)
	lam_seq <- exp(seq(log(max_lam), log(1e-4 * max_lam),, nlam))
	all_lam <- t(apply(cbind(lam), 1, function(u) exp(seq(log(u), log(u * 1e-6),, nlam))))
	
	mod_all <- vector('list', nlam)
	browser()
	for(i in 1:nlam){
		print(paste('Iteration:', i))
		if(i == 1)
			B_start <- NULL
		else
			B_start <- mod
		
		mod <- smtlasso(x, y, grp, lam_seq[i], rho, ngam, max.iter,
				eps, method, intercept,
				parallel, pos, singled, 
				scale_lambda, scale_rho, 
				gam_eps, ridge_df,
				B_start = B_start)
		
		mod_all[[i]] <- mod
	}
	
	N <- nrow(y)
	
	score <- unlist(lapply(mod_all, function(b){
#						rss_all <- colSums(y - x %*% b)^2/(N - colSums(b != 0))
#						bdf <- sum(rowSums(b) != 0)
#						bdf <- sum(b != 0)
#						bdf <- colSums(b != 0)
#						rss <- sum((y - x %*% b)^2)/(N - bdf)
#						sum(colSums((y - x %*% b)^2)/(N - bdf))
#						rss <- sum((y - x %*% b)^2)/(N - bdf)
#						sum(rowSums(b != 0)))
						bdf <- sum(rowSums(b) != 0)
#						bdf <- sum(b != 0)
						rss <- sum((y - x %*% b)^2)/(N - bdf)
						N * log(rss) + sum(bdf) * log(N)
#						N * log(sum(rss)) + sum(bdf) * log(N)
					}))
	
}