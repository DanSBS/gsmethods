#include <Rcpp.h>
//#include </home/samarov/cpplibs/armadillo-3.800.2/include/armadillo>
#include <C:\cpplibs\armadillo-4.320.0\include\armadillo>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
using namespace std;
using namespace Rcpp;
//using namespace RcppArmadillo;
using namespace arma;

/*
 * This code performs the multitask lasso.
 */

template <typename T>
double sgn(const T& val) {
	return (T(0) < val) - (val < T(0));
}

// [[Rcpp::export]]
SEXP sgl(SEXP R_C, SEXP R_D, SEXP R_B, const int& maxIter,
		SEXP R_lambda1, SEXP R_lambda2,
		const int& nl,
		const double& eps,
		const int& singled,
		const int& pos,
		const int& dslices){

	NumericMatrix Rcpp_C(R_C);
	NumericVector Rcpp_D(R_D);
	NumericVector Rcpp_B(R_B);
	NumericVector Rcpp_lambda1(R_lambda1);
	NumericMatrix Rcpp_lambda2(R_lambda2);
	mat C(Rcpp_C.begin(), Rcpp_C.nrow(), Rcpp_C.ncol(), false);
	mat lambda2(Rcpp_lambda2.begin(), Rcpp_lambda2.nrow(), Rcpp_lambda2.ncol(), false);
	cube D(Rcpp_D.begin(), C.n_rows, C.n_rows, dslices, false);
	cube B(Rcpp_B.begin(), C.n_rows, C.n_cols, nl, false);

	// Different parameters w/ dimension information
	int k, p;
	k = C.n_cols;
	p = C.n_rows;

	//Group sparse parameter
	vec lambda1(Rcpp_lambda1.begin(), p);

	// Initial values for coefficients B. So each row here corresponds to
	// an observed, multivariate Y_k = n x 1 vector. In the main loop of
	// the algorithm we are going to be looping through the columns which
	// correspond to the contribution of the j^th coefficient B_j to
	// the observed outcome.
	arma::mat Bold(p,k);
	Bold = B.slice(0);
	arma::mat Bdiffs(p,k);
	arma::mat Ds(p,p);

	// Initialize matrix A
	arma::vec A(k);
	A.fill(0);

	// Initialize temp vector
	double sA;
	double tmpA;
	double tmpA0;
	double l2;
	double seps;
	double seps_old;
	double den;
	double s;
	double b2;
	seps = eps + 1;
	int iter;

	//================== Beginning of main loop ====================
	for(int nli = 0; nli < nl; nli++){
		// Get new lamdba

		// Reset check for convergence
		seps = eps + 1;
		seps_old = seps;
		iter = 0;
		if(nli > 0){
			Bold = B.slice(nli-1);
		}
		while((iter < maxIter) & (seps > eps)){
			iter += 1;
			//cout << 777777 << endl;
			for(int j = 0; j < p; j++){
				if((sum(sum(abs(B.slice(nli).row(j)))) <= 1e-10) & (iter > 1)){
					//cout << 999999 << endl;
					continue;
				}
				// Initializing A
				//cout << 8888 << endl;
				for(int c = 0; c < k; c++){
					l2 = lambda2(nli, c);

					if(singled == 1){
						Ds = D.slice(0);
					}
					else{
						Ds = D.slice(c);
					}

					tmpA0 = C(j, c) - dot( Ds.row(j), Bold.col(c)) +
							Ds(j, j) * Bold(j, c);
					tmpA = abs(tmpA0) - l2;
					//cout << tmpA0 << endl;

					if(tmpA0 < 0)
						s = -1.0;
					else
						s = 1.0;
					//s = (double)sgn(tmpA0);
					//cout << s << endl;
					if((pos == 1) & (s < 0)){
						s = 0.0;
					}

					if(tmpA > 0){
						A(c) = s * tmpA;
					}
					else{
						A(c) = 0.0;
					}
				}

				sA = sqrt(sum(sum(square(A))));
				//cout << sA << endl;
				if((double)sA <= lambda1(j)){
					B.slice(nli).row(j).fill(0);
				}
				else{

					for(int kk = 0; kk < k; kk++){
						if(singled == 1){
							Ds = D.slice(0);
						}
						else{
							Ds = D.slice(kk);
						}

						/*if(iter == 1)
							b2 = 1;
						else*/
						//b2 = sqrt(sum(sum(pow(B.slice(nli).row(j),2))));
						b2 = sqrt(sum(sum(square(A))));
						if(b2 > 1e-7){
							if((Ds(j,j) != 0) | (lambda1(j) > 0))
								B(j,kk,nli) = A(kk)/(Ds(j,j)+lambda1(j)/b2);
							else
								B(j,kk,nli) = 0.0;
						}
						else{
							b2 = 1.0;
							if(iter > 1){
								B(j, kk, nli) = 0.0;
							}
							else{
								if((Ds(j,j) != 0) | (lambda1(j) > 0))
									B(j,kk,nli) = A(kk)/(Ds(j,j)+lambda1(j)/b2);
								else
									B(j,kk,nli) = 0.0;
								//B(j,kk,nli) = A(kk)/(Ds(j,j)+lambda1(j)/b2);
							}
						}

						den = Bold(j,kk);
						if(abs(den) > 1e-7){
							Bdiffs(j, kk) = (den - B(j,kk,nli))/den;
						}
						else{
							Bdiffs(j, kk) = 0.0;
						}
						Bold(j,kk) = B(j,kk,nli);

					}
				}
			}

			//Bdiffs = abs(Bold - B.slice(nli));
			seps = mean(mean(abs(Bdiffs)));
			if(abs(seps - seps_old) < eps){
				//cout << 11111 << endl;
				seps = eps;
			}
			seps_old = seps;
		}

		//B.slice(nli) = (1 + lambda) * B.slice(nli);
	}

	// Next we need find the best lambda according the Mallow's Cp (this isn't actually
	// the only choice, but we use it here for now)
	// cout << B.subcube(span(1),span(),span()).n_slices << endl;
	return Rcpp::wrap(B);
}
