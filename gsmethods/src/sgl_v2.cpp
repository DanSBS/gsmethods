#include <Rcpp.h>
//#include </home/samarov/cpplibs/armadillo-3.800.2/include/armadillo>
#include <C:\cpplibs\armadillo-4.320.0\include\armadillo>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <omp.h>
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
		const int& dslices,
		SEXP R_y, SEXP R_x){
	omp_set_num_threads( 8 );
	NumericMatrix Rcpp_C(R_C);
	NumericVector Rcpp_D(R_D);
	NumericVector Rcpp_B(R_B);
	NumericVector Rcpp_lambda1(R_lambda1);
	NumericMatrix Rcpp_lambda2(R_lambda2);

	NumericMatrix Rcpp_x(R_x);
	mat x(Rcpp_x.begin(), Rcpp_x.nrow(), Rcpp_x.ncol(), false);
	NumericMatrix Rcpp_y(R_y);
	mat y(Rcpp_y.begin(), Rcpp_y.nrow(), Rcpp_y.ncol(), false);

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
	//mat Ds(p,p);

	// Initialize matrix A
	arma::vec A(k);
	A.fill(0);

	// Initialize temp vector
	double sA;
	//double tmpA;
	//double tmpA0;
	//double l2;
	double seps;
	double seps_old;
	//double den;
	//double s;
	//double b2;
	seps = eps + 1;
	int iter;
	mat betai(1, p);
	double frymat = sqrt(sum(sum(square(y))));
	//mat Ds = D.slice(0);

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
			for(int j = 0; j < p; j++){
				if((sum(sum(abs(B.slice(nli).row(j)))) <= 1e-10) & (iter > 1)){
					continue;
				}
				vec A(k);
				A.fill(0);
				betai = Bold.row(j);
				double bden = sqrt(sum(sum(square(betai))));

#pragma omp parallel
				{
#pragma omp for
					for(int c = 0; c < k; c++){
						if((abs(B(j,c,nli)) <= 1e-10) & (iter > 1)){
							continue;
						}
						double l2 = lambda2(nli, c);

						//if(singled == 1){
						//	mat Ds = D.slice(0);
						//}
						//else{
						//	mat Ds = D.slice(c);
						//}
						int cc;
						if(singled == 1){
							cc = 0;
						}
						else{
							cc = c;
						}
						mat Ds = D.slice(cc);

						double s;
						double tmpA0 = C(j, c) - dot( Ds.row(j), Bold.col(c)) +
								Ds(j, j) * Bold(j, c);
						double tmpA = abs(tmpA0) - l2;

						if(tmpA0 < 0)
							s = -1.0;
						else
							s = 1.0;
						if((pos == 1) & (s < 0)){
							s = 0.0;
						}

						if(tmpA > 0){
							A(c) = s * tmpA;
						}
						else{
							A(c) = 0.0;
						}


						if(bden > 1e-7){
							if((Ds(j,j) != 0) | (lambda1(j) > 0))
								B(j,c,nli) = A(c)/(Ds(j,j)+lambda1(j)/bden);
							else
								B(j,c,nli) = 0.0;
						}
						else{
							bden = 1.0;
							if(iter > 1){
								B(j, c, nli) = 0.0;
							}
							else{
								if((Ds(j,j) != 0) | (lambda1(j) > 0))
									B(j,c,nli) = A(c)/(Ds(j,j)+lambda1(j)/bden);
								else
									B(j,c,nli) = 0.0;
							}
						}

						double den = Bold(j,c);
						if(abs(den) > 1e-7){
							Bdiffs(j, c) = (den - B(j,c,nli))/den;
						}
						else{
							Bdiffs(j, c) = 0.0;
						}
						Bold(j,c) = B(j,c,nli);

					}
				}

				sA = sqrt(sum(sum(square(A))));

				if((double)sA <= lambda1(j)){
					B.slice(nli).row(j).fill(0);
				}
			}
			seps = sqrt(sum(sum(square(y.t() - x * B.slice(nli))))) / frymat;
			if(abs(seps - seps_old) < eps){
				//cout << seps << endl;
				seps = eps;
			}
			seps_old = seps;

		}

	}

	// Next we need find the best lambda according the Mallow's Cp (this isn't actually
	// the only choice, but we use it here for now)
	// cout << B.subcube(span(1),span(),span()).n_slices << endl;
	return Rcpp::wrap(B);
}
