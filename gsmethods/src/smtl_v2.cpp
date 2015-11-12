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

class FindMRes
{
public:
	unsigned int m;
	double tot;
	void set_values(const unsigned int&, const double&);
};

void FindMRes::set_values(const unsigned int& a, const double& b){
	m = a;
	tot = b;
}

FindMRes findM (const uvec& sIndx, const mat& aA,
		const double& lam, const unsigned int& k)
{

	unsigned int m;
	m = 1;
	double tot, tot0;
	tot = aA(sIndx(0)) - lam;
	tot0 = tot - 1.0;
	int ik;
	ik = (int)k - 1;

	for(int kk = 0; kk < ik; kk++){
		tot0 = tot;
		tot = ((double)m * tot + aA(sIndx(m))) / ((double)m + 1.0);
		if(tot > tot0)
			m += 1;
		else{
			tot = tot0;
			break;
		}
	}

	FindMRes out;
	out.set_values(m, tot);

	return out;
}

arma::uvec rankvec(const arma::uvec& sInd, const int& k){
	arma::uvec r(k);
	for(int i = 0; i < k; i++){
		r(sInd(i)) = i;
	}
	return r;
}

// [[Rcpp::export]]
SEXP smtl(SEXP R_C, SEXP R_D, SEXP R_B, const int& maxIter,
		SEXP R_lambda, SEXP R_lambda1, SEXP R_lambda2,
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
	NumericVector Rcpp_lambda(R_lambda);
	NumericVector Rcpp_lambda1(R_lambda1);
	NumericMatrix Rcpp_lambda2(R_lambda2);
	mat C(Rcpp_C.begin(), Rcpp_C.nrow(), Rcpp_C.ncol(), false);
	mat lambda2(Rcpp_lambda2.begin(), Rcpp_lambda2.nrow(), Rcpp_lambda2.ncol(), false);
	cube D(Rcpp_D.begin(), C.n_rows, C.n_rows, dslices, false);
	cube B(Rcpp_B.begin(), C.n_rows, C.n_cols, nl, false);

	NumericMatrix Rcpp_x(R_x);
	mat x(Rcpp_x.begin(), Rcpp_x.nrow(), Rcpp_x.ncol(), false);
	NumericMatrix Rcpp_y(R_y);
	mat y(Rcpp_y.begin(), Rcpp_y.nrow(), Rcpp_y.ncol(), false);
	// Different parameters w/ dimension information

	int k, p;
	k = C.n_cols;
	p = C.n_rows;

	//Ridge parameter
	//vec lambda(Rcpp_lambda.begin(), p);
	vec lambda(Rcpp_lambda.begin(), k);
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
	//arma::vec A(k);
	//A.fill(0);
	//arma::vec aA(k);
	//aA.fill(0);

	// Initialize temp vector
	double sA;
	//double tmpA;
	//double tmpA0;
	//double l2;
	double seps;
	double seps_old;
	//double den;
	//double s;
	seps = eps + 1;
	arma::uvec sInd(k);
	arma::uvec rInd(k);
	FindMRes mres;
	int iter;
	double frymat = sqrt(sum(sum(square(y))));
	//cout << lambda(0) << endl;
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

				// Initializing A
#pragma omp parallel
				{
#pragma omp for
					for(int c = 0; c < k; c++){
						double l2 = lambda2(nli, c);

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

						s = (double)sgn(tmpA0);

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
				}

				sA = sum(sum(abs(A)));

				if((double)sA <= lambda1(j)){
					B.slice(nli).row(j).fill(0);
				}
				else{
					vec aA = abs(A);

					sInd = sort_index(aA, 1);
					rInd = rankvec(sInd, k);
					mres = findM(sInd, aA, lambda1(j), k);

					// mres.m is the max 'm'
					// mres.tot is the associated total value
#pragma omp parallel
					{
#pragma omp for
						for(int c = 0; c < k; c++){
							int cc;
							if(singled == 1){
								cc = 0;
							}
							else{
								cc = sInd(cc);
							}
							mat Ds = D.slice(cc);
							double lam_r = lambda(0);
							if(rInd(sInd(c)) >= (mres.m - 1)){
								if((Ds(j,j) != 0) | (lam_r > 0))
									B(j,sInd(c),nli) = A(sInd(c))/(Ds(j,j)+lam_r);
								else
									B(j,sInd(c),nli) = 0.0;
							}
							else{
								if((Ds(j,j) != 0) | (lam_r > 0))
									B(j,sInd(c),nli) = mres.tot / (Ds(j,j)+lam_r);
								else
									B(j,sInd(c),nli) = 0.0;
							}
							double den = Bold(j,sInd(c));
							if(abs(den) > 1e-7){
								Bdiffs(j, sInd(c)) = (den - B(j,sInd(c),nli))/den;
							}
							else{
								Bdiffs(j, sInd(c)) = 0.0;
							}
							Bold(j,sInd(c)) = B(j,sInd(c),nli);

						}
					}
				}
			}

			//Bdiffs = abs(Bold - B.slice(nli));
			//seps = mean(mean(abs(Bdiffs)));
			seps = sqrt(sum(sum(square(y.t() - x * B.slice(nli))))) / frymat;
			if(abs(seps - seps_old) < eps){
				//cout << 11111 << endl;
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
