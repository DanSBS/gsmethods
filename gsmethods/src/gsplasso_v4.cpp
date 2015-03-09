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

class spmodel
{
public:
	mat gam;
	mat spindx;
	void set_values(mat&, mat&);
};

void spmodel::set_values(mat& a, mat& b){
	gam = a;
	spindx = b;
}

vec mdist(cube x, cube xaug,
		int nr,
		int nc,
		int ns){

	mat xmat(nr * nc, ns);
	mat xaugmat(nr * nc, ns);
	mat xdiff2(nr * nc, ns);
	mat x2(nr * nc, ns);
	mat xmi(nr * nc, 1);
	mat xmai(nr * nc, 1);
	mat xdi(nr * nc, 1);
	int n;
	double mx;
	double mxi;
	double w;
	vec dist(nr * nc);
	vec edistsc(nr * nc);

	n = nr * nc;

	xmat = reshape(mat(x.memptr(), x.n_elem, 1, false), n, ns);
	mx = max(max(xmat));
	xaugmat = reshape(mat(xaug.memptr(), xaug.n_elem, 1, false), n, ns);

	xdiff2 = pow(xmat - xaugmat, 2);

	for(int i = 0; i < n; i++){
		xmi = xmat.row(i);
		xmai = xaugmat.row(i);
		xdi = xdiff2.row(i);
		mxi = max(max(xmai));
		w = (mx - mxi)/mxi;
		rowvec x2tmp(ns);


		x2tmp.elem(find((xmi > xmai) == (xmai > 0))) =
				w * xdi.elem(find((xmi > xmai) == (xmai > 0)));
		x2tmp.elem(find((xmi > xmai) == (xmai <= 0))) =
				(1/w)*xdi.elem(find((xmi > xmai) == (xmai <= 0)));

		x2tmp.elem(find((xmi <= xmai) == (xmai > 0))) =
				xdi.elem(find((xmi <= xmai) == (xmai > 0)));
		x2tmp.elem(find((xmi <= xmai) == (xmai <= 0))) =
				xdi.elem(find((xmi <= xmai) == (xmai <= 0)));

		x2.row(i) = x2tmp;
	}

	dist = exp(-1 * sum(x2, 1));
	edistsc = dist/max(max(dist));

	return edistsc;

}

spmodel spweights(cube x,
		mat w,
		const int& nr,
		const int& nc,
		const int& ns,
		const int& k,
		const double& thresh){
	//int cores;
	//cores = 8;
	//omp_set_num_threads(cores);
	// Matrix that will contain spatial weights
	mat gam(nr * nc, pow(2 * k + 1, 2));
	gam.fill(0);

	// Augmented matrix
	cube xaug(nr + 2 * k, nc + 2 * k, ns);
	xaug.fill(0);

	// Get denominator values for angle
	vec xnorm = sqrt(sum(pow(reshape(mat(x.memptr(), x.n_elem, 1, false), nr * nc, ns),2), 1));

	// Augmenting with neighborhood k
	xaug.tube(span(k, nr + k - 1), span(k, nc + k - 1)) = x;

	// Temp cube for generating neighborhood information
	cube xtmp(nr, nc, ns);
	cube xtxtmp(nr, nc, ns);

	// Vectors for computing the angle
	vec num(nr * nc);
	vec den(nr * nc);
	vec rat(nr * nc);

	// Index keeping track of neighbohood
	int indx = 0;

	// Neighborhood index
	mat nbraug(nr + 2 * k, nc + 2 * k);
	nbraug.fill(-1);
	mat nbr(nr, nc);
	int ii = 0;
	for(int i = 0; i < nc; i++){
		for(int j = 0; j < nr; j++){
			nbr(j, i) = ii;
			ii += 1;
		}
	}

	nbraug.submat(span(k, nr + k - 1), span(k, nc + k - 1)) = nbr;
	cube spindx(nr, nc, pow(2 * k + 1, 2));
	spindx.fill(0);
	mat smat(nr, nc);
	//#pragma omp parallel for schedule(static)
	for(int j = 0; j < (2 * k + 1); j++){
		for(int i = 0; i < (2 * k + 1); i++){

			// Grabbing i,j'th neighor information
			xtmp = xaug.tube(span(i, nr + i - 1), span(j, nc + j - 1));
			// Computation for numerator


			xtxtmp = x % xtmp;
			num = sum(reshape(mat(xtxtmp.memptr(), xtxtmp.n_elem, 1, false), nr * nc, ns), 1);

			// Denominator
			den = sqrt(sum(pow(reshape(mat(xtmp.memptr(), xtmp.n_elem, 1, false), nr * nc, ns), 2), 1)) % xnorm;

			// Angle
			rat = num/den;

			// Remove den values w/ 0
			rat.elem(find(den == 0)).zeros();
			// Remove angles less than threshold

			// Set to 0 those elements whose angle is too large
			rat.elem(find(rat < thresh)).zeros();

			//rat = mdist(x, xtmp, nr, nc, ns) * w(i, j);
			rat.elem(find(rat < thresh)).zeros();
			// Mult. spatial weights
			gam.col(indx) = rat * w(i, j);
			//gam.col(indx) = rat;


			// Spatial index
			smat = nbraug.submat(span(i, nr + i - 1), span(j, nc + j - 1));
			spindx.slice(indx) = smat;

			// Iterate 1
			indx += 1;
		}
	}

	mat spindxmat = reshape(mat(spindx.memptr(), spindx.n_elem, 1, false), nr * nc, pow(2 * k + 1, 2));

	spmodel out;
	out.set_values(gam, spindxmat);

	return out;

}


// [[Rcpp::export]]
SEXP gsplasso(SEXP R_x, SEXP R_y,
		SEXP R_B,
		const int& k,
		const int& pos,
		SEXP R_lams,
		SEXP R_lamg,
		const double& rho,
		SEXP R_w,
		const int& nls,
		const int& nr,
		const int& nc,
		const int& ns,
		const double& eps,
		const int& maxIter,
		const double& thresh,
		SEXP R_ys){
	//int cores;
	//cores = 8;
	//omp_set_num_threads(cores);
	//=================================================
	// x: Design matrix
	// y: Cube of responses
	// k: spatial neighborhood to use
	// pos: Whether positivity constraints should be imposed
	// lams: Element-wise penalty
	// lamg: Group penalty
	// rho: Spatial penalty
	// w: Weight matrix
	// nls: Number of lambda values for the sparsity component
	omp_set_num_threads( 8 );
	// Getting data structured
	// Design matrix (e.g. endmembers)
	NumericMatrix Rcpp_x(R_x);
	mat x(Rcpp_x.begin(), Rcpp_x.nrow(), Rcpp_x.ncol(), false);
	// Hyperspectral data cube
	NumericVector Rcpp_y(R_y);
	cube y(Rcpp_y.begin(), nr, nc, ns, false);
	// Convert y to a matrix
	mat ymat = reshape(mat(y.memptr(), y.n_elem, 1, false), nr * nc, ns);
	// Spatial decay matrix
	NumericMatrix Rcpp_w(R_w);
	mat w(Rcpp_w.begin(), Rcpp_w.nrow(), Rcpp_w.ncol(), false);
	// Cross-covariance
	mat C = x.t() * ymat.t() / x.n_rows;
	//cout << max(max(C)) << endl;
	//cout << min(min(C)) << endl;
	// Covariance
	mat D = x.t() * x / x.n_rows;
	// Matrix containing coefficients
	NumericVector Rcpp_B(R_B);
	cube B(Rcpp_B.begin(), C.n_rows, C.n_cols, nls, false);
	// Smoothed Hyperspectral data cube
	NumericVector Rcpp_ys(R_ys);
	cube ys(Rcpp_ys.begin(), nr, nc, ns, false);

	// Sparse regularization parameters, cube
	NumericMatrix Rcpp_lams(R_lams);
	mat lams(Rcpp_lams.begin(), Rcpp_lams.nrow(), Rcpp_lams.ncol(), false);



	// Group regularization parameter
	NumericVector Rcpp_lamg(R_lamg);
	vec lamg(Rcpp_lamg.begin(), C.n_rows);

	// Grabbing dimensional parameters
	int n, p, pk;
	n = C.n_cols;
	p = C.n_rows;
	pk =  pow(2 * k + 1, 2);

	// Object for storing spatial weight and indeces
	spmodel spw;

	// Get spatial weights and index
	spw = spweights(ys, w, nr, nc, ns, k, thresh);
	mat gam = spw.gam;
	mat spindx = spw.spindx;

	// Matrices for storing coefficient results during iterations
	mat Bold(p,n);
	Bold = B.slice(0);
	mat Bdiffs(p,n);

	double sA;
	//double tmpA;
	//double tmpA0;
	//double ls;
	double seps;
	double seps_old;
	//double den;
	//double s;
	//double bden;
	int iter;

	mat gami(1, gam.n_cols);
	mat wtszero(1,n);
	//uvec rnd(p);
	double frymat = sqrt(sum(sum(square(ymat))));
	mat spindxi(1, spindx.n_cols);
	mat betai(1, p);
	//int spiwk;
	//int j;

	for(int nli = 0; nli < nls; nli++){
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
			vec rord = randu<vec>(p);
			uvec jind = sort_index(rord);
			for(int jj = 0; jj < p; jj++){
				int j = jind(jj);
				if((sum(sum(abs(B.slice(nli).row(j)))) <= 1e-10) & (iter > 1)){
					continue;
				}
				vec A(n);
				A.fill(0);
				wtszero.fill(0);
				betai = Bold.row(j);
				double bden = sqrt(sum(sum(square(betai))));
#pragma omp parallel
				{

#pragma omp for
					for(int c = 0; c < n; c++){
						if((abs(B(j,c,nli)) <= 1e-10) & (iter > 1)){
							continue;
						}
						double ls = lams(nli, c);
						mat gami = gam.row(c);
						mat spindxi = spindx.row(c);
						//betai = Bold.row(j);
						double wta, wtg, s;
						wta = 0.0;
						wtg = 0.0;

						for(int wk = 0; wk < pk; wk++){
							int spiwk = spindxi(wk);
							if(spiwk >= 0){
								wta += gami(wk) * betai(spiwk);
								wtg += gami(wk);
							}
						}

						if(wta > 0.0){
							wta = wta/wtg;
							wtszero(c) = 1.0;
						}

						double tmpA0 = C(j, c) - dot( D.row(j), Bold.col(c)) +
								D(j, j) * Bold(j, c) + rho * wta;
						double tmpA = abs(tmpA0) - ls;

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
							if((D(j,j) != 0) | (lamg(j) > 0))
								B(j,c,nli) = A(c)/(D(j,j)+lamg(j)/bden + wtszero(c)*rho);
							else
								B(j,c,nli) = 0.0;
						}
						else{
							bden = 1.0;
							if(iter > 1){
								B(j, c, nli) = 0.0;
							}
							else{
								if((D(j,j) != 0) | (lamg(j) > 0))
									B(j,c,nli) = A(c)/(D(j,j)+lamg(j)/bden + wtszero(c)*rho);
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
				if((double)sA <= lamg(j)){
					B.slice(nli).row(j).fill(0);
				}
			}

			//seps = mean(mean(abs(Bdiffs)));
			seps = sqrt(sum(sum(square(ymat.t() - x * B.slice(nli))))) / frymat;
			if(abs(seps - seps_old) < eps){
				seps = eps;
			}

			seps_old = seps;
		}
	}

	return Rcpp::wrap(B);

}
