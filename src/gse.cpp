#include "gse.h"
using namespace Rcpp;
using namespace arma;
using namespace std;


/***************************************************/
/*               Function prototypes               */
/***************************************************/
vec rho1(vec x);
double rho1p(double x);
double solve_scales( vec maj, vec cc1, double tol, int miter );
mat pmd_adj(mat x, mat sigma0, mat sigmak, bool equalsig, mat mu, umat miss_group_unique, uvec miss_group_counts, vec tuning_const_group);
double scales(mat x, mat sigma0, mat sigmak, bool equalsig, mat mu, umat miss_group_unique, uvec miss_group_counts, vec tuning_const_group, double tol, int miter);
mat iterS( mat x, mat sigma0, mat sigmak, bool equalsig, mat mu, double sk, umat miss_group_unique, uvec miss_group_counts, vec muning_const_group);
SEXP GSE(SEXP X, SEXP N, SEXP P, SEXP Mu0, SEXP S0, SEXP Tol, SEXP Maxiter, SEXP Tol_scale, SEXP Miter_scale, SEXP Miss_group_unique, SEXP Miss_group_counts, SEXP Tuning_const_group, SEXP Print_step);


/***************************************************/
/*                Main GSE function                */
/***************************************************/
SEXP GSE(SEXP X, SEXP N, SEXP P, SEXP Mu0, SEXP S0, SEXP Tol, SEXP Maxiter, SEXP Tol_scale, SEXP Miter_scale, SEXP Miss_group_unique, SEXP Miss_group_counts, SEXP Tuning_const_group, SEXP Print_step)
{
	try{	
	// Input argument
	mat x = as<mat>(X);
	int p = as<int>(P);
	mat mu0 = as<mat>(Mu0);
	mat Sigma0 = as<mat>(S0);
	double tol = as<double>(Tol);
	int maxiter = as<int>(Maxiter);
	double tol_scale = as<double>(Tol_scale);
	int miter_scale = as<int>(Miter_scale);
	umat miss_group_unique = as<umat>(Miss_group_unique);
	uvec miss_group_counts = as<uvec>(Miss_group_counts);
	vec tuning_const_group = as<vec>(Tuning_const_group);
	int print_step = as<int>(Print_step);

	// Error code
	int error_code = 0; // 0 = no error, 1 = non-positive scale, 2 = non-positive definite est scatter
	
	// Basic setup
	double stilde0 = scales(x, Sigma0, Sigma0, true, mu0, miss_group_unique, miss_group_counts, tuning_const_group, tol_scale, miter_scale);
	double ep = 1.0;
	int iter = 0;

	// Proceed only if no error
	if( stilde0 <= 0) error_code = 1; 
	if( error_code == 0){
		Sigma0 = stilde0 * Sigma0;
		mat Omega = Sigma0;
		stilde0 = 1; 
		
		// Start iteration
		mat mu1(1,p);
		mat Sigma1(p,p);
		do
		{
			iter++;
			mat iter_res = iterS(x, Omega, Sigma0, false, mu0, stilde0, miss_group_unique, miss_group_counts, tuning_const_group);
			mu1 = iter_res.row(0);
			iter_res.shed_row(0);
			mat Stilde1 = iter_res;
			
			// Ensure symmetric
			Stilde1 = symmatl(Stilde1);
			
			// Check to see if it's positive definite
			mat Stilde1_chol(p,p);
			bool error_code_chol = chol(Stilde1_chol, Stilde1);
			if( !error_code_chol ) error_code = 2; 
			
			if( error_code == 0){
				double s1 = scales(x,Stilde1,Stilde1,true,mu1,miss_group_unique,miss_group_counts,tuning_const_group, tol_scale, miter_scale);
				Sigma1 = s1 * Stilde1;	

				double stilde1 = scales(x,Omega,Stilde1,false,mu1,miss_group_unique,miss_group_counts,tuning_const_group, tol_scale, miter_scale);
				ep = fabs(1 - stilde1/stilde0 );

				// Print out intermediate steps if necessary
				if( print_step == 2 ) Rcout << "iter=" << iter << "; tol=" << ep << "; scale=" << stilde1 << std::endl;

				// Updates
				Sigma0 = Sigma1;
				mu0 = mu1;
				stilde0 = stilde1;
			}
		}
		while( (ep > tol) && (error_code == 0) && (iter <= maxiter) );
	}
	
	return List::create( Named("S")=Sigma0, Named("mu")=mu0, Named("stilde0")=stilde0,
		Named("iter")=iter, Named("ep")=ep, Named("error_code")=error_code);

	} catch( std::exception& __ex__ ){
		forward_exception_to_r( __ex__ );
	} catch(...){
		::Rf_error( "c++ exception " "(unknown reason)" );
	}
	
}






/***************************************************/
/*               Iterative step                    */
/***************************************************/
mat iterS( mat x, mat sigma0, mat sigmak, bool equalsig, mat mu, double sk, umat miss_group_unique, uvec miss_group_counts, vec tuning_const_group)
{ 
	try{		
		/*********************************/	
		/* Some declaration of variables */
		int n_counts = miss_group_unique.n_rows;
		//int n = x.n_rows;
		unsigned int p = x.n_cols;
	
		uvec pp_grp = sum(miss_group_unique, 1);  // vector of observed variables for each missingness group
		int rowid_start = 0;

		double sweight1 = 0.0;				// sum of weight 1, w
		double sweight2 = 0.0;				// sum of weight 2, w X w*

		mat mu1 = zeros<mat>(1,p);			// mu_(k+1)
		mat Sigma1 = zeros<mat>(p,p);			// Sigma^tilde_(k+1)

		// things help to cross check (can be deleted)
		//mat xpredlist(n, p);				// n X p matrix of best predicted values 

		// Start computing
		for(int i = 0; i < n_counts; i++){
			mat mu_nonmiss(1, pp_grp(i));
			mat sigma0_nonmiss( pp_grp(i), pp_grp(i) );
			mat sigmak_nonmiss( pp_grp(i), pp_grp(i) );
			mat sigmaku_nonmiss( p, pp_grp(i) );
			mat xi( miss_group_counts(i) , pp_grp(i) );
			mat xpred(1, p);
			mat Ck = sigmak;
			int rowid_end = rowid_start + miss_group_counts(i) - 1;
			if( pp_grp(i) < p ){
				int mm = 0;
				for(unsigned int j=0; j<p; j++){
					int nn=mm;
					if(miss_group_unique(i,j) == 1){
						for(unsigned int k=j; k<p; k++){
							if( miss_group_unique(i,k) == 1 ){
								sigma0_nonmiss(mm, nn) = sigma0(j,k);
								sigma0_nonmiss(nn, mm) = sigma0(k,j);
								sigmak_nonmiss(mm, nn) = sigmak(j,k);
								sigmak_nonmiss(nn, mm) = sigmak(k,j);
								nn++;
							}
						}
						xi.col(mm) = x( span(rowid_start, rowid_end ), j );
						sigmaku_nonmiss.col(mm) = sigmak.col(j);	
						mu_nonmiss(0, mm) = mu(0, j);
						mm++;
					}
				}
			} else{
				sigmak_nonmiss = sigmak;
				sigma0_nonmiss = sigma0;
				xi = x.rows(rowid_start, rowid_end);
				mu_nonmiss = mu;
				//xpredlist.rows(rowid_start, rowid_end) = x.rows(rowid_start, rowid_end);
			}

			mat A = ones<mat>( pp_grp(i), pp_grp(i));
			mat diagA = diagmat(A);
			mat sigmak_nonmiss_inv = solve( sigmak_nonmiss, diagA );
			mat betas = sigmaku_nonmiss * sigmak_nonmiss_inv;
			if( pp_grp(i) < p){
				Ck = sigmaku_nonmiss * sigmak_nonmiss_inv * trans(sigmaku_nonmiss);
			}
			Ck = sigmak - Ck;

			double ppi = (double) pp_grp(i);
			
			// Added May 10, 2012, if fixed pt, no need to calculate dee or de0, as they are the same and ratio is always 1
			double dee = 1, de0 = 1, dee_ratio = 1;
			if( !equalsig ){
				dee = det( sigmak_nonmiss );
				dee = pow(dee, 1/ppi );
				de0 = det( sigma0_nonmiss );
				de0 = pow(de0, 1/ppi );		
				dee_ratio = dee / de0;
			}
			
			for(unsigned int m = 0; m < miss_group_counts(i); m++){
				mat xii = xi.row(m) - mu_nonmiss;
				double md = as_scalar(xii * sigmak_nonmiss_inv * trans(xii));
				double md2 = md*dee_ratio / (tuning_const_group(i) * sk );
				double weight1 = rho1p( md2 ) * dee_ratio;
				double weight2 = weight1 * (md / ppi);

				sweight1 += weight1;
				sweight2 += weight2;
				if( pp_grp(i) < p ){
					xpred = trans(betas * trans(xii)) + mu;
				} else{
					xpred = x.row(rowid_start + m);					
				}
				mu1 += weight1 * xpred;
				xpred -= mu;
				Sigma1 += weight1 * (trans(xpred) * xpred) + weight2 * Ck;
			}
			rowid_start = rowid_start + miss_group_counts(i);
		}

		mu1 = mu1/sweight1;
		Sigma1 = Sigma1/sweight2;

		mat res = join_cols(mu1, Sigma1);
		
		return res; 
	} catch( std::exception& __ex__ ){
		forward_exception_to_r( __ex__ );
	} catch(...){
		::Rf_error( "c++ exception " "(unknown reason)" );
	}
}


/*****************************************************************/
/* Compute scales or S* for given scaled md and tuning constants */
/*****************************************************************/
double scales(mat x, mat sigma0, mat sigmak, bool equalsig, mat mu, umat miss_group_unique, uvec miss_group_counts, vec tuning_const_group, double tol, int miter)
{
	try{
		mat res = pmd_adj(x, sigma0, sigmak, equalsig, mu, miss_group_unique, miss_group_counts, tuning_const_group);
		vec partialVec = res.col(0);
		vec cc1 = res.col(1);
		double s0 = solve_scales( partialVec, cc1, tol, miter );
		return s0; 
	} catch( std::exception& __ex__ ){
		forward_exception_to_r( __ex__ );
	} catch(...){
		::Rf_error( "c++ exception " "(unknown reason)" );
	}
}



mat pmd_adj(mat x, mat sigma0, mat sigmak, bool equalsig, mat mu, umat miss_group_unique, uvec miss_group_counts, vec tuning_const_group)
{
	try{
		int n_counts = miss_group_unique.n_rows;
		int n = x.n_rows;
		unsigned int p = x.n_cols;

		vec partialVec(n);
		vec cc1(n);
		uvec pp = sum(miss_group_unique, 1);
		int rowid_start = 0;

		for(int i = 0; i < n_counts; i++){
			mat mu_nonmiss(1, pp(i));
			mat sigma0_nonmiss( pp(i), pp(i) );
			mat sigmak_nonmiss( pp(i), pp(i) );
			mat xi( miss_group_counts(i) , pp(i) );
			int rowid_end = rowid_start + miss_group_counts(i) - 1;
			if( pp(i) < p ){
				int mm = 0;
				for(unsigned int j=0; j<p; j++){
					int nn=mm;
					if(miss_group_unique(i,j) == 1){
						for(unsigned int k=j; k<p; k++){
							if( miss_group_unique(i,k) == 1 ){
								sigma0_nonmiss(mm, nn) = sigma0(j,k);
								sigma0_nonmiss(nn, mm) = sigma0(k,j);
								sigmak_nonmiss(mm, nn) = sigmak(j,k);
								sigmak_nonmiss(nn, mm) = sigmak(k,j);
								nn++;
							}
						}
						xi.col(mm) = x( span(rowid_start, rowid_end ), j );
						mu_nonmiss(0, mm) = mu(0, j);
						mm++;
					}
				}
			} else{
				sigmak_nonmiss = sigmak;
				sigma0_nonmiss = sigma0;
				xi = x.rows(rowid_start, rowid_end);
				mu_nonmiss = mu;
			}

			mat A = ones<mat>( pp(i), pp(i));
			mat diagA = diagmat(A);
			mat sigmak_nonmiss_inv = solve( sigmak_nonmiss, diagA );

			double ppi = (double) pp(i);
			double aabb_ratio = 1;
			if( !equalsig ){
				double aa = det( sigmak_nonmiss );
				aa = pow( aa, 1/ppi );
				double bb = det( sigma0_nonmiss );
				bb = pow( bb, 1/ppi );
				aabb_ratio = aa / bb;
			}
			
			for(unsigned int m = 0; m < miss_group_counts(i); m++){
				mat xii = xi.row(m) - mu_nonmiss;
				double md = as_scalar(xii * sigmak_nonmiss_inv * trans(xii));
				partialVec(rowid_start + m) = md * aabb_ratio / tuning_const_group(i);
				cc1(rowid_start + m) = tuning_const_group(i);
			}
			rowid_start = rowid_start + miss_group_counts(i);
		}	
		
		mat res( n, 2);
		res.col(0) = partialVec;
		res.col(1) = cc1;
		
		return res; 
	} catch( std::exception& __ex__ ){
		forward_exception_to_r( __ex__ );
	} catch(...){
		::Rf_error( "c++ exception " "(unknown reason)" );
	}
}


// May 17, 2013: Add stopping criteria for solving the scale
double solve_scales( vec maj, vec cc1, double tol, int miter )
{
	int n = maj.n_elem;
	uvec maj_ind = sort_index( maj );
	vec maj_sort(n);
	vec cc1_sort(n);
	vec cc1_sort_cum(n);
	double scc1_half = sum(cc1)/2;
	double scc1_tmp = 0;
	int mid_ind = 0;
	bool found = false;
	for(int i = 0; i < n; i++){
		maj_sort(i) = maj( maj_ind(i) );
		cc1_sort(i) = cc1( maj_ind(i) );
		scc1_tmp += cc1_sort(i);
		cc1_sort_cum(i) = scc1_tmp;
		if( cc1_sort_cum(i) >= scc1_half  && !found){
			mid_ind = i;
			found = true;
		}
	}
	
	double s0 = maj_sort(mid_ind);
	double eps = 1.0;
	int iter= 0;
	while( eps > tol){
		vec zz;
		zz = rho1(maj_sort/s0);
		vec zzcc1 = cc1_sort % zz;
		double s1 = s0 * sum(zzcc1) / scc1_half;
		eps = fabs(s1 - s0)/s0;
		s0 = s1;
		iter++;
		if( iter > miter ) break;
	}
	return s0;
}

/***************************************************/
/*              Vectorized  Rho                    */
/***************************************************/
// Rho 1 = Tukey Bisquare
// Rho 2 = almost linear
vec rho1(vec x) {
	int n = x.n_elem;
	mat xmat(x);
	mat xxtmp = ones<mat>(n, 2);
	xxtmp.col(0) = xmat;
	vec xx = min(xxtmp, 1);
	vec yy = 1 - xx + (pow(xx, 2)/3);
	vec zz = xx % yy;
	zz = 3 * zz;
	return zz; 
}

/***************************************************/
/*                 Rho Prime                       */
/***************************************************/
double rho1p(double x){
	double z;
	if( x < 1 ){
		z = x;
	} else{
		z = 1;
	}
	double y = 3 - 6 * z + 3 * pow(z, 2);
	return y; 
}