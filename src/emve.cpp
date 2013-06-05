#include "emve.h"
using namespace Rcpp;
using namespace arma;
using namespace std;

/***************************************************/
/*                  Some constants                 */
/***************************************************/
#define MAX_RESAMPLE_ITER 50
#define ABS_MIN_RCOND 0.0000000001
//#define MIN_RCOND 0.0001

/***************************************************/
/*               Function prototypes               */
/***************************************************/
double rcondCov(mat A);
SEXP partial_mahalanobis(SEXP X_mu_diff, SEXP X_nonmiss, SEXP Sigma);
SEXP emve_resamp(SEXP X, SEXP X_nonmiss, SEXP nRes, SEXP nSubSize);
SEXP emve_scale_missing( SEXP Sigma, SEXP Miss_group_unique, SEXP Miss_group_counts );

/***************************************************/
/*                Called functions                 */
/***************************************************/
SEXP emve_resamp(SEXP X, SEXP X_nonmiss, SEXP nRes, SEXP nSubSize, SEXP MinRcondition)
{
	try{
		mat x = as<mat>(X);
		mat x_nonmiss = as<mat>(X_nonmiss);
		int nResample = as< int >(nRes);
		int nSubsampleSize = as< int >(nSubSize);
		double minRcondition = as<double>(MinRcondition);
		int n = x.n_rows;
		unsigned int p = x.n_cols;

		// Some initialization before resampling
		int i, j;
		int subsample_ind_i; 								// a subsample index
		mat subsample( nSubsampleSize, p); 					// subsample of size n
		mat subsample_nonmiss( nSubsampleSize, p);
		mat subsample_ind_list( nResample, nSubsampleSize); // a list of index of each good conditioned subsample
		double subsample_rcond = 0;
		
		// Start resampling
		int k_res_iter;
		int step = (int)( MAX_RESAMPLE_ITER / 5 );
		double minRcond = minRcondition;
		double curMinRcond = minRcondition;
		
		
		// Aug 4, 2012
		// Set seed
		srand( 99 );
		
		for( i=0; i < nResample; i++){
			k_res_iter = 0;
			subsample_rcond = 0;
			minRcond = minRcondition;
			while( subsample_rcond < minRcond && k_res_iter < MAX_RESAMPLE_ITER){
				// Obtain a sample
				//subsample.clear();
				
				vec subsample_ind = randu<vec>( nSubsampleSize );

				for(j = 0; j < nSubsampleSize; j++){
					subsample_ind_i = (int)floor(subsample_ind(j) * n); // random int from 0 to n
					if( subsample_ind_i == n ) subsample_ind_i -= 1;
					subsample.row(j) = x.row( subsample_ind_i );
					subsample_nonmiss.row(j) = x_nonmiss.row( subsample_ind_i );
					subsample_ind_list(i, j) = subsample_ind_i + 1; // index in R starts at 1 not 0
				}
				// Check if subsample contain completely missing column
				rowvec b = sum( subsample_nonmiss );
				if( b.min() > 0 ){	
					subsample_rcond = rcondCov( subsample );
					if( k_res_iter % step == 0 && minRcond > ABS_MIN_RCOND){
						minRcond = minRcond / 10;
					}
					k_res_iter++;
					if( minRcond < curMinRcond ){
						curMinRcond = minRcond;
					}
				}
			}
			if( k_res_iter >= MAX_RESAMPLE_ITER )
				throw std::range_error("Fail to obtain an initial estimate of EMVE due to singular subsample matrix.");
		}
		return List::create( Named("subsample_index_list")=subsample_ind_list, Named("curMinRcond")=curMinRcond );
	} catch( std::exception& __ex__ ){
		forward_exception_to_r( __ex__ );
	} catch(...){
		::Rf_error( "c++ exception " "(unknown reason)" );
	}
}


SEXP fast_partial_mahalanobis(SEXP X_mu_diff, SEXP Sigma, SEXP Miss_group_unique, SEXP Miss_group_counts)
{
	try{
		mat x_mu_diff = as<mat>(X_mu_diff);
		mat sigma = as<mat>(Sigma);
		umat miss_group_unique = as<umat>(Miss_group_unique);
		uvec miss_group_counts = as<uvec>(Miss_group_counts);
		
		int n_counts = miss_group_unique.n_rows;
		int n = x_mu_diff.n_rows;
		unsigned int p = x_mu_diff.n_cols;

		vec partialVec(n);
		uvec pp = sum(miss_group_unique, 1);
		int rowid_start = 0;
		for(int i = 0; i < n_counts; i++){
			mat sigma_nonmiss( pp(i), pp(i) );
			mat xi( miss_group_counts(i) , pp(i) );
			int rowid_end = rowid_start + miss_group_counts(i) - 1;
			if( pp(i) < p ){
				int mm = 0;
				for(unsigned int j=0; j<p; j++){
					int nn=mm;
					if(miss_group_unique(i,j) == 1){
						for(unsigned int k=j; k<p; k++){
							if( miss_group_unique(i,k) == 1 ){
								sigma_nonmiss(mm, nn) = sigma(j,k);
								sigma_nonmiss(nn, mm) = sigma(k,j);
								nn++;
							}
						}
						xi.col(mm) = x_mu_diff( span(rowid_start, rowid_end ), j );
						mm++;
					}
				}
			} else{
				sigma_nonmiss = sigma;
				xi = x_mu_diff.rows(rowid_start, rowid_end);
			}

			mat A = ones<mat>( pp(i), pp(i));
			mat diagA = diagmat(A);
			mat sigma_nonmiss_inv = solve( sigma_nonmiss, diagA );

			for(unsigned int m = 0; m < miss_group_counts(i); m++){
				mat xii = xi.row(m);
				partialVec(rowid_start + m) = as_scalar(xii * sigma_nonmiss_inv * trans(xii));
			}
			rowid_start = rowid_start + miss_group_counts(i);
		}	
		return wrap(partialVec); 
	} catch( std::exception& __ex__ ){
		forward_exception_to_r( __ex__ );
	} catch(...){
		::Rf_error( "c++ exception " "(unknown reason)" );
	}
}

// Computes the reciprical condition number of a matrix
// The function makes use of Armadillo solve function to invert a matrix approximately using 
// a fast algorithm
double rcondCov(mat A){
	try{
		int p = A.n_cols;
		mat covA = cov(A);
		vec covAeval = eig_sym(covA);
		double rconditionNum = covAeval(0)/covAeval(p-1);
		if( rconditionNum < 0) rconditionNum = 0.0;
		return rconditionNum;
	} catch( std::exception& __ex__ ){
		forward_exception_to_r( __ex__ );
	} catch(...){
		::Rf_error( "c++ exception " "(unknown reason)" );
	}
}


SEXP emve_scale_missing( SEXP Sigma, SEXP Miss_group_unique, SEXP Miss_group_counts )
{
	try{
		using namespace Rcpp;
		using namespace arma;

		mat sigma = as<mat>(Sigma);
		umat miss_group_unique = as<umat>(Miss_group_unique);
		uvec miss_group_counts = as<uvec>(Miss_group_counts);
		
		int n_counts = miss_group_unique.n_rows;
		unsigned int p = sigma.n_cols;
		
		mat sigma_nonmiss;

		double a = 0;
		double val, sign;
		uvec pp = sum( miss_group_unique, 1);
		for(int i = 0; i < n_counts; i++){
			mat sigma_nonmiss( pp(i), pp(i) );
			if( pp(i) < p ){
				int mm = 0;
				for(unsigned int j=0; j<p; j++){
					int nn=mm;
					if(miss_group_unique(i,j) == 1){
						for(unsigned int k=j; k<p; k++){
							if( miss_group_unique(i,k) == 1 ){
								sigma_nonmiss(mm, nn) = sigma(j,k);
								sigma_nonmiss(nn, mm) = sigma(k,j);
								nn++;
							}
						}
						mm++;
					}
				}
			} else{
				sigma_nonmiss = sigma;
			}

			
			log_det(val, sign, sigma_nonmiss);
			a += val * miss_group_counts(i);
		}	

		double ppp = as_scalar(sum(pp % miss_group_counts));
		a = exp(-1 * a / ppp);

		return List::create( Named("a")=a );
	} catch( std::exception& __ex__ ){
		forward_exception_to_r( __ex__ );
	} catch(...){
		::Rf_error( "c++ exception " "(unknown reason)" );
	}
}

