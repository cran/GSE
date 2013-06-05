#ifndef _EMVE_H
#define _EMVE_H
#include <RcppArmadillo.h>
RcppExport SEXP emve_resamp(SEXP X, SEXP X_nonmiss, SEXP nRes, SEXP nSubSize, SEXP MinRcondition);
RcppExport SEXP emve_scale_missing( SEXP Sigma, SEXP Miss_group_unique, SEXP Miss_group_counts );
RcppExport SEXP fast_partial_mahalanobis(SEXP X_mu_diff, SEXP Sigma, SEXP Miss_group_unique, SEXP Miss_group_counts);
#endif
