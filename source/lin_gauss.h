#ifndef SOURCE_LINGAUSS_H_
#define SOURCE_LINGAUSS_H_

// *** source/lin_gauss.h ***
// Author: Kevin Wolz, date: 06/2018
//
// This program defines the class "LinGauss" describing linear Gaussian models.

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <stdlib.h>
#include <sys/stat.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h> 
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_integration.h>
#include "utils.h"


class LinGauss {
 public:		
	unsigned int datdim_;
	unsigned int pardim_;
	// Ensures that g_matrix_, covar_matrix_, prior_matrix_, prior_param_values_, 
	// data_mean_values_ are initialized and != 0.
	bool initialized_;
	// Ensures that cover_matrix_ is inverted and stored in inv_cover_matrix_.
	bool inverted_;
	// Ensures that this instance of LinGauss is initialized and fisher_matrix_
	// is calculated.
	bool fisherized_;	
	gsl_matrix* g_matrix_; 						// design matrix of the model
	gsl_matrix* covar_matrix_;				// data covariance matrix
	gsl_matrix* inv_covar_matrix_; 		// inverse data covariance matrix
	gsl_matrix* fisher_matrix_; 			// Fisher matrix
	gsl_matrix* prior_matrix_; 			 	// prior precision matrix
	gsl_vector* prior_param_values_; 	// prior fiducial parameter values
	gsl_vector* param_values_; 				// ML estimator values
	gsl_vector* data_mean_values_; 		// mean of the data vector
	
	LinGauss(unsigned int num_data, unsigned int num_params);
	void DisallocateMembers();
	void CalculateFisherMatrix();
	void CalculateMlEstimator(gsl_vector* data);
	double CalculateChisquared(gsl_vector* data, gsl_vector* params);
};
#endif // SOURCE_LINGAUSS_H_
