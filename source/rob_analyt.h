#ifndef SOURCE_ROBANALYT_H_
#define SOURCE_ROBANALYT_H_

// *** source/rob_analyt.h ***
// Author: Kevin Wolz, date: 08/2018
//
// Class RobAnalyt: analytical algorithm to calculate distribution in the WPL.

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
#include <gsl/gsl_randist.h>
#include <gsl/gsl_integration.h>

#include "utils.h"
#include "print_to_file.h"


class RobAnalyt {	
 public:	
 	gsl_matrix* gaussian_quadratic_matrix_;		// quadratic contraction with data
	unsigned int num_steps_;		 							// number of PDF evaluation points
	double lower_limit_;				 							// lower limit of PDF calculation
	double upper_limit_;				 							// upper limit of PDF calculation
	
	RobAnalyt(int num_steps,
	 					double upper_limit,
	 					double lower_limit,
	          unsigned int datdim, 							// number of data points
	          unsigned int pardim, 							// number of models parameters
	          gsl_matrix* g_matrix, 						// model design matrix
	          gsl_matrix* inv_covar_matrix, 		// inverse data covariance matrix
	          gsl_matrix* fisher_matrix_tot, 		// Fisher matrix of total data
						gsl_matrix* fisher_matrix_1, 			// Fisher matrix of first subset
						gsl_matrix* fisher_matrix_2, 			// Fisher matrix of second subset
						gsl_matrix* prior_matrix,
						gsl_vector* prior_param_values, 	// prior fiducial parameter values
						gsl_vector* data_mean_values); 		// mean of data vector

	int CalculateAnalytWeakPrior(std::string& filepath);

	
 private:
	unsigned int datdim_; 							
	unsigned int pardim_; 										
	gsl_matrix* g_matrix_tot_; 				
	gsl_matrix* g_matrix_1_; 			
	gsl_matrix* g_matrix_2_; 			
	gsl_matrix* inv_covar_matrix_tot_; 	
	gsl_matrix* inv_covar_matrix_1_; 		
	gsl_matrix* inv_covar_matrix_2_;
	gsl_matrix* fisher_matrix_tot_;
	gsl_matrix* fisher_matrix_1_;
	gsl_matrix* fisher_matrix_2_;
	gsl_matrix* prior_matrix_;
	gsl_vector* prior_param_values_; 
	gsl_vector* data_mean_values_;
	
	double delta_y_; 											// constant shift of Robustness (1)
	double delta_y_wpl_;									// constant shift of Robustness (WPL, 1)
	double rob_constant_;									// constant shift of Robustness (2)
	gsl_vector* linear_gaussian_vector_; 	// linear contraction with data
	gsl_matrix* covar_matrix_tot_; 				// covariance matrix of total data
};

#endif // SOURCE_ROBANALYT_H_
