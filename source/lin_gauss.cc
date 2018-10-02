// *** source/lin_gauss.cc ***
// Author: Kevin Wolz, date: 06/2018
//
// This program defines the class "LinGauss" describing linear Gaussian models.

#include "lin_gauss.h"


LinGauss::LinGauss(unsigned int datdim, unsigned int pardim)
	  : datdim_(datdim), pardim_(pardim) {
	inv_covar_matrix_ = gsl_matrix_calloc(datdim, datdim); 
	covar_matrix_ = gsl_matrix_calloc(datdim, datdim);
	g_matrix_ = gsl_matrix_alloc(pardim, datdim);
	fisher_matrix_ = gsl_matrix_alloc(pardim, pardim);
	prior_matrix_ = gsl_matrix_alloc(pardim, pardim);
	prior_param_values_ = gsl_vector_alloc(pardim);
	param_values_ = gsl_vector_alloc(pardim);
	data_mean_values_ = gsl_vector_calloc(datdim);
	initialized_ = false;
	inverted_ = false;
	fisherized_ = false;
}


// A manual "destructor" that does not destruct the object but frees space.
void LinGauss::DisallocateMembers() {
	gsl_matrix_free(covar_matrix_);
	gsl_matrix_free(inv_covar_matrix_);
	gsl_matrix_free(g_matrix_);
	gsl_matrix_free(fisher_matrix_);
	gsl_matrix_free(prior_matrix_);
	gsl_vector_free(param_values_);
	gsl_vector_free(data_mean_values_);
	gsl_vector_free(prior_param_values_);
}


void LinGauss::CalculateFisherMatrix() {
	int d = datdim_;
	int p = pardim_;

	if (initialized_ == false) {
	  std::cout<<"@ "<<__FUNCTION__<<": current instance of LinGauss is not \
	  						initialized_. Please initialize to non-zero value. "<<std::endl;
	  abort();
	}
	gsl_matrix* intermed_matrix = gsl_matrix_alloc(d, p);
	if (inverted_ == false) {
		InvertGslMatrix(covar_matrix_, d, inv_covar_matrix_);
		inverted_ = true;
	}	

	// Multiplies inv_covar_matrix_ * g_matrix_^t.
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, inv_covar_matrix_,
	               g_matrix_, 0.0, intermed_matrix);
	// Multiplies g * (inv_covar_matrix_ * g^t).
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, g_matrix_, intermed_matrix,
	               0.0, fisher_matrix_);
	fisherized_ = true;
	gsl_matrix_free(intermed_matrix);
}


void LinGauss::CalculateMlEstimator(gsl_vector* data_vector) {
	int d = datdim_;
	int p = pardim_;

	if (data_vector->size != d) {
	  std::cout<<"@ "<<__FUNCTION__<<": vector 'data_vector' does not match \
	         dimensions. Please allocate correctly. "<<std::endl;
	  abort();
	}
	if (initialized_ == false)	{
	  std::cout<<"@ "<<__FUNCTION__<<": current instance of LinGauss is not \
	         initialized_. Please initialize to non-zero value. "<<std::endl;
	  abort();
	}
	if (inverted_ == false) {
		InvertGslMatrix(covar_matrix_, d, inv_covar_matrix_);
		inverted_ = true;
	}	
	
	gsl_matrix* fisherin = gsl_matrix_alloc(p, p);
	gsl_vector* D = gsl_vector_alloc(p);
	gsl_matrix* intermed = gsl_matrix_alloc(d, p);
	
	// Multiplies inv_covar_matrix_ * g_matrix_^t.
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, inv_covar_matrix_, g_matrix_,
	               0.0, intermed);
	// Multiplies g_matrix_ * (inv_covar_matrix_ * g_matrix_^t).
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, g_matrix_, intermed, 0.0,
	               fisher_matrix_);
	fisherized_ = true;
	// Multiplies (inv_covar_matrix_ * g_matrix_^t)^t * data_vector.
	gsl_blas_dgemv (CblasTrans, 1.0, intermed, data_vector, 0.0, D);
	InvertGslMatrix(fisher_matrix_, p, fisherin);	
	// Multiplies fisherin * D.
	gsl_blas_dgemv (CblasNoTrans, 1.0, fisherin, D, 0.0, param_values_);

	gsl_matrix_free(fisherin);
	gsl_vector_free(D);
	gsl_matrix_free(intermed);
}


// Calculates the chi-squared value of a given dataset with given parameters for
// the linear Gaussian model.
double LinGauss::CalculateChisquared(gsl_vector* data_vector,
                                     gsl_vector* param_vector) {
	int d = datdim_;
	int p = pardim_;
	
	if (data_vector->size != d) {
	  std::cout<<"@ "<<__FUNCTION__<<": vector 'data_vector' does not match \
	         dimensions. Please allocate correctly. "<<std::endl;
	  abort();
	}
	if (param_vector->size != p) {
	  std::cout<<"@ "<<__FUNCTION__<<": vector 'param_vector' does not match \
	         dimensions. Please allocate correctly. "<<std::endl;
	  abort();
	}
	if (initialized_ == false)	{
	  std::cout<<"@ "<<__FUNCTION__<<": current instance of LinGauss is not \
	         initialized. Please initialize to non-zero value. "<<std::endl;
	  abort();
	}
	if (inverted_ == false) {
		InvertGslMatrix(covar_matrix_, d, inv_covar_matrix_);
		inverted_ = true;
	}
	
	gsl_vector* hilf_vector = gsl_vector_alloc(d);
	// param_vector * g_matrix
	gsl_blas_dgemv(CblasTrans, 1.0, g_matrix_, param_vector, 0.0, hilf_vector);

	// g_matrix_ * params_vector - data_vector
	gsl_vector_sub(hilf_vector, data_vector); 

	double chisquared = GslVec1MatVec2(hilf_vector, inv_covar_matrix_, 
	                                   hilf_vector, d);
	gsl_vector_free(hilf_vector);
	
	return chisquared;
}
