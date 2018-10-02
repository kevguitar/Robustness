// *** source/rob_analyt.cc ***
// Author: Kevin Wolz, date: 08/2018
//
// Class RobAnalyt: analytical algorithm to calculate distribution in the WPL.

#include "rob_analyt.h"


// Calculates the mean ML estimator vector, the gaussian_quadratic_matrix_ and
// the linear_gaussian_vector_.
RobAnalyt::RobAnalyt(int num_steps, 					// number of PDF evaluation points
	 									 double upper_limit, 			// upper limit of PDF calculation
	 									 double lower_limit, 			// lower limit of PDF calculation
	          				 unsigned int datdim, 		// number of data points
	          				 unsigned int pardim, 		// number of models parameters
	          				 gsl_matrix* g_matrix, 							// model design matrix
	          				 gsl_matrix* inv_covar_matrix, 			
	         					 gsl_matrix* fisher_matrix_tot,	// Fisher matrix of all data
										 gsl_matrix* fisher_matrix_1, 	// Fisher matrix of set 1
										 gsl_matrix* fisher_matrix_2, 	// Fisher matrix of set 2
										 gsl_matrix* prior_matrix,
										 gsl_vector* prior_param_values,// prior fiducial parameters
										 gsl_vector* data_mean_values)     // mean of data vector
    : num_steps_(num_steps), upper_limit_(upper_limit),
      lower_limit_(lower_limit),
      datdim_(datdim), pardim_(pardim), g_matrix_tot_(g_matrix),
      fisher_matrix_tot_(fisher_matrix_tot), fisher_matrix_1_(fisher_matrix_1),
      fisher_matrix_2_(fisher_matrix_2), prior_matrix_(prior_matrix),
      prior_param_values_(prior_param_values),
      data_mean_values_(data_mean_values) {
	if (datdim%2 != 0) {
	  std::cout<<"@ "<<__FUNCTION__<<": 'datdim' is not an even number. Please \
	  						change accordingly. ";
	  abort(); // !!!	IMPORTANT: datdim = dim(total data)	!!!
	}
	if ((g_matrix->size1 != pardim) || (g_matrix->size2 != datdim))	{
	  std::cout<<"@ "<<__FUNCTION__<<": matrix 'g_matrix' does not match \
	  						dimensions. Please allocate correctly. ";
	  abort();
	}
	if ((inv_covar_matrix->size1 != datdim) 
			|| (inv_covar_matrix->size2 != datdim)) {
	  std::cout<<"@ "<<__FUNCTION__<<": matrix 'inv_covar_matrix' does not match \
								dimensions. Please allocate correctly. ";
	  abort();
	}
	if ((fisher_matrix_tot->size1 != pardim) 
			|| (fisher_matrix_tot->size2 != pardim)) {
	  std::cout<<"@ "<<__FUNCTION__<<": matrix 'fisher_matrix_tot' does not \
								match dimensions. Please allocate correctly. ";
	  abort();
	}
	if ((fisher_matrix_1->size1 != pardim) 
			|| (fisher_matrix_1->size2 != pardim)) {
	  std::cout<<"@ "<<__FUNCTION__<<": matrix 'fisher_matrix_1' does not match \
								dimensions. Please allocate correctly. ";
	  abort();
	}
	if ((fisher_matrix_2->size1 != pardim) 
			|| (fisher_matrix_2->size2 != pardim)) {
	  std::cout<<"@ "<<__FUNCTION__<<": matrix 'fisher_matrix_2' does not match \
								dimensions. Please allocate correctly. ";
	  abort();
	}
	if ((prior_matrix->size1 != pardim) || (prior_matrix->size2 != pardim))	{
	  std::cout<<"@ "<<__FUNCTION__<<": matrix 'prior_matrix' does not match \
	  						dimensions. Please allocate correctly. ";
	  abort();
	}
	if (prior_param_values->size != pardim) {
	  std::cout<<"@ "<<__FUNCTION__<<": vector 'prior_param_values' does not \
	  						match dimensions. Please allocate correctly. ";
	  abort();
	}
	if (data_mean_values->size != datdim) {
	  std::cout<<"@ "<<__FUNCTION__<<": vector 'data_mean_values' does not match \
	  						dimensions. Please allocate correctly. ";
	  abort();
	}
	if ((gsl_matrix_isnull(g_matrix) == 1) 
			|| (gsl_matrix_isnull(inv_covar_matrix) == 1) 
			|| (gsl_matrix_isnull(prior_matrix) == 1)) {
		std::cout<<"@ "<<__FUNCTION__<<": the matrices 'g_matrix', \
								'inv_covar_matrix' and 'prior_matrix' are not initialized to \
								values != 0. Please do so before proceeding. ";
		abort(); 
	}
	if ((gsl_matrix_isnull(fisher_matrix_tot) == 1) 
			|| (gsl_matrix_isnull(fisher_matrix_1) == 1)
			|| (gsl_matrix_isnull(fisher_matrix_2) == 1))	{
		std::cout<<"@ "<<__FUNCTION__<<": the matrices 'fisher_matrix_tot', \
								'fisher_matrix_1' and 'fisher_matrix_2' are not initialized to \
								values != 0. Please do so before proceeding. ";
		abort(); 
	}

	// Class members of RobAnalyt 
	linear_gaussian_vector_ = gsl_vector_alloc(datdim);
	gaussian_quadratic_matrix_ = gsl_matrix_alloc(datdim, datdim);
	g_matrix_1_ = gsl_matrix_alloc(pardim, datdim/2);		
	g_matrix_2_ = gsl_matrix_alloc(pardim, datdim/2);	
	inv_covar_matrix_1_ = gsl_matrix_alloc(datdim/2, datdim/2);	
	inv_covar_matrix_2_ = gsl_matrix_alloc(datdim/2, datdim/2);
	covar_matrix_tot_ = gsl_matrix_alloc(datdim, datdim);
	
	//	Intermediate quantities
	gsl_matrix* inv_fisher_matrix_1 = gsl_matrix_alloc(pardim, pardim); // IL1
	gsl_matrix* inv_fisher_matrix_2 = gsl_matrix_alloc(pardim, pardim); // IL2
	gsl_matrix* inv_fisher_matrix_tot = gsl_matrix_alloc(pardim, pardim); // IL
	gsl_matrix* total_fisher_1 = gsl_matrix_calloc(pardim, pardim); // F1
	gsl_matrix* total_fisher_2 = gsl_matrix_calloc(pardim, pardim); // F2
	gsl_matrix* total_fisher_tot = gsl_matrix_calloc(pardim, pardim); // F
	gsl_matrix* inv_total_fisher_1 = gsl_matrix_alloc(pardim, pardim); // IF1
	gsl_matrix* inv_total_fisher_2 = gsl_matrix_alloc(pardim, pardim); // IF2
	gsl_matrix* inv_total_fisher_tot = gsl_matrix_alloc(pardim, pardim); // IF
	gsl_matrix* t_prime_matrix = gsl_matrix_alloc(datdim/2, datdim/2); // T1

	gsl_matrix* inv_gaussian_quadratic_matrix 
			= gsl_matrix_alloc(datdim, datdim); // IT
	gsl_matrix* gaussian_quadratic_matrix_trans 
			= gsl_matrix_alloc(datdim, datdim); // Tt
	gsl_matrix* ILg1 = gsl_matrix_alloc(pardim,datdim/2);
	gsl_matrix* LgSig1 = gsl_matrix_alloc(datdim/2, pardim);
	gsl_vector* parmean1 = gsl_vector_alloc(pardim);
	
	for(int i=0; i<pardim; i++) {
	  for(int j=0; j<datdim; j++)	{
		  if (j<datdim/2)	{
		    gsl_matrix_set(g_matrix_1_, i, j, gsl_matrix_get(g_matrix, i, j));
		  } else {
		    gsl_matrix_set(g_matrix_2_, i, (j - datdim/2), 
		    	             gsl_matrix_get(g_matrix ,i ,j));
		  }  	 
		}
	}
	
	// Calculation of Fisher and T-matrices.
	for(int i=0; i<pardim; i++) {
		for(int j=0; j<datdim; j++)	{
			if (j<datdim/2 ) {
				gsl_matrix_set(g_matrix_1_, i, j, gsl_matrix_get(g_matrix, i, j));
			} else {
			  gsl_matrix_set(g_matrix_2_, i, (j - datdim/2),
			                 gsl_matrix_get(g_matrix, i, j));	
			}                  
		}
	}
	for( int i=0; i<datdim; i++ )		{
		for( int j=0; j<datdim; j++ )	{
			if ((j<datdim/2) && (i<datdim/2))	{
			  gsl_matrix_set(inv_covar_matrix_1_, i, j,
			                 gsl_matrix_get(inv_covar_matrix, i, j));
			} else if ((j>=datdim/2) && (i>=datdim/2)) {
			  gsl_matrix_set(inv_covar_matrix_2_, (i - datdim/2), (j - datdim/2), 
			                 gsl_matrix_get(inv_covar_matrix, i, j));
			}	 
		}
	}		
	
	gsl_matrix_memcpy(total_fisher_1, fisher_matrix_1);
	gsl_matrix_add(total_fisher_1, prior_matrix);				 		// total_fisher_1
	gsl_matrix_memcpy(total_fisher_2, fisher_matrix_2);
	gsl_matrix_add(total_fisher_2, prior_matrix);						// total_fisher_2
	gsl_matrix_memcpy(total_fisher_tot, fisher_matrix_tot);
	gsl_matrix_add(total_fisher_tot, prior_matrix);				  // total_fisher_tot
	
	if (DeterminantGslMatrix(total_fisher_1, pardim) == 0.0) {
		std::cout<<__FILE__<<", l. "<<__LINE__<<": Matrix is singular: "<<std::endl;
		PrintGslMatrix(total_fisher_1,"total_fisher_1");
		abort();
	}
	if (DeterminantGslMatrix(total_fisher_2, pardim) == 0.0) {
		std::cout<<__FILE__<<", l. "<<__LINE__<<": Matrix is singular: "<<std::endl;
		PrintGslMatrix(total_fisher_2,"total_fisher_2");
		abort();
	}
	if (DeterminantGslMatrix(total_fisher_tot, pardim) == 0.0) {
		std::cout<<__FILE__<<", l. "<<__LINE__<<": Matrix is singular: "<<std::endl;
		PrintGslMatrix(total_fisher_tot,"total_fisher_tot");
		abort();
  }
	if (DeterminantGslMatrix(fisher_matrix_1, pardim) == 0.0) {
		std::cout<<__FILE__<<", l. "<<__LINE__<<": Matrix is singular: "<<std::endl;
		PrintGslMatrix(fisher_matrix_1,"fisher_matrix_1");
		abort();
  }
	if (DeterminantGslMatrix(fisher_matrix_2, pardim) == 0.0) {
		std::cout<<__FILE__<<", l. "<<__LINE__<<": Matrix is singular: "<<std::endl;
		PrintGslMatrix(fisher_matrix_2,"fisher_matrix_2");
		abort();
	}
	if (DeterminantGslMatrix(fisher_matrix_tot, pardim) == 0.0) {
		std::cout<<__FILE__<<", l. "<<__LINE__<<": Matrix is singular: "<<std::endl;
		PrintGslMatrix(fisher_matrix_tot,"fisher_matrix_tot");
		abort();
	}
	InvertGslMatrix(total_fisher_1, pardim, 
						inv_total_fisher_1);			// inv_total_fisher_1
	InvertGslMatrix(total_fisher_2, pardim,
					  inv_total_fisher_2);			// inv_total_fisher_2
	InvertGslMatrix(total_fisher_tot, pardim,
						inv_total_fisher_tot);		// inv_total_fisher_tot
	InvertGslMatrix(fisher_matrix_1, pardim,
						inv_fisher_matrix_1);			// inv_fisher_matrix_1
	InvertGslMatrix(fisher_matrix_2, pardim,
				 		inv_fisher_matrix_2);			// inv_fisher_matrix_2
	InvertGslMatrix(fisher_matrix_tot, pardim,
						inv_fisher_matrix_tot);		// inv_fisher_matrix_tot
	
	
	gsl_vector* Sigd1 = gsl_vector_alloc(datdim/2);
	gsl_vector* gSigd1 = gsl_vector_alloc(pardim);
	gsl_vector* data_mean_vector_1 = gsl_vector_alloc(datdim/2);

	for (int i=0; i<datdim/2; i++) {
	  double cur = gsl_vector_get(data_mean_values,i);
	  gsl_vector_set(data_mean_vector_1, i, cur);					// data_mean_vector_1
	}
	gsl_blas_dgemv(CblasNoTrans, 1.0, inv_covar_matrix_1_, data_mean_vector_1,
								 0.0, Sigd1);
	gsl_blas_dgemv(CblasNoTrans, 1.0, g_matrix_1_, Sigd1, 0.0, gSigd1);
	gsl_blas_dgemv(CblasNoTrans, 1.0, inv_fisher_matrix_1, gSigd1, 0.0, 
								 parmean1);															// parmean1
	
	gsl_matrix* gSiginv1 = gsl_matrix_alloc(pardim, datdim/2);
	gsl_matrix* ILgSiginv1 = gsl_matrix_alloc(pardim, datdim/2);
	gsl_matrix* IFP1 = gsl_matrix_alloc(pardim, pardim);
	gsl_matrix* IFPhalf = gsl_matrix_alloc(pardim, pardim);
	gsl_matrix* Ihalf = gsl_matrix_alloc(pardim, pardim);
	gsl_matrix* halfbracket = gsl_matrix_alloc(pardim, pardim);
	gsl_matrix* fullbracket = gsl_matrix_alloc(pardim, pardim);
	gsl_matrix* firstpart = gsl_matrix_alloc(datdim/2, pardim);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, g_matrix_1_,
								 inv_covar_matrix_1_, 0.0, gSiginv1);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, inv_fisher_matrix_1,
								 gSiginv1, 0.0, ILgSiginv1);	
	gsl_matrix_set_identity(Ihalf);
	gsl_matrix_scale(Ihalf, -0.5);	
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 0.5, inv_total_fisher_tot,	
								 prior_matrix, 0.0, IFPhalf);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, inv_total_fisher_1,
								 prior_matrix, 0.0, IFP1);
	gsl_matrix_memcpy(halfbracket, Ihalf);
	gsl_matrix_sub(halfbracket, IFPhalf);
	gsl_matrix_memcpy(fullbracket, halfbracket);
	gsl_matrix_add(fullbracket, IFP1);
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, gSiginv1, fullbracket, 0.0,
	               firstpart);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, firstpart, ILgSiginv1, 0.0,
	               t_prime_matrix);
	BlockDiagFromTwoGslMatrices(t_prime_matrix, t_prime_matrix, 
															gaussian_quadratic_matrix_, datdim, datdim);
															// gaussian_quadratic_matrix_
	
	gsl_matrix_free(gSiginv1);
	gsl_matrix_free(ILgSiginv1);
	gsl_matrix_free(IFP1);
	gsl_matrix_free(IFPhalf);
	gsl_matrix_free(Ihalf);
	gsl_matrix_free(halfbracket);
	gsl_matrix_free(fullbracket);
	gsl_matrix_free(firstpart);
	
	if (DeterminantGslMatrix(inv_covar_matrix, datdim) == 0.0) {
		std::cout<<__FILE__<<", l. "<<__LINE__<<": Matrix is singular: "<<std::endl;
		PrintGslMatrix(inv_covar_matrix, "inv_covar_matrix");
		abort();
	}	
	InvertGslMatrix(inv_covar_matrix, datdim, covar_matrix_tot_);	 
									// covar_matrix_tot_

	// Calculation of delta_y_ and delta_y_wpl_.
	delta_y_ = log(DeterminantGslMatrix(total_fisher_1, pardim)
								 * DeterminantGslMatrix(total_fisher_2, pardim)
								 / DeterminantGslMatrix(total_fisher_tot, pardim)
								 / DeterminantGslMatrix(prior_matrix, pardim));  
								 // delta_y_
	delta_y_wpl_ = log(DeterminantGslMatrix(fisher_matrix_1, pardim)
								 		 * DeterminantGslMatrix(fisher_matrix_2, pardim)
								 		 / DeterminantGslMatrix(fisher_matrix_tot, pardim)
								 	 	 / DeterminantGslMatrix(prior_matrix, pardim));  
								 	 	 // delta_y_wpl

	
	gsl_matrix* innerbracket = gsl_matrix_alloc(pardim, pardim);
	gsl_matrix* Linnerbracket = gsl_matrix_alloc(pardim, pardim);
	gsl_matrix* outterbracket = gsl_matrix_alloc(pardim, pardim);
	gsl_matrix_memcpy(innerbracket, inv_total_fisher_tot);
	gsl_matrix_sub(innerbracket, inv_total_fisher_1);	
	gsl_matrix* innerP = gsl_matrix_alloc(pardim, pardim);
	gsl_matrix* ginnerP = gsl_matrix_alloc(datdim/2, pardim);
	gsl_matrix* ICginnerP = gsl_matrix_alloc(datdim/2, pardim);
	gsl_vector* uprime = gsl_vector_alloc(datdim/2);

	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, innerbracket,
								 prior_matrix, 0.0, innerP);
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, g_matrix_1_,
								 innerP, 0.0, ginnerP);
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 2.0, inv_covar_matrix_1_, ginnerP,
								 0.0, ICginnerP);
	gsl_blas_dgemv(CblasNoTrans, 1.0, ICginnerP, prior_param_values, 0.0, uprime);
	AppendTwoGslVectors(uprime, uprime,
						 linear_gaussian_vector_);		// linear_gaussian_vector_
	
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 2.0, fisher_matrix_1,
								 innerbracket, 0.0, Linnerbracket);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Linnerbracket, prior_matrix,
								 0.0, outterbracket);

	// First part of the calculation of "rob_constant_"               
	rob_constant_ = - 0.5 * GslVec1MatVec2(prior_param_values, outterbracket,
	                               			 prior_param_values, pardim) 
									+ 0.5 * delta_y_;				
	
	gsl_matrix_free (innerbracket);
	gsl_matrix_free (Linnerbracket);
	gsl_matrix_free (outterbracket);
	gsl_matrix_free (innerP);
	gsl_matrix_free (ginnerP);
	gsl_matrix_free (ICginnerP);
	gsl_vector_free (uprime);
	
  // Symmetrizes gaussian_quadratic_matrix_.
	gsl_matrix_memcpy(gaussian_quadratic_matrix_trans,
									  gaussian_quadratic_matrix_);
	gsl_matrix_transpose(gaussian_quadratic_matrix_trans);
	gsl_matrix_add(gaussian_quadratic_matrix_, gaussian_quadratic_matrix_trans);				
	gsl_matrix_scale(gaussian_quadratic_matrix_,
									 0.5);	      // gaussian_quadratic_matrix_
	
	// Deallocates intermediate quantities.
	gsl_matrix_free (inv_fisher_matrix_1);
	gsl_matrix_free (inv_fisher_matrix_2);
	gsl_matrix_free (inv_fisher_matrix_tot);
	gsl_matrix_free (total_fisher_1);
	gsl_matrix_free (total_fisher_2);
	gsl_matrix_free (total_fisher_tot);
	gsl_matrix_free (inv_total_fisher_1);
	gsl_matrix_free (inv_total_fisher_2);
	gsl_matrix_free (inv_total_fisher_tot);
	gsl_matrix_free (t_prime_matrix);
	gsl_matrix_free (gaussian_quadratic_matrix_trans);
	gsl_vector_free (Sigd1);
	gsl_vector_free (gSigd1);
	gsl_vector_free (data_mean_vector_1);
	gsl_vector_free (parmean1);
}

							
// Calculates the analytical PDF and CDF of the Robustness in the weak prior
// limit (WPL).
int RobAnalyt::CalculateAnalytWeakPrior(std::string& filepath) {
	// Opens the data file to write to.
	std::string datapath("data/");
	datapath.append(Timestamp());
	datapath.append("_21.rana");
  filepath.append(datapath);
  std::ofstream datafile(datapath, std::ofstream::out | std::ofstream::app);

  // Calculates the data and writes them to the file.
	for (int i=0; i<num_steps_; i++) {
		int p = 1e6*i / num_steps_;
		if ((p%10000) == 0)	{
		  std::cout<<"Evaluating analytical distribution in weak prior limit... "
		           <<p/10000<<" percent completed."<<std::endl;
		}
		// Note: xi is measured in units of log-robustness R.
		double xi = lower_limit_ + i*(upper_limit_ - lower_limit_)/num_steps_;
		double qi = - 2*xi + delta_y_wpl_;
		double pdf_cur = 2*gsl_ran_chisq_pdf(qi, pardim_);

		WriteData(datafile, xi, pdf_cur);
	}
	return 1;
}

