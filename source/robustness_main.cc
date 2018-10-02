// *** source/robustness_main.cc ***
// Author: Kevin Wolz, date: 08/2018
//
// Numerical sampling of the Robustness distribution (GP and WPL).
// int main() manages the choice of the operation mode (analytical WPL and 
// numerical WPL and GP algorithms), parsing of inital configurations, logging
// of data and metadata, and plotting of the results.

#include "robustness_main.h"


// Accepts Fisher determinants and chisquared values and calculates the
// robustness in the weak prior limit.
double RobustnessWeakPriorLimit(double chi_squared_1,    // chi1
                                double chi_squared_2,    // chi2
                                double chi_squared_tot,  // chitot
                                double det_fisher_1,     // det1
                                double det_fisher_2,     // det2
                                double det_fisher_tot,   // dettot
                                double det_prior) {      // detP
	double first = -0.5 * log(det_prior);
	double second = 0.5 * log(det_fisher_1 * det_fisher_2 / det_fisher_tot);
	double third = -0.5 * (chi_squared_tot - chi_squared_1 - chi_squared_2); 
	
	return first + second + third;
}


// Calculates the robustness for general Gaussian linear models.
// Accepts Fisher matrices, chisquared values and ML parameter estimators coming
// from the three data subsets, and the prior matrix and the prior parameters.
double RobustnessGeneralPrior(double chi_squared_1,           // chi1
                              double chi_squared_2,           // chi2
                              double chi_squared_tot,         // chitot
                              gsl_matrix* fisher_matrix_1,  	// fish1
                              gsl_matrix* fisher_matrix_2,  	// fish2
                              gsl_matrix* prior_matrix,       // P
                              gsl_vector* param_values_1,     // theta1
                              gsl_vector* param_values_2,     // theta2
                              gsl_vector* param_values_tot,   // thetatot
                              gsl_vector* prior_param_values, // thetaP
                              unsigned int pardim) {          // dim
	if ((fisher_matrix_1->size1 != pardim) 
			|| (fisher_matrix_1->size2 != pardim)) {
	  std::cout<<"@ "<<__FUNCTION__<<": fisher_matrix_1 does not have the \
								right dimensions. Please provide a suitable gsl_matrix. ";
	  abort();	
	}
	if ((fisher_matrix_2->size1 != pardim) 
		  || (fisher_matrix_2->size2 != pardim)) {
	  std::cout<<"@ "<<__FUNCTION__<<": fisher_matrix_2 does not have the \
								right dimensions. Please provide a suitable gsl_matrix. ";
	  abort();
	}
	if ((prior_matrix->size1 != pardim) || (prior_matrix->size2 != pardim)) {
	  std::cout<<"@ "<<__FUNCTION__<<": prior_matrix does not have the right \
	  						dimensions. Please provide a suitable gsl_matrix. ";
	  abort();
	}
	if (param_values_1->size != pardim)	{
	  std::cout<<"@ "<<__FUNCTION__<<": param_values_1 does not have the \
	  						right dimensions. Please provide a suitable gsl_vector. ";
	  abort();
	}
	if (param_values_2->size != pardim)	{
	  std::cout<<"@ "<<__FUNCTION__<<": param_values_2 does not have the \
	  						right dimensions. Please provide a suitable gsl_vector. ";
	  abort();
	}
	if (param_values_tot->size != pardim)	{
	  std::cout<<"@ "<<__FUNCTION__<<": param_values_tot does not have the \
	  						right dimensions. Please provide a suitable gsl_vector. ";
	  abort();
	  }
	if (prior_param_values->size != pardim)	{ 
	  std::cout<<"@ "<<__FUNCTION__<<": prior_param_values does not have the \
	   						right dimensions. Please provide a suitable gsl_vector. ";
	  abort();
	}
	
	gsl_matrix* L1 = gsl_matrix_alloc(pardim, pardim);
	gsl_matrix* L2 = gsl_matrix_alloc(pardim, pardim);
	gsl_matrix* Ltot = gsl_matrix_alloc(pardim, pardim);
	gsl_matrix* F1 = gsl_matrix_alloc(pardim, pardim);
	gsl_matrix* F2 = gsl_matrix_alloc(pardim, pardim);
	gsl_matrix* Ftot = gsl_matrix_alloc(pardim, pardim);
	gsl_matrix* Fin1 = gsl_matrix_alloc(pardim, pardim);
	gsl_matrix* Fin2 = gsl_matrix_alloc(pardim, pardim);
	gsl_matrix* Fintot = gsl_matrix_alloc(pardim, pardim);
	gsl_matrix* LFP1 = gsl_matrix_alloc(pardim, pardim);
	gsl_matrix* LFP2 = gsl_matrix_alloc(pardim, pardim);
	gsl_matrix* LFPtot = gsl_matrix_alloc(pardim, pardim);
	gsl_vector* diff1 = gsl_vector_alloc(pardim);
	gsl_vector* diff2 = gsl_vector_alloc(pardim);
	gsl_vector* difftot = gsl_vector_alloc(pardim);
	
	gsl_matrix_memcpy(L1, fisher_matrix_1);
	gsl_matrix_memcpy(L2, fisher_matrix_2);
	gsl_matrix_memcpy(F1, L1);
	gsl_matrix_add(F1, prior_matrix);
	gsl_matrix_memcpy(F2, L2);
	gsl_matrix_add(F2, prior_matrix);
	gsl_matrix_memcpy(Ltot, L1);
	gsl_matrix_add(Ltot, L2);
	gsl_matrix_memcpy(Ftot, Ltot);
	gsl_matrix_add(Ftot, prior_matrix);
	if (DeterminantGslMatrix(F1, pardim) == 0.0) {
		std::cout<<__FILE__<<", l. "<<__LINE__<<": Matrix is singular: "
		         <<std::endl;
		PrintGslMatrix(F1,"F1");
		abort();
	}
	InvertGslMatrix(F1, pardim, Fin1);
	if (DeterminantGslMatrix(F2, pardim) == 0.0) {
		std::cout<<__FILE__<<", l. "<<__LINE__<<": Matrix is singular: "<<std::endl;
		PrintGslMatrix(F2,"F2");
		abort();
	}
	InvertGslMatrix(F2, pardim, Fin2);
	if (DeterminantGslMatrix(Ftot, pardim) == 0.0) {
	  std::cout<<__FILE__<<", l. "<<__LINE__<<": Matrix is singular: "<<std::endl;
	  PrintGslMatrix(Ftot,"Ftot");
	  abort();
	}
	InvertGslMatrix(Ftot, pardim, Fintot);
	MultiplyThreeGslMatrices(L1, Fin1, prior_matrix, pardim, LFP1);
	MultiplyThreeGslMatrices(L2, Fin2, prior_matrix, pardim, LFP2);
	MultiplyThreeGslMatrices(Ltot, Fintot, prior_matrix, pardim, LFPtot);
	gsl_vector_memcpy(diff1, param_values_1);
	gsl_vector_add(diff1, prior_param_values);
	gsl_vector_memcpy(diff2, param_values_2);
	gsl_vector_add(diff2, prior_param_values);
	gsl_vector_memcpy(difftot, param_values_tot);
	gsl_vector_add(difftot, prior_param_values);
	
	double det_prior = DeterminantGslMatrix(prior_matrix, pardim);
	double detF1 = DeterminantGslMatrix(F1, pardim);
	double detF2 = DeterminantGslMatrix(F2, pardim);
	double detFtot = DeterminantGslMatrix(Ftot, pardim);
	double chiptot = GslVec1MatVec2(difftot, LFPtot, difftot, pardim);
	double chip1 = GslVec1MatVec2(diff1, LFP1, diff1, pardim);
	double chip2 = GslVec1MatVec2(diff2, LFP2, diff2, pardim);
	
	double first = 0.5 * log(detF1 * detF2 / detFtot / det_prior) ; 
	double second = -0.5 * (chi_squared_tot - chi_squared_1 - chi_squared_2); 
	double third = -0.5 * (chiptot - chip1 - chip2);
	
	return first + second + third;
}


// Performs a repeated random sampling of data based on the models provided.
// Calculates the robustness for each sample and stores the result of all
// repetitions inside a vector. 
int NumericalSamplingOfRobustness(unsigned int repetitions,
                                  unsigned int datdim,
                                  bool is_weak_prior_limit,
                                  LinGauss model1,
                                  LinGauss model2,
                                  LinGauss model_tot,            // modelt
                                  std::string& filepath) { 
	if (model1.fisherized_ == false)	{
	  std::cout<<"@ "<<__FUNCTION__<<": instance 'model1' of LinGauss is not \
	              fisherized. Please do so before proceeding. ";
	  abort();
	}
	if (model2.fisherized_ == false)	{
	  std::cout<<"@ "<<__FUNCTION__<<": instance 'model2' of LinGauss is not \
	              fisherized. Please do so before proceeding. ";
	  abort();
	}
	if (model_tot.fisherized_ == false)	{
	  std::cout<<"@ "<<__FUNCTION__<<": instance 'model_tot' of LinGauss is not \
	              fisherized. Please do so before proceeding. ";
	  abort();
	}
	if (((model1.covar_matrix_)->size1 != datdim/2) || 
	    ((model1.covar_matrix_)->size2 != datdim/2)) {
	  std::cout<<"@ "<<__FUNCTION__<<": model1.covar_matrix_ has the wrong \
								dimensions. Please allocate correctly. ";
	  abort();	
	}
	if (((model2.covar_matrix_)->size1 != datdim/2) ||
	    ((model2.covar_matrix_)->size2 != datdim/2)) {
	  std::cout<<"@ "<<__FUNCTION__<<": model2.covar_matrix_ has the wrong \
								dimensions. Please allocate correctly. ";
	  abort();
	}
	if (((model_tot.covar_matrix_)->size1 != datdim) 
			|| ((model_tot.covar_matrix_)->size2 != datdim)) {
	  std::cout<<"@ "<<__FUNCTION__<<": model_tot.covar_matrix_ has the wrong \
	  						dimensions. Please allocate correctly. ";
	  abort();	
	}
	
	// Assumption: dimension of dat1 = dimension of dat2.
	gsl_vector* dat1 = gsl_vector_alloc(datdim/2);
	gsl_vector* dat2 = gsl_vector_alloc(datdim/2);
	gsl_vector* datt = gsl_vector_alloc(datdim);
	gsl_matrix* L1 = gsl_matrix_alloc(datdim/2,datdim/2);
	gsl_matrix* L2 = gsl_matrix_alloc(datdim/2,datdim/2);
	gsl_matrix_memcpy(L1, model1.covar_matrix_);
	gsl_matrix_memcpy(L2, model2.covar_matrix_);
	gsl_linalg_cholesky_decomp1(L1);
	gsl_linalg_cholesky_decomp1(L2);
	
	const gsl_rng_type * T;
  gsl_rng * r;
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  // Opens the data file to write to.
	std::string datapath("./data/");
	datapath.append(Timestamp());
  if (is_weak_prior_limit) datapath.append("_20.rnum");
  else datapath.append("_30.rnum");
  filepath.append(datapath);
  std::ofstream datafile(datapath, std::ofstream::out | std::ofstream::app);

	// Samples the data for calculating the RobustnessGeneralPrior...	
	for (int i=0; i<repetitions; i++)	{
		int p = 1e6*i/repetitions;
		if ((p%10000) == 0)	{
		  std::cout<<"Sampling numerical distribution... "<<p/10000<<" percent \
		              completed."<<std::endl;
		}
		
		// Performs random sampling of data vectors from gaussian distribution.
		gsl_ran_multivariate_gaussian(r, model1.data_mean_values_, L1, dat1);
		gsl_ran_multivariate_gaussian(r, model2.data_mean_values_, L2, dat2);
		AppendTwoGslVectors(dat1, dat2, datt);
		
		// Calculate the best-fit parameters for each dataset.
		model1.CalculateMlEstimator(dat1);
		model2.CalculateMlEstimator(dat2);
		model_tot.CalculateMlEstimator(datt);

		// Calculates the chi-squared of each dataset, using the respective ML 
		// estimators.
		double chi_squared_1 = model1.CalculateChisquared(dat1,
																											model1.param_values_);
		double chi_squared_2 = model2.CalculateChisquared(dat2,
																										  model2.param_values_);
		double chit = model_tot.CalculateChisquared(datt,
			                                          model_tot.param_values_);
		double rob;
		

		if (is_weak_prior_limit) {
			// Calculates and writes the numerical_output of the robustness.
			rob = RobustnessWeakPriorLimit(chi_squared_1, chi_squared_2, chit, 
			                            DeterminantGslMatrix(model1.fisher_matrix_,
			                                        model1.pardim_),
										              DeterminantGslMatrix(model2.fisher_matrix_,
										                          model2.pardim_),
										              DeterminantGslMatrix(model_tot.fisher_matrix_,
										                          model_tot.pardim_),
										              DeterminantGslMatrix(model_tot.prior_matrix_,
										                          model_tot.pardim_));
		} else {
			rob = RobustnessGeneralPrior(chi_squared_1, chi_squared_2, chit,
			                             model1.fisher_matrix_,
			                             model2.fisher_matrix_,
			                             model_tot.prior_matrix_,
			                             model1.param_values_,
			                             model2.param_values_,
			                             model_tot.param_values_,
			                             model_tot.prior_param_values_,
			                             model_tot.pardim_);				
		}
    // For numerical sampling we don't need the analytical PDF.
		double pdf = 0;

		WriteData(datafile, rob, pdf);
	}	
	std::cout<<std::endl;
	datafile.close();
		
	gsl_vector_free(dat1);
	gsl_vector_free(dat2);
	gsl_vector_free(datt);
	gsl_matrix_free(L1);
	gsl_matrix_free(L2);
	gsl_rng_free(r);

  return 1;	
}
				

/*
________________________________________________________________________________

MAIN
________________________________________________________________________________

*/
int main (int argc, char* argv[])	{
	// Parsing of the operation_mode
  if (argc != 2) {
    std::cerr<<"Usage: "<<argv[0]<<" operation_mode"<<std::endl;
    std::cerr<<"Allowed arguments are:"<<std::endl;
    std::cerr<<"  20 - numerical sampling of robustness"
             <<" for multiple parameters (weak prior limit)"<<std::endl;
    std::cerr<<"  21 - analytical distribution of robustness"
             <<" for multiple parameters (weak prior limit)"<<std::endl;
    std::cerr<<"  30 - numerical sampling of robustness"
             <<" for multiple parameters";

    return 1;
  }

  // Initial operational parameters.
	unsigned datdim;		
	unsigned pardim; 
	unsigned num_sampled_data;	 
	double lower_limit_analyt_distribution;			
	double upper_limit_analyt_distribution;
	unsigned num_steps_analyt_distribution;						 								 									
																																										
	int operation_mode = atoi(argv[1]);			 
	std::cout<<"You have chosen operation mode "<<operation_mode<<"."<<std::endl;																			
	
  // Parses and stores information about the dimensions.
  std::string filepath("init/");
  if (operation_mode == 20) {
    filepath.append("wpl_num.init");
  } else if (operation_mode == 21) {
    filepath.append("wpl_ana.init");
  } else if (operation_mode == 30) {
    filepath.append("gp_num.init");
  }
  PreParseFromFile(filepath, datdim, pardim, num_sampled_data, 
                   num_steps_analyt_distribution, 
                   lower_limit_analyt_distribution,
                   upper_limit_analyt_distribution);

	// Initializes the Gaussian linear models for all three data (sub-)sets.
	LinGauss model_1(datdim/2 , pardim);
	LinGauss model_2(datdim/2, pardim);
	LinGauss model_tot(datdim, pardim);

  // Parses and stores the content of the vectors and matrices.
  PostParseFromFile(filepath, model_1.covar_matrix_, model_1.data_mean_values_, 
                    model_1.prior_param_values_, model_1.prior_matrix_,
                    model_1.g_matrix_
                   );
  gsl_matrix_memcpy(model_2.g_matrix_, model_1.g_matrix_);
  gsl_matrix_memcpy(model_2.covar_matrix_, model_1.covar_matrix_);
  gsl_matrix_memcpy(model_2.prior_matrix_, model_1.prior_matrix_);
  gsl_matrix_memcpy(model_tot.prior_matrix_, model_1.prior_matrix_);
  gsl_vector_memcpy(model_2.data_mean_values_, model_1.data_mean_values_);
  gsl_vector_memcpy(model_2.prior_param_values_, model_1.prior_param_values_);
  gsl_vector_memcpy(model_tot.prior_param_values_, model_1.prior_param_values_);

	// Creates data covariance matrix of model_tot.	
	BlockDiagFromTwoGslMatrices(model_1.covar_matrix_, model_2.covar_matrix_,
						 model_tot.covar_matrix_, datdim, datdim);	
	// Creates data mean of model_tot.
	AppendTwoGslVectors(model_1.data_mean_values_, model_2.data_mean_values_,
						 model_tot.data_mean_values_);
	// Creates g_matrix of model_tot.
	ConcatenateTwoGslMatrices(model_1.g_matrix_, model_2.g_matrix_,
													  model_tot.g_matrix_, 0);
	try {
    gsl_linalg_cholesky_decomp1(model_1.covar_matrix_);       
    gsl_matrix_memcpy(model_1.covar_matrix_, model_2.covar_matrix_);
  }
  catch (...) {
    std::cerr<<"covar_matrix_ is not positive-definite. ";
    
    return 1;
  }

  // Checks whether Sigma is positive-definite: if not, this throws an error.
	gsl_linalg_cholesky_decomp1(model_1.covar_matrix_);				
	gsl_matrix_memcpy(model_1.covar_matrix_, model_2.covar_matrix_);
	
	model_1.initialized_ = true;
	model_2.initialized_ = true;
	model_tot.initialized_ = true;
	model_1.CalculateFisherMatrix();
	model_2.CalculateFisherMatrix();
	model_tot.CalculateFisherMatrix();

  // Information printed to stdout.
  PrintGslMatrix(model_1.fisher_matrix_, "Fisher matrix");
  PrintGslVector(model_1.param_values_, "ML estimators");
  double detF= DeterminantGslMatrix(model_1.fisher_matrix_, pardim);
  double detP = DeterminantGslMatrix(model_1.prior_matrix_, pardim);
  std::cout<<"FoM ratio (prior vs Fisher):"<<std::endl;
  std::cout<<sqrt(detP/detF)<<std::endl;
	
	// Creates an instance "analyt" of the class "RobAnalyt".
	RobAnalyt analyt(num_steps_analyt_distribution,
									 upper_limit_analyt_distribution,
									 lower_limit_analyt_distribution, datdim, pardim,
									 model_tot.g_matrix_, model_tot.inv_covar_matrix_,
									 model_tot.fisher_matrix_, model_1.fisher_matrix_,
									 model_2.fisher_matrix_, model_tot.prior_matrix_,
									 model_tot.prior_param_values_, model_tot.data_mean_values_);

	std::string datapath;

  if (operation_mode == 20) {
		// Robustness distribution in nd weak prior case (numerical).
		int success = NumericalSamplingOfRobustness(num_sampled_data, datdim, true,
                                                model_1, model_2, model_tot,
                                                datapath);

		// Opens and writes to the metadata file.
		std::string metapath("metadata/");
		metapath.append(Timestamp());
		metapath.append("_20.rmeta");
  	std::ofstream metafile(metapath, std::ofstream::out | std::ofstream::app);

		WriteMetadata(metafile, datdim, model_1.covar_matrix_, pardim,
                  model_1.prior_param_values_, model_1.prior_matrix_,
                  model_1.g_matrix_, num_sampled_data);
		metafile.close();

		// Opens and writes to the logs file.
  	std::ofstream logfile("logs/robustness.log", std::ofstream::out 
  												| std::ofstream::app);
		WriteLogs(logfile, metapath, 20, datapath, success);
		logfile.close();

	} else if (operation_mode == 21) {
		// Robustness distribution in nd weak prior case (analytical).
    int success = analyt.CalculateAnalytWeakPrior(datapath);

    // Opens and writes to the metadata file.
    std::string metapath("metadata/");
    metapath.append(Timestamp());
    metapath.append("_21.rmeta");
    std::ofstream metafile(metapath, std::ofstream::out | std::ofstream::app);

    WriteMetadata(metafile, datdim, model_1.covar_matrix_, pardim,
                  model_1.prior_param_values_, model_1.prior_matrix_,
                  model_1.g_matrix_, 0);
    metafile.close();

    // Opens and writes to the logs file.
    std::ofstream logfile("logs/robustness.log", std::ofstream::out 
                          | std::ofstream::app);
    WriteLogs(logfile, metapath, 21, datapath, success);
    logfile.close();

	} else if (operation_mode == 30) {
		// Robustness distribution in nd case (numerical).
		int success = NumericalSamplingOfRobustness(num_sampled_data, datdim, false,
                                                model_1, model_2, model_tot,
                                                datapath);

		// Opens and writes to the metadata file.
		std::string metapath("metadata/");
		metapath.append(Timestamp());
		metapath.append("_30.rmeta");
  	std::ofstream metafile(metapath, std::ofstream::out | std::ofstream::app);

		WriteMetadata(metafile, datdim, model_1.covar_matrix_, pardim,
                  model_1.prior_param_values_, model_1.prior_matrix_,
                  model_1.g_matrix_, num_sampled_data);
		metafile.close();

		// Opens and writes to the logs file.
  	std::ofstream logfile("logs/robustness.log", std::ofstream::out 
  												| std::ofstream::app);
		WriteLogs(logfile, metapath, 30, datapath, success);
		logfile.close();

	} else {
		std::cout<<"No valid operation_mode selected.";
		abort();
	}

	model_1.DisallocateMembers();
	model_2.DisallocateMembers();
	model_tot.DisallocateMembers();																																		

  // Plots the distribution with python.
	std::string plot_command = (std::string)"python source/plot_robustness.py " 
                             + datapath;
	system(plot_command.c_str());
	
	return 0;
}
