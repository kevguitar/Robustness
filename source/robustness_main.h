#ifndef SOURCE_ROBUSTNESSMAIN_H_
#define SOURCE_ROBUSTNESSMAIN_H_

// *** source/robustness_main.h ***
// Author: Kevin Wolz, date: 08/2018
//
// Numerical sampling of the Robustness distribution (GP and WPL).

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <stdlib.h>
#include <sys/stat.h>
#include <ctime>
#include <time.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h> 
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "utils.h"
#include "lin_gauss.h"
#include "rob_analyt.h"
#include "print_to_file.h"
#include "parse_from_init_file.h"



double RobustnessWeakPriorLimit(double chi_squared_1,     // subset 1 
                                double chi_squared_2,     // subset 2
                                double chi_squared_tot,   // total data
                                double det_fisher_1,      // determinants
                                double det_fisher_2,
                                double det_fisher_tot,
                                double det_prior);

double RobustnessGeneralPrior(double chi_squared_1,
                              double chi_squared_2,
                              double chi_squared_tot,
                              gsl_matrix* fisher_matrix_1,
                              gsl_matrix* fisher_matrix_2,
                              gsl_matrix* prior_matrix,
                              gsl_vector* param_values_1, // fiducial parameters
                              gsl_vector* param_values_2, 
                              gsl_vector* param_values_tot,
                              gsl_vector* prior_param_values,// prior parameters
                              unsigned int pardim);      // number of parameters

int NumericalSamplingOfRobustness(unsigned int repetitions, // number of samples
                                  unsigned int datdim,  // number of data points
                                  bool is_weak_prior_limit,
                                  LinGauss model1,          // subset 1
                                  LinGauss model2,          // subset 2
                                  LinGauss model_tot,       // total data
                                  std::string& filepath);   // path of data logs


#endif // SOURCE_ROBUSTNESSMAIN_H_
