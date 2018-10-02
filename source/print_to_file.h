#ifndef SOURCE_PRINTTOFILE_H_
#define SOURCE_PRINTTOFILE_H_

// *** release/print_to_file.h ***
// Author: Kevin Wolz, date: 08/2018
//
// Utilities to write the results obtained by robustness_main into files. 

#include <iostream>
#include <iomanip> 
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
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "utils.h"
#include "md5hash.h"


// Appends the robustness (and, if analytical output, pdf) value to a file.
void WriteData(std::ofstream& datafile, double& robustness, double& pdf);

// Appends the initial conditions of a run (e.g. dimensions of data, covariance
// matrix etc.) to a file.
void WriteMetadata(std::ofstream& metafile, const unsigned& datdim, 
                   gsl_matrix* covar_matrix, const unsigned& pardim, 
                   gsl_vector* prior_param_values, gsl_matrix* prior_matrix,
                   gsl_matrix* g_matrix, const unsigned& num_samplings);

// Appends the metadata about a run (e.g. data hash ID, filepaths, timestamp, 
// etc) to a file.
void WriteLogs(std::ofstream& logfile, std::string& meta_filepath, 
               const int& op_mode, std::string& data_filepath, 
               const int& success);

#endif // SOURCE_PRINTTOFILE_H_
