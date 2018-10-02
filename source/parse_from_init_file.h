#ifndef SOURCE_PARSEFROMINITFILE_H_
#define SOURCE_PARSEFROMINITFILE_H_

// *** source/parse_from_init_file.h ***
// Author: Kevin Wolz, date: 08/2018
//
// Parses initial configurations (scalars, vectors and matrices) from .init 
// file. The latter has to have a specific format; in particular, it has to
// start with a semicolon and is does not have to end with a semicolon. More
// details can be inferred from the files located in ../init.

#include<iostream> 
#include<fstream> // filestream
#include<string>
#include<sstream> // stringstream
#include<algorithm> // remove_if
#include<vector>

#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>


// This defines a template for splitting strings by delimiters. It outputs the
// result inside an Out object (e.g. std::vector).
template<typename Out>
void split(const std::string &s, char delim, Out result) {
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, delim)) *(result++) = item;
}

std::vector<double> split(const std::string &s, char delim);
void PostParseFromFile(std::string& path, gsl_matrix* data_covar_matrix,
                       gsl_vector* data_mean_vector, 
                       gsl_vector* prior_params_vector, gsl_matrix* prior_matrix,
                       gsl_matrix* g_matrix);
void PreParseFromFile(std::string& path, unsigned& data_dim, 
                      unsigned& params_dim, unsigned& number_of_samplings,
                      unsigned& number_of_bins, double& lower_limit,
                      double& upper_limit);

#endif // SOURCE_PARSEFROMINITFILE_H_