// *** source/print_to_file.cc ***
// Author: Kevin Wolz, date: 08/2018
//
// Utilities to write the results obtained by robustness_main into files. 

#include "print_to_file.h"


// Writes data to file specified by std::ofstream datafile. Data consists of two
// columns, robustness and PDF, separated by a space. Requires filestream to be
// open().
void WriteData(std::ofstream& datafile, double& robustness, double& pdf) {
  if (!datafile.is_open()) { 
    std::cout<<"Datafile is not open.";
    abort();
  }
  datafile<<robustness<<" "<<pdf<<std::endl;
} 


// Writes metadata to file specified by std::ofstream metafile. Delimits every
// by a semicolon. Requires filestream to be open().
void WriteMetadata(std::ofstream& metafile, const unsigned& datdim, 
                   gsl_matrix* covar_matrix, const unsigned& pardim, 
                   /*gsl_vector* fid_param_values,*/
                   gsl_vector* prior_param_values, gsl_matrix* prior_matrix,
                   gsl_matrix* g_matrix, const unsigned& num_samplings) {
  if (!metafile.is_open()) {
    std::cout<<"Metafile is not open.";
    abort(); 
  }
  
  metafile<<"datdim "<<datdim<<";"<<std::endl;

  metafile<<"covar_matrix"<<std::endl;
  for (int i=0; i<(covar_matrix->size1); i++) { 
    for (int j=0; j<(covar_matrix->size2); j++) {
      metafile<<std::setprecision(5)<<gsl_matrix_get(covar_matrix, i, j)<<" ";
    }
    metafile<<std::endl;
  }
  metafile<<";"<<std::endl;

  metafile<<"pardim "<<pardim<<";"<<std::endl;

  metafile<<"prior_param_values"<<std::endl;
  for (int i=0; i<(prior_param_values->size); i++) {
    metafile<<std::setprecision(5)<<gsl_vector_get(prior_param_values, i)<<" ";
  }
  metafile<<";"<<std::endl;

  metafile<<"prior_matrix"<<std::endl;
  for (int i=0; i<(prior_matrix->size1); i++) { 
    for (int j=0; j<(prior_matrix->size2); j++) {
      metafile<<std::setprecision(5)<<gsl_matrix_get(prior_matrix, i, j)<<" ";
    }
    metafile<<";"<<std::endl;
  }

  metafile<<"g_matrix"<<std::endl;
  for (int i=0; i<(g_matrix->size1); i++) { 
    for (int j=0; j<(g_matrix->size2); j++) {
      metafile<<std::setprecision(5)<<gsl_matrix_get(g_matrix, i, j)<<" ";
    }
    metafile<<";"<<std::endl;
  }

  metafile<<"num_samplings "<<num_samplings<<";"<<std::endl;
}


// Writes information about the runs to a file represented by the open
// filestream logfile.
void WriteLogs(std::ofstream& logfile, std::string& meta_filepath, 
               const int& op_mode, std::string& data_filepath, 
               const int& success) {
  if (!logfile.is_open()) {
    std::cout<<"Logfile is not open.";
    abort(); 
  }

  std::string hash_id = HashFile(data_filepath);
  std::string timestamp = Timestamp();

  logfile<<hash_id<<" "<<timestamp<<" "<<op_mode<<" "<<data_filepath<<" "
         <<meta_filepath<<" "<<success<<std::endl;
}
