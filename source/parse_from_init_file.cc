// *** source/parse_from_init_file.h ***
// Author: Kevin Wolz, date: 08/2018
//
// Parses initial configurations (scalars, vectors and matrices) from .init 
// file. The latter has to have a specific format; in particular, it has to
// start with a semicolon and is does not have to end with a semicolon. More
// details can be inferred from the files located in ../init.

#include "parse_from_init_file.h"


// Splits a string by a certain delimiter and outputs the elements inside a
// std::vector<double>.
std::vector<double> split(const std::string &s, char delim) {
  std::vector<std::string> string_vector;
  split(s, delim, std::back_inserter(string_vector));
  std::vector<double> double_vector(string_vector.size());
  std::transform(string_vector.begin(), string_vector.end(),
                 double_vector.begin(), [](const std::string& val) {
    return std::stod(val);
  });
  
  return double_vector;
}

// Opens the file located at "path", reads and splits every non-empty line 
// (except scalars) at the first occurence of the delimiter "=" into two
// strings, "key" and "value". Prints out key and value via std::cout, and 
// stores the vector and matrix quantities.
void PostParseFromFile(std::string& path, gsl_matrix* data_covar_matrix,
                       gsl_vector* data_mean_vector, 
                       gsl_vector* prior_params_vector, gsl_matrix* prior_matrix,
                       gsl_matrix* g_matrix) {
  std::ifstream inFile;

  inFile.open(path);
  if (!inFile) {
    std::cerr << "Unable to open file "<<path<<". "<<std::endl;
    exit(1);   // Calls system to stop.
  }
  std::string line_content;
  std::string delimiter("=");
  std::string comment("//");
  std::string bracket("{");
  std::string two_bracket("{{");
  std::string row_end("},{");
  std::string two_close("}}");

  while (std::getline(inFile, line_content, ';')) {
    // Excludes comments and empty lines.
    auto pos = line_content.find(comment);
    if (pos != std::string::npos) line_content.erase(pos);
    if (!line_content.empty() /*&&
        line_content.compare(0, comment.length(), comment)*/) {
      // Removes empty spaces from line string and resets the container size.
      line_content.erase(remove_if(line_content.begin(),
                         line_content.end(), isspace), line_content.end());

      std::string key = line_content.substr(0, line_content.find(delimiter));
      std::string value = line_content.substr(line_content.find(delimiter) + 1,
                                              line_content.length());
      
      std::vector<std::vector<double>> matrix;

      if (!value.compare(0, 2, two_bracket)) { // 2d array
        std::cout<<key<<":"<<std::endl;
        
        std::string rest(value.substr(2, value.length()));
        
        if (rest.find(row_end)) {   // 2d array with more than one line
          while (rest.find(row_end) <= rest.find(two_close)) {
            std::vector<double> row_vec = split(
                rest.substr(0, rest.find(row_end)), ','
                );
            matrix.push_back(row_vec);
            std::string temp(rest.substr(rest.find(row_end) + 3, rest.length()));
            rest = temp;
          }      
        }
        std::vector<double> row_vec = split(
            rest.substr(0, rest.length() - 2), ','
            );
        matrix.push_back(row_vec);

        for (int rows = 0; rows<matrix.size(); rows++) {
          for (int cols = 0; cols<matrix[rows].size() - 1; cols++) {
            std::cout<<matrix[rows][cols]<<" ";
          }
          std::cout<<matrix[rows][matrix[rows].size() - 1]<<" "<<std::endl;
        }
        if (key == "data_covar_array") {
          for (int i = 0; i<data_covar_matrix->size1; i++) {
            for (int j = 0; j<data_covar_matrix->size2; j++) {
              gsl_matrix_set(data_covar_matrix, i, j, matrix[i][j]);
            }
          }
        }
        if (key == "g_matrix_array") {
          for (int i = 0; i<g_matrix->size1; i++) {
            for (int j = 0; j<g_matrix->size2; j++) {
              gsl_matrix_set(g_matrix, i, j, matrix[i][j]);
            }
          }
        }
        if (key == "prior_matrix") {
          for (int i = 0; i<prior_matrix->size1; i++) {
            for (int j = 0; j<prior_matrix->size2; j++) {
              gsl_matrix_set(prior_matrix, i, j, matrix[i][j]);
            }
          }
        }

      } else if (!value.compare(0, 1, bracket)) { // 1d array
        std::cout<<key<<":"<<std::endl;
        std::vector<double> vec = split(value.substr(1, value.length() - 1),
                                        ',');
        for (int i=0; i<vec.size(); i++) std::cout<<"  "<<vec[i];
        std::cout<<std::endl;

        if (key == "data_mean_array") {
          for (int i = 0; i<data_mean_vector->size; i++) {
              gsl_vector_set(data_mean_vector, i, vec[i]);
            }
        }
        if (key == "prior_params_array") {
          for (int i = 0; i<prior_params_vector->size; i++) {
              gsl_vector_set(prior_params_vector, i, vec[i]);
            }
        }

      } else { // scalars
        // Do nothing.
      }
    }
  }  
}


// Works like PostParseFromFile, but reads and stores only scalar quantites.
void PreParseFromFile(std::string& path, unsigned& data_dim, 
                      unsigned& params_dim, unsigned& number_of_samplings,
                      unsigned& number_of_bins, double& lower_limit,
                      double& upper_limit) {
  std::ifstream inFile;

  inFile.open(path);
  if (!inFile) {
    std::cerr << "Unable to open file "<<path<<". "<<std::endl;
    exit(1);   // Calls system to stop.
  }
  std::string line_content;
  std::string delimiter("=");
  std::string comment("//");
  std::string bracket("{");
  std::string two_bracket("{{");
  std::string row_end("},{");
  std::string two_close("}}");

  while (std::getline(inFile, line_content, ';')) {
    // Excludes comments and empty lines.
    auto pos = line_content.find(comment);
    if (pos != std::string::npos) line_content.erase(pos);
    if (!line_content.empty()) {
      // Removes empty spaces from line string and resets the container size.
      line_content.erase(remove_if(line_content.begin(),
                         line_content.end(), isspace), line_content.end());

      std::string key = line_content.substr(0, line_content.find(delimiter));
      std::string value = line_content.substr(line_content.find(delimiter) + 1,
                                              line_content.length());
      
      std::vector<std::vector<double>> matrix;

      if (!value.compare(0, 2, two_bracket)) { // 2d array
        // Do nothing.

      } else if (!value.compare(0, 1, bracket)) { // 1d array
        // Do nothing.

      } else { // scalars
        std::cout<<key<<":"<<std::endl;
        std::size_t size = value.length();
        double scalar = std::stod(value, &size);
        std::cout<<"  "<<scalar<<std::endl;

        if (key == "pardim") {
          params_dim = (unsigned) scalar;
        }
        if (key == "datdim") {
          data_dim = (unsigned) scalar;
        }
        if (key == "num_sampled_data") {
          number_of_samplings = (unsigned) scalar;
        }
        if (key == "num_steps_analyt_distribution") {
          number_of_bins = (unsigned) scalar;
        }
        if (key == "lower_limit") {
          lower_limit = (double) scalar;
        }
        if (key == "upper_limit") {
          upper_limit = (double) scalar;
        }
      }
    }
  }  
}
