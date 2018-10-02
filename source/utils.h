#ifndef SOURCE_UTILS_H_
#define SOURCE_UTILS_H_

// *** source/utils.h ***
// Author: Kevin Wolz, date: 08/2018
//
// Various utilities used in the code.

#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <string>
#include <unistd.h>
#include <algorithm>
#include <stdlib.h>
#include <sys/stat.h>
#include <time.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_cdf.h>


std::string Timestamp();

double get_signum(double x);

void PrintGslVector(gsl_vector* vector);

void PrintGslVector(gsl_vector* vector, std::string label);

void PrintGslMatrix(gsl_matrix* matrix);

void PrintGslMatrix(gsl_matrix* matrix, std::string label);

bool file_exists(const std::string& filename);

void AppendTwoGslVectors(gsl_vector* vector1, gsl_vector* vector2, 
                         gsl_vector* new_vector);

void ConcatenateTwoGslMatrices(gsl_matrix* mmatrix1, gsl_matrix* matrix2,
                               gsl_matrix* new_matrix, bool vertical_axis);

void BlockDiagFromTwoGslMatrices(gsl_matrix* matrix1, gsl_matrix* matrix2,
                                 gsl_matrix* new_matrix, unsigned int dimnew1,
                                 unsigned int dimnew2);

void SubtractTwoDoubleVectors(const std::vector<double> &vector1,
                              const std::vector<double> &vector2,
                              std::vector<double> &difference);

double GslVec1MatVec2(gsl_vector* vector1, gsl_matrix* matrix,
                      gsl_vector* vector2, int dim);

void InvertGslMatrix(gsl_matrix* tobeinverted, unsigned int dim, 
                     gsl_matrix* inverted);

double DeterminantGslMatrix(gsl_matrix* square_matrix, unsigned int dim);

double Vec1MatVec2(std::vector<double> vector1, gsl_matrix* matrix,
                   std::vector<double> vector2, int dim);

double TraceGslMatrix(gsl_matrix* square_matrix, int dim);

void MultiplyFourGslMatrices(gsl_matrix* matrix1, gsl_matrix* matrix2,
                             gsl_matrix* matrix3, gsl_matrix* matrix4,
                             int dim, gsl_matrix* resulting_matrix);

void MultiplyThreeGslMatrices(gsl_matrix* matrix1, gsl_matrix* matrix2,
                              gsl_matrix* matrix3, int dim, 
                              gsl_matrix* resulting_matrix);

#endif // SOURCE_UTILS_H_
