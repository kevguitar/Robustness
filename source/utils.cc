// *** source/utils.cc ***
// Author: Kevin Wolz, date: 08/2018
//
// Various utilities used in the code.

#include "utils.h"


// Creates a string of the current timestamp in the format "23-Aug-18_14-04-36".
std::string Timestamp() {
  time_t rawtime;
  struct tm * timeinfo;
  char timestamp [20];
  time (&rawtime);
  timeinfo = localtime (&rawtime);
  strftime(timestamp, 20, "%d-%b-%g_%H-%M-%S", timeinfo);

  return timestamp;
}


double get_signum(double x) {
	if (x > 0) return 1;
	if (x < 0) return -1;
	return 0;
}


void PrintGslVector(gsl_vector* vector) {
	int dim = vector->size;
	for (int x=0; x<dim; x++) {
			std::cout<<std::setprecision(4)<<std::setw(10)<<gsl_vector_get(vector, x)
      <<" ";
	}
	std::cout<<std::endl;
}

// Prints a gsl_vector* to std::cout, labeled by label.
void PrintGslVector(gsl_vector* vector, std::string label) {
  std::cout<<label<<":"<<std::endl;
  int dim = vector->size;
  for (int x=0; x<dim; x++) {
    std::cout<<std::setprecision(4)<<std::setw(10)<<gsl_vector_get(vector, x)
    <<" ";
  }
  std::cout<<std::endl;
}


void PrintGslMatrix(gsl_matrix* matrix) {
	int dim1 = matrix->size1;
	int dim2 = matrix->size2;
	for (int y=0; y<dim1; y++) {
		for (int x=0; x<dim2; x++)	{
			std::cout<<std::setprecision(4)<<std::setw(10)
      <<gsl_matrix_get(matrix, y, x)<<" ";
		}
		std::cout<<std::endl;
	}
}

// Prints a gsl_matrix* to std::cout, labeled by label.
void PrintGslMatrix(gsl_matrix* matrix, std::string label) {
	std::cout<<label<<":"<<std::endl;
	int dim1 = matrix->size1;
	int dim2 = matrix->size2;
	for (int y=0; y<dim1; y++) {
		for (int x=0; x<dim2; x++)	{
			std::cout<<std::setprecision(4)<<std::setw(10)
      <<gsl_matrix_get(matrix, y, x)<<" ";
		}
		std::cout<<std::endl;
	}
}


// Checks whether a file of a specific name exists inside a directory.
bool file_exists(const std::string& filename) {
  struct stat buffer;
  return (stat(filename.c_str(), &buffer) == 0);
}


// Concatenates two gsl_vectors*, creating a vector of length equal to
// the sum of the lenghts of the original vectors. The concatenated gsl_vector* 
// new_vector must be allocated externally with matching size.
void AppendTwoGslVectors(gsl_vector* vector1, gsl_vector* vector2, 
                         gsl_vector* new_vector) {
	int d1 = vector1->size;
	int d2 = vector2->size;
	int dn = new_vector->size;
	gsl_vector* intermed = gsl_vector_calloc(dn);
	
	if ((d1 + d2) != dn)	{
		std::cout<<"@ "<<__FUNCTION__<<": new_vector does not have the right \
                dimension. Please allocate correctly.";
	  abort();
	}
	for (int i=0; i<d1; i++) {
		double hilf = gsl_vector_get(vector1, i);
		gsl_vector_set(intermed, i, hilf);
	}
	for (int i=0; i<d2; i++) {
		double hilf = gsl_vector_get(vector2, i);
		gsl_vector_set (intermed, d1 + i, hilf);
	}
	gsl_vector_memcpy(new_vector, intermed);
}


// Concatenate two gsl_matrices vertically or horizontally. The concatenated
// matrix mnew must be allocated externally, and with matching sizes. 
// Note: vertical_axis = 1 refers to vertical concatenation, and 
//       vertical_axis = 0 refers to horizontal concatenation.
void ConcatenateTwoGslMatrices(gsl_matrix* matrix1, gsl_matrix* matrix2,
                               gsl_matrix* new_matrix, bool vertical_axis) {
	int m11 = matrix1->size1;
	int m12 = matrix1->size2;
	int m21 = matrix2->size1;
	int m22 = matrix2->size2;
	int mn1 = new_matrix->size1;
	int mn2 = new_matrix->size2;
	gsl_matrix* intermed = gsl_matrix_calloc(mn1, mn2);
	
	if (vertical_axis) {
		if ((m11 + m21) != mn1)	{
			std::cout<<"@ "<<__FUNCTION__<<": new_matrix does not have the right \
                  y-dimension. Note that vertical_axis = "<<vertical_axis<<". \
                  Please allocate correctly.";
			abort();
		}
		if ((m12 != m22) || (m22 != mn2))	{
			std::cout<<"@ "<<__FUNCTION__<<": Provided matrices do not have equal \
			            x-dimensions. Note that vertical_axis = "<<vertical_axis<<". \
                  Please allocate correctly.";
		  abort();
		}
		for (int j=0; j<mn2; j++) {
			for (int i=0; i<m11; i++) {
				double hilf = gsl_matrix_get(matrix1, i, j);
				gsl_matrix_set (intermed, i, j, hilf);
			}
			for (int i=0; i<m21; i++) {
				double hilf = gsl_matrix_get(matrix2, i, j);
				gsl_matrix_set(intermed, m11 + i, j, hilf);
			}
		}
	} else {
		if (((m12 + m22) != mn2))	{
			std::cout<<"@ "<<__FUNCTION__<<": new_matrix does not have the right \
                  x-dimension. Note that vertical_axis = "<<vertical_axis<<". \
                  Please allocate correctly.";
			abort();
		}
		if ((m11 != m21) || (m21 != mn1))	{
			std::cout<<"@ "<<__FUNCTION__<<": Provided matrices do not have equal \
			            y-dimensions. Note that vertical_axis = "<<vertical_axis<<". \
                  Please allocate correctly.";
		  abort();
		}
		
		for (int i=0; i<mn1; i++) {
			for (int j=0; j<m12; j++) {
				double hilf = gsl_matrix_get(matrix1, i, j);
				gsl_matrix_set (intermed, i, j, hilf);
			}
			for (int j=0; j<m22; j++) {
				double hilf = gsl_matrix_get(matrix2, i, j);
				gsl_matrix_set (intermed, i, m12+j, hilf);
			}
		}
	}	
	gsl_matrix_memcpy(new_matrix, intermed);
}


// Takes the gsl_matrices matrix1 and matrix2 and puts them into a 
// block-diagonal gsl_matrix new_matrix which must be allocated externally.
void BlockDiagFromTwoGslMatrices(gsl_matrix* matrix1, gsl_matrix* matrix2,
                                 gsl_matrix* new_matrix, unsigned int dimnew1,
                                 unsigned int dimnew2) {
	if ((new_matrix->size1 != dimnew1) || (new_matrix->size2 != dimnew2)) {
	  std::cout<<"@ "<<__FUNCTION__<<": Block matrix has wrong dimensions. You \
                have to allocate correctly. ";
	  abort();
	}
	gsl_matrix* intermed = gsl_matrix_calloc(dimnew1, dimnew2);
	int m11 = matrix1->size1;
	int m12 = matrix1->size2;
	int m21 = matrix2->size1;
	int m22 = matrix2->size2;
	
	
	if ((m11 + m21) < dimnew1) {
		std::cout<<"@ "<<__FUNCTION__<<": dimnew1 < dim1(matrix1) + dim1(matrix2).";
		abort();
	}
	if ((m12 + m22) < dimnew2) {
		std::cout<<"@ "<<__FUNCTION__<<": dimnew2 < dim2(matrix1) + dim2(matrix2).";
		abort();
	}
  // Requires dimensions dimnew1 >= dim1(matrix1) + dim1(matrix2), 
  // dimnew2 >= dim2(matrix1) + dim2(matrix2).
	for (int i=0; i<m11; i++)	{
		for (int j=0; j<m12; j++)	{
			double hilf = gsl_matrix_get(matrix1, i, j);
			gsl_matrix_set(intermed, i, j, hilf);
		}
	}	
	for (int i=0; i<m21; i++)	{
		for (int j=0; j<m22; j++)	{
			double hilf = gsl_matrix_get(matrix2, i, j);
			gsl_matrix_set(intermed, m11 + i, m12 + j, hilf);
		}
	}				
	gsl_matrix_memcpy(new_matrix, intermed);
	gsl_matrix_free (intermed);
}


void SubtractTwoDoubleVectors(const std::vector<double> &vector1,
                              const std::vector<double> &vector2,
                              std::vector<double> &difference) {
	size_t n = vector1.size();
  difference.resize(n);
	if (n == vector2.size())	{
		for (size_t i = 0; i < n; ++i) {
			difference[i] = vector1[i] - vector2[i];	
		};
	}
	else { 
	  std::cout<<"@ "<<__FUNCTION__<<": vector sizes don't match. "; abort();
	}		
}

// Fully contracts two gsl_vectors and a gsl_matrix, outputting the result as a
// double.
double GslVec1MatVec2(gsl_vector* vector1, gsl_matrix* matrix,
                      gsl_vector* vector2, int dim) {
  gsl_vector* hilf = gsl_vector_alloc(dim);
  double scalar_result;
  // matrix*vector2
  gsl_blas_dgemv(CblasNoTrans, 1.0, matrix, vector2, 0.0, hilf);  
  gsl_blas_ddot(hilf, vector1, &scalar_result);
  gsl_vector_free(hilf);
  return scalar_result;
}


// Takes the matrix "tobeinverted", which is of dimension "dim", and inverts it.
// The result is stored in an externally provided gsl_matrix "inverted".
void InvertGslMatrix(gsl_matrix* tobeinverted, unsigned int dim, 
                     gsl_matrix* inverted) {
  // Do not allocate space for the resulting matrix within this function!
  // The code would compile, and store the matrix here - but you want to use
  // it outside this function. However, the outside allocated matrix will then
  // not be the one you just calculated!
  gsl_matrix* hilf = gsl_matrix_alloc(dim,dim);

  // Use memcopy, such that the content of the matrix is copied. Else, the 
  // pointers would be set equal.
  gsl_matrix_memcpy(hilf, tobeinverted);
  int signum;
  gsl_permutation * perm = gsl_permutation_alloc(dim);
  gsl_linalg_LU_decomp(hilf, perm, &signum);
  gsl_linalg_LU_invert(hilf, perm, inverted);
  gsl_matrix_free(hilf);
  gsl_permutation_free(perm);
}


// Calculates the determinant of a square gsl_matrix. The latter needs to be 
// externally allocated with sizes (dim, dim).
double DeterminantGslMatrix(gsl_matrix* square_matrix, unsigned int dim) {
  gsl_matrix* hilf = gsl_matrix_alloc(dim, dim);
  gsl_matrix_memcpy(hilf, square_matrix); 
  int signum;
  gsl_permutation * perm = gsl_permutation_alloc(dim);
  gsl_linalg_LU_decomp(hilf, perm, &signum);
  double det = gsl_linalg_LU_det(hilf, signum);
  gsl_permutation_free(perm);
  gsl_matrix_free(hilf);
  return det;
}


// Fully contracts two c++ vectors of doubles and a gsl_matrix*. The latter one
// needs to be externally allocated with size (dim, dim).
double Vec1MatVec2(std::vector<double> vector1, gsl_matrix* square_matrix,
                   std::vector<double> vector2, int dim) {
  double res = 0;
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      res += vector1[i] * gsl_matrix_get(square_matrix, i, j) * vector2[j];
    }
  }
  return res;
}


// Calculates the trace of a square gsl_matrix*. The latter one needs to be
// externally allocated with size (dim, dim).
double TraceGslMatrix(gsl_matrix* square_matrix, int dim) {
  double tr = 0.0;
  for (int i = 0; i < dim; i++) {
    tr += gsl_matrix_get(square_matrix, i, i);
  }
  return tr;
}

// Calculates the sum of four square gsl_matrices of size (dim, dim) and stores
// the result in an externally allocated matrix.
void MultiplyFourGslMatrices(gsl_matrix* matrix1, gsl_matrix* matrix2,
                             gsl_matrix* matrix3, gsl_matrix* matrix4,
                             int dim, gsl_matrix* resulting_matrix) {
  gsl_matrix* intermed1 = gsl_matrix_calloc(dim, dim);
  gsl_matrix* intermed2 = gsl_matrix_calloc(dim, dim);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, matrix3, matrix4, 0.0,
                 intermed1);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, matrix2, intermed1, 0.0, 
                 intermed2);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, matrix1, intermed2, 0.0,
                 resulting_matrix);
  gsl_matrix_free(intermed1);
  gsl_matrix_free(intermed2);
}

// Calculates the sum of three square gsl_matrices of size (dim, dim) and stores
// the result in an externally allocated matrix.
void MultiplyThreeGslMatrices(gsl_matrix* matrix1, gsl_matrix* matrix2,
                              gsl_matrix* matrix3, int dim, 
                              gsl_matrix* resulting_matrix) {
  gsl_matrix* intermed1 = gsl_matrix_calloc(dim, dim);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, matrix2, matrix3, 0.0,
                 intermed1);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, matrix1, intermed1, 0.0,
                 resulting_matrix);
  gsl_matrix_free(intermed1);
}

