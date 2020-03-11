/***** written by Thomas Schoenemann as a private person, March 2020 *****/

#ifndef MATRIX_DETERMINANT_HH
#define MATRIX_DETERMINANT_HH

#include "matrix.hh"

//matrix determinant via LU-decomposition. It works for any mathematical T thanks to internally working with doubles.
// But since determinants only exist for square matrices, we do not make this a member.
template<typename T, typename ST>
long double matrix_determinant(const Math2D::Matrix<T,ST>& m, const double fullrank_tolerance = 0.0);

/****** implementation ****/

//matrix determinant via LU-decomposition
template<typename T, typename ST>
long double matrix_determinant(const Math2D::Matrix<T,ST>& m, const double fullrank_tolerance)
{
  if (m.xDim() != m.yDim()) {
    INTERNAL_ERROR << " No determinant for non-square matrix \"" << m.name() << "\" of size "
                   << m.xDim() << "x" << m.yDim() << ". Exiting..." << std::endl;
  }
  
  const ST dim = m.xDim();
  Math2D::Matrix<double,ST> lu(dim,dim);
  
  Math1D::Vector<ST,ST> true_row(dim);
  for (ST i = 0; i < dim; i++)
    true_row[i] = i;

  long double det = 1.0;

  for (ST r = 0; r < dim; r++) {
    
    //update u
    for (ST c = 0; c < r; c++) {
 
      const T* l_row = lu.row_ptr(c);
      
      T sum = m(r, true_row[c]);
      for (ST cc = 0; cc < c; cc++)
        sum -= lu(r, cc) * l_row[cc];
      
      lu(r, c) = sum;
    }
    
    //update l (save final div)
    T best_val = 0.0;
    ST arg_best = r;
    for (ST c = r; c < dim; c++) {

      const T* l_row = lu.row_ptr(c);
     
      T sum = m(r, true_row[c]);
      for (ST cc = 0; cc < r; cc++)
        sum -= lu(r, cc) * l_row[cc];
      
      lu(r, c) = sum;
      
      if (Makros::abs(sum) > best_val) {
        best_val = Makros::abs(sum);
        arg_best = c;
      }
    }
    
    if (fabs(best_val) <= fullrank_tolerance)
      return 0.0;
    
    if (arg_best != r) {
      
      std::swap_ranges(lu.row_ptr(r),lu.row_ptr(r)+r+1,lu.row_ptr(arg_best));
      
      std::swap(true_row[r], true_row[arg_best]);
      det = -det;
    }
    
    //apply final div
    const T div = lu(r,r);
    for (ST c = r+1; c < dim; c++)
      lu(r, c) /= div;
    
    det *= div;
  }
  
  return det;
}

#endif