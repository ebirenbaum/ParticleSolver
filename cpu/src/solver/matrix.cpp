/* ---------------------------------------------- 
   file: SparseMatrix.cpp
   auth: Travis Fischer
   acct: tfischer
   date: Spring 2008

   Provides basic functionality for a sparse matrix of 
   arbitrary dimensions
   ---------------------------------------------- */

#include "matrix.h"

#include <math.h>

// Returns the transpose of this matrix
SparseMatrix SparseMatrix::getTranspose() const {
   SparseArray m;
   
   SparseColArrayConstIterator j;
   SparseArrayConstIterator i;
   
   for(i = m_elements.begin(); i != m_elements.end(); ++i)
      for(j = i->second.begin(); j != i->second.end(); ++j)
         m[j->first][i->first] = j->second;
   
   // note m and n are switched
   return SparseMatrix(m_n, m_m, m);
}

// Returns the Frobenius Norm of this matrix
double SparseMatrix::getFrobeniusNorm() {
   double squaredNorm = 0.0;
   
   SparseColArrayConstIterator j;
   SparseArrayConstIterator i;
   
   for(i = m_elements.begin(); i != m_elements.end(); ++i)
      for(j = i->second.begin(); j != i->second.end(); ++j)
         squaredNorm += j->second * j->second;
   
   return sqrt(squaredNorm);
}

// Multiplies this MxN SparseMatrix by the dense mxn matrix 'rhs' and stores 
// the dense result in 'out' which should be preallocated to hold Mxn doubles
void SparseMatrix::multiply(double *out, const double *rhs, int m, int n) {
   
   memset(out, 0, sizeof(double) * m_m * n);
   SparseColArrayConstIterator j;
   SparseArrayConstIterator i;
   
   for(i = m_elements.begin(); i != m_elements.end(); ++i) {
      const int row = i->first;
      
      for(j = i->second.begin(); j != i->second.end(); ++j) {
         const int col = j->first;
         
         for(int k = n; k--;)
            out[row * n + k] += j->second * rhs[col * n + k];
      }
   }
}

// Multiplies this MxN SparseMatrix by the dense mxn matrix 'rhs' and stores 
// the result in the sparse Mxn matrix 'out'
void SparseMatrix::multiply(SparseMatrix &out, const double *rhs, int m, int n) {

   SparseColArrayConstIterator j;
   SparseArrayConstIterator i;

   for(i = m_elements.begin(); i != m_elements.end(); ++i) {
      const int row = i->first;

      for(j = i->second.begin(); j != i->second.end(); ++j) {
         const int col = j->first;

         for(int k = n; k--;) {
            double value = j->second * rhs[col * n + k];

            if (value != 0)
               out.setValue(row, k, out.getValue(row, k) + value);
         }
      }
   }
}

