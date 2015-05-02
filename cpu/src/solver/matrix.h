/* ---------------------------------------------- 
   file: SparseMatrix.h
   auth: Travis Fischer
   acct: tfischer
   date: Spring 2008

   Provides basic functionality for a sparse matrix of 
   arbitrary dimensions
   ---------------------------------------------- */

#ifndef MATRIX_H
#define MATRIX_H

#include <memory>
#include <iostream>
#include <sstream>
#include <vector>
#include <map>
using namespace std;

typedef map<int, double>               SparseColArray;
typedef SparseColArray::iterator       SparseColArrayIterator;
typedef SparseColArray::const_iterator SparseColArrayConstIterator;

typedef map<int, SparseColArray >   SparseArray;
typedef SparseArray::iterator       SparseArrayIterator;
typedef SparseArray::const_iterator SparseArrayConstIterator;

class SparseMatrix {
   public:
      typedef shared_ptr<SparseMatrix> Ref;

      // Constructors / Initialization
      // -----------------------------

      // Constructs an invalid matrix temporarily. This constructor lets you add
      // SparseMatrix as a member of your Solver class.
      //
      // If you initialize a SparseMatrix using this function, you *MUST* call
      // reset(m, n) before using any other method on this SparseMatrix, or the
      // results will be undefined.
      inline SparseMatrix();

      // Constructs an m by n sparse matrix
      inline SparseMatrix(int m, int n);
      
      inline SparseMatrix(int m, int n, const SparseArray &elements);
      
      // Copy Constructor
      inline SparseMatrix(const SparseMatrix &m);
      
      inline void reset(int m, int n);

      inline void reset();
      
      // Static convenience constructors to generate common matrices
      // -----------------------------------------------------------

      // Generates an empty sparse matrix of the desired dimension
      static inline SparseMatrix zero(int m, int n = -1) {
         if (n == -1)
            n = m; // square
         
         return SparseMatrix(m, n);
      }
      
      // Generates a sparse nxn identity matrix
      static inline SparseMatrix identity(int n) {
         SparseArray a;
         
         for(int i = n; i--;)
            a[i][i] = 1.0;
         
         return SparseMatrix(n, n, a);
      }
      
      
      // Accessor Operators
      // ------------------
      
      // Returns the number of rows in this SparseMatrix
      inline int getM() const;
      inline int numRows() const;
      // Returns the number of cols in this SparseMatrix
      inline int getN() const;
      inline int numCols() const;
      
      inline SparseArray &getData();
      inline const SparseArray &getData() const;
      
      
      // Equality Operators
      // ------------------
      inline       bool       operator==(const SparseMatrix &m) const;
      inline       bool       operator!=(const SparseMatrix &m) const;
      
      
      // Mutator Operators
      // -----------------
      inline       SparseMatrix &operator =(const SparseMatrix &m);
      inline       SparseMatrix &operator*=(const SparseMatrix &m);
      inline       SparseMatrix &operator+=(const SparseMatrix &m);
      inline       SparseMatrix &operator-=(const SparseMatrix &m);
      
      // Scalar Mutator Operators
      inline       SparseMatrix &operator*=(const double scale);
      inline       SparseMatrix &operator/=(const double scale);
      
      
      // Arithmetic Operators
      // -----------------
      inline       SparseMatrix  operator* (const SparseMatrix &rhs) const;
      inline       SparseMatrix  operator+ (const SparseMatrix &rhs) const;
      inline       SparseMatrix  operator- (const SparseMatrix &rhs) const;
      
      inline       SparseMatrix  operator* (const double scale) const;
      inline       SparseMatrix  operator/ (const double scale) const;
      
      
      // More Complex Functionality
      // --------------------------
      
      // Multiplies this MxN SparseMatrix by the dense mxn matrix 'rhs' and stores 
      // the dense result in 'out' which should be preallocated to hold Mxn doubles
      void multiply(double *out, const double *rhs, int m, int n);
      
      // Multiplies this MxN SparseMatrix by the dense mxn matrix 'rhs' and stores 
      // the result in the sparse Mxn matrix 'out'
      void multiply(SparseMatrix &out, const double *rhs, int m, int n);
      
      // Returns whether or not this SparseMatrix is square
      inline bool isSquare() const;
      
      // Returns whether or not this SparseMatrix is diagonal
      inline bool isDiagonal() const;
      
      // Returns the number of elements currently in this SparseMatrix
      inline int getSize() const;
      
      // Returns the full size in bytes of this SparseMatrix if it were to be 
      // expanded in non-sparse form
      inline int fullSize() const;
      
      // Stores a (row-major) full version of this SparseMatrix in paramater out, 
      // which should be a preallocated buffer of fullSize() bytes
      inline void toFull(double *out) const;
      
      // Sets the value in this SparseMatrix as specified
      inline SparseMatrix &setValue(int row, int col, double value);
      
      inline void setIdentityColumn(int col);
      
      // Returns whether or not this SparseMatrix contains a value for the 
      // given cell
      inline bool   hasValue(int row, int col);
      
      // Returns the value stored at (row, col) or zero if none exists
      inline double getValue(int row, int col);
      
      inline void clearRow(int row);
      
      // Returns the transpose of this matrix
      SparseMatrix getTranspose() const;
      
      // Returns the Frobenius Norm of this matrix
      double getFrobeniusNorm();
      
      // Cleans up a matrix (0's out entries that are less than epsilon)
      void cleanup();
      
      // Saves an image representation of this SparseMatrix to the file given, 
      // with black pixels representing data and white pixels representing 
      // sparse regions
      bool save(const std::string &filename) const;

      void printMatrix(int precision, bool sparse) {
          std::ostringstream oss;
          oss << "%." << precision << "f ";
          std::string format = oss.str();
          printf("\n");

          if (!sparse) {
              for (int r = 0; r < numRows(); r++) {
                  for (int c = 0; c < numCols(); c++) {
                      printf(format.data(), getValue(r,c));
                  }
                  printf("\n");
                  fflush(stdout);
              }
          } else {
              for (int r = 0; r < numRows(); r++) {
                  for (int c = 0; c < numCols(); c++) {
                      if (getValue(r,c) == 0.0 || getValue(r,c) == -0.0) {
                          continue;
                      }
                      printf("(%i,%i) ", r, c);
                      printf(format.data(), getValue(r,c));
                      printf("\n");
                      fflush(stdout);
                  }
              }
          }
      }
      
   protected:
      // Sparse representation of data
      SparseArray m_elements;
      
      // Dimensions
      int m_m, m_n;
};


// Extra operators where SparseMatrix is on right-hand side
// -----------------------------------------------------
inline       SparseMatrix  operator* (const double scale, const SparseMatrix &rhs);


/* Include inline implementations */
#include "matrix.inl"

#endif

