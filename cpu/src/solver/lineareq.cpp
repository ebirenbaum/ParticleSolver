/**<!-------------------------------------------------------------------->
   @file   LinearEquation.cpp
   @author Travis Fischer (fisch0920@gmail.com)
   @date   Spring 2009
   
   @brief
      Wrapper around UMFPack for efficiently solving sparse linear systems of 
   the standard form Ax=b, including a cache of the factorization of A=LU such 
   that repeated invocations of solving Ax=b with the same A=LU can be 
   performed in nearly linear time by back-substitution of LUx=b.
   
   @note for more info on UMFPack (the internal linear solver), see 
      http://www.cise.ufl.edu/research/sparse/umfpack
   <!-------------------------------------------------------------------->**/

#include "lineareq.h"
#include <UMFPACK/umfpack.h>
#include <algorithm>

struct SparseMatrixElement {
   unsigned row;
   unsigned col;
   double value;
   
   inline SparseMatrixElement() {
      row   = 0;
      col   = 0;
      value = 0;
   }

   inline SparseMatrixElement(unsigned r, unsigned c, double v) {
      row   = r;
      col   = c;
      value = v;
   }

   inline bool operator< (const SparseMatrixElement& e) const {
      if (col < e.col)
         return true;
      if (col > e.col)
         return false;

      return (row < e.row);
   }

   inline bool operator==(const SparseMatrixElement &e) const {
      return (row == e.row && col == e.col && value == e.value);
   }
};

void LinearData::clean() {
   umfpack_di_free_symbolic(&symbolic);
   umfpack_di_free_numeric(&numeric);
   if(Ap) delete[] Ap;
   if(Ai) delete[] Ai;
   if(Ax) delete[] Ax;

}

void LinearData::init(unsigned n_, unsigned nElements_) {
   clean();
   
   n = n_;
   nElements = nElements_;
   
   Ap = new int[n + 1];
   Ai = new int[nElements];
   Ax = new double[nElements];
}

bool LinearData::init(const SparseMatrix &A) {
   const SparseArray &elements = A.getData();
   vector<SparseMatrixElement> temp;
   SparseColArrayConstIterator j;
   SparseArrayConstIterator i;
   
   for(i = elements.begin(); i != elements.end(); ++i) {
      for(j = i->second.begin(); j != i->second.end(); ++j) {
         SparseMatrixElement element(i->first, j->first, j->second);
         temp.push_back(element);
      }
   }
   
   std::sort(temp.begin(), temp.end());
   init(A.getN(), temp.size());
   
   Ap[0]	= 0;
   Ap[1]	= 0;
   
   unsigned cur   = 1;
   unsigned index = 0;
   
   // convert the sparse data into the format umfpack expects
   for(unsigned i = 0; i < temp.size(); ++i) {
      const SparseMatrixElement &element = temp[i];
      
      Ax[index] = element.value;
      Ai[index] = element.row;
      ++index;
      
      if (element.col == cur) {
         ++cur;
         Ap[cur] = Ap[cur - 1];
      }
      
      ++Ap[cur];
   }
   
   int status = UMFPACK_OK;
   bool ret   = true;
   
   status = umfpack_di_symbolic(n, n, Ap, Ai, Ax, &symbolic, NULL, NULL);
   if (status != UMFPACK_OK) {
      ret = false;
      printf("umfpack_di_symbolic failed: %d\n", status);
   }
   
   status = umfpack_di_numeric(Ap, Ai, Ax, symbolic, &numeric, NULL, NULL);
   if (status != UMFPACK_OK) {
      ret = false;
      printf("umfpack_di_numeric failed: %d\n", status);
   }
   
   return ret;
}

void LinearEquation::setA(const SparseMatrix *A) {
   m_A = A;
   m_dirty = true;
   
}

bool LinearEquation::solve(const double *b, double *x) {
   if (m_dirty) {
      if (!m_data.init(*m_A))
         return false;
      
      m_dirty = false;
   }
   
   int ret = umfpack_di_solve(UMFPACK_A, m_data.Ap, m_data.Ai, m_data.Ax, x, b, 
                              m_data.numeric, NULL, NULL);

   if (ret != UMFPACK_OK)
       printf("umfpack_di_solve failed: %d\n", ret);

   return ret == UMFPACK_OK;
}

