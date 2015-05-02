/**<!-------------------------------------------------------------------->
   @class  LinearEquation
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

#ifndef LINEAREQ_H
#define LINEAREQ_H

#include "matrix.h"

/**
 * @brief
 *    Underlying data encapsulating UMFPack LU decomposition of a sparse 
 * matrix A, used internally by LinearEquation (you should not need to touch 
 * this)
 */
struct LinearData {
   int     n;
   int     nElements;
   
   int    *Ap;
   int    *Ai;
   double *Ax;
   
   void   *symbolic;
   void   *numeric;
   
   inline LinearData()
      : n(0), nElements(0), Ap(NULL), Ai(NULL), Ax(NULL), 
        symbolic(NULL), numeric(NULL)
   { }
   
   ~LinearData() {
      clean();
   }
   
   void clean();
   
   bool init(const SparseMatrix &A);
   void init(unsigned n_, unsigned nElements_);
};

class LinearEquation {
   public:
      typedef shared_ptr<LinearEquation> Ref;

      inline LinearEquation(const SparseMatrix *A)
         : m_A(A), m_dirty(true)
      {

      }
      
      LinearEquation()
         : m_A(NULL), m_dirty(true)
      { }
      
      /**
       * @brief
       *    Solves a linear equation of the standard form Ax=b, with x and b 
       * passed in and A having previously been set either in this 
       * LinearEquation's constructor or via the setA mutator.
       * 
       * Will compute the LU decomposition of A if this is the first call to 
       * solve since A has been modified (if A's underlying data has changed 
       * since the last call to solve, LinearEquation has no way of knowing 
       * this, and so if you want to dirty the cache, you must explicitly call 
       * setA to force the recomputation of A's LU cache, encapsulated in 
       * the LinearData struct)
       * 
       * @htmlonly
       * Example usage:
       * <pre><code>
       *    LinearEquation solver;
       *    double b[5], x[5];
       *    // fill in b...
       *    
       *    SparseMatrix A(5, 5);
       *    // fill in A
       *    
       *    // setup solver and solve Ax=b
       *    solver.setA(A);
       *    bool solved = solver.solve(b, x);
       *    ASSERT(solved);
       *    
       *    // solve Ax=c using the same factorization of A
       *    double c[5];
       *    // fill in c
       *    
       *    solved = solver.solve(c, x);
       *    ASSERT(solved);
       *    // ...
       * </code></pre>
       * @endhtmonly
       * 
       * @note if A is nxn (must be square for Ax=b to be well-defined), b and 
       *    x are both assumed to be preallocated to hold n doubles)
       * 
       * @returns whether or not Ax=b was successfully solved, and upon 
       *    success, the solution vector x will be filled in appropriately
       */
      virtual bool solve(const double *b, double *x);
      
      /**
       * @brief
       *    Sets 'A' in Ax=b and "dirties" the underlying LinearData cache 
       * which stores the LU factorization of A such that the next call to 
       * solve will implicitly factor the current A before solving
       * 
       * @note because setA dirties the factorization cache, you are advised 
       *    to change A only when absolutely necessary, as factoring A=LU is 
       *    the most time-consuming part of solving a sparse Ax=b
       * 
       * @note A must be a square matrix for Ax=b to be solvable
       */
      virtual void setA(const SparseMatrix *A);
      
      /**
       * @returns the currently set 'A' which will be used in calls to 
       *    solve(b, x) for Ax=b
       * 
       * @note getLinearData optionally returns the cached LU decomposition 
       *    of A in UMFPack format iff the cache is not dirty (iff solve has 
       *    been called since the last call to setA)
       */
      inline const SparseMatrix *getA() const {
         return m_A;
      }
      
      /**
       * @returns the underlying UMFPack cache, including the most recently 
       *    computed LU decomposition of A during a call to solve
       * 
       * @note this method is provided for convenience only and its use is in 
       *    no way necessary to solve Ax=b
       */
      inline const LinearData &getLinearData() const {
         return m_data;
      }
      
   protected:
      const SparseMatrix *m_A;
      
      LinearData m_data;
      bool       m_dirty;
};

#endif // LINEAREQ_H

