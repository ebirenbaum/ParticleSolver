#ifndef SOLVER_H
#define SOLVER_H

#include "matrix.h"
#include "particle.h"
#include "lineareq.h"

#define RELAXATION_PARAMETER 1.

class Solver
{
public:
    Solver();
    virtual ~Solver();

    SparseMatrix m_invM, m_JT, m_A;
    double *m_b, *m_gamma, *m_dp;
    int *m_counts;
    int m_nParts, m_nCons;
    LinearEquation m_eq;

    int getCount(int idx);

    void setupM(QList<Particle *> *particles, bool contact = false);
    void setupSizes(int numParts, QList<Constraint *> *constraints);
    void solveAndUpdate(QList<Particle *> *particles, QList<Constraint *> *constraints, bool stabile = false);
};

#endif // SOLVER_H
