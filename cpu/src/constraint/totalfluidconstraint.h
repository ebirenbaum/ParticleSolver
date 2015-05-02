#ifndef TOTALFLUIDCONSTRAINT_H
#define TOTALFLUIDCONSTRAINT_H

//#define H 3.5
//#define H6 1838.265625
//#define H9 78815.6386719

//#define H 5.
//#define H6 15625.
//#define H9 1953125.

//#define H 4.
//#define H6 4096.
//#define H9 262144.

#define H 2.
#define H2 4.
#define H6 64.
#define H9 512.

// Epsilon in gamma correction denominator
#define RELAXATION .01

// Pressure terms
#define K_P .1
#define E_P 4
#define DQ_P .2

// Fluid-solid coupling constant
#define S_SOLID 0.

#include "particle.h"

class TotalFluidConstraint : public Constraint
{
public:
    TotalFluidConstraint(double density, QList<int> *particles);
    virtual ~TotalFluidConstraint();

    void project(QList<Particle *> *estimates, int *counts);
    void draw(QList<Particle *> *particles);

    double evaluate(QList<Particle *> *estimates);
    glm::dvec2 gradient(QList<Particle *> *estimates, int respect);
    void updateCounts(int *counts);

    double poly6(double rlen);
    glm::dvec2 spikyGrad(const glm::dvec2 &r, double rlen);
    glm::dvec2 grad(QList<Particle *> *estimates, int k, int j);

private:
    double p0;
    QList<int> ps;
    QList<int> *neighbors;
    glm::dvec2 *deltas;
    QHash<int, double> lambdas;
};

#endif // TOTALFLUIDCONSTRAINT_H
