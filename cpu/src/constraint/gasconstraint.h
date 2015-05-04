#ifndef GASCONSTRAINT_H
#define GASCONSTRAINT_H

#define H 2.
#define H2 4.
#define H6 64.
#define H9 512.

// USE H = 4, density = .5, look for the lattice

// Epsilon in gamma correction denominator
#define RELAXATION .01

// Pressure terms
#define K_P .1
#define E_P 4
#define DQ_P .1

// Fluid-solid coupling constant
#define S_SOLID .5

#include "particle.h"
#include <QSet>

class GasConstraint : public Constraint
{
public:
    GasConstraint(double density, QList<int> *particles, bool open);
    virtual ~GasConstraint();

    void project(QList<Particle *> *estimates, int *counts);
    void draw(QList<Particle *> *particles);

    double evaluate(QList<Particle *> *estimates);
    glm::dvec2 gradient(QList<Particle *> *estimates, int respect);
    void updateCounts(int *counts);

    double poly6(double rlen);
    glm::dvec2 spikyGrad(const glm::dvec2 &r, double rlen);
    glm::dvec2 grad(QList<Particle *> *estimates, int k, int j);

    void addParticle(Particle *p, int index);

private:
    double p0;
    QList<int> ps;
    QList<int> *neighbors;
    int numParticles;
    glm::dvec2 *deltas;
    QHash<int, double> lambdas;
    bool m_open;
};

#endif // GASCONSTRAINT_H
