#ifndef BOUNDARYCONSTRAINT_H
#define BOUNDARYCONSTRAINT_H

#include "particle.h"

// Collision with the boundaries of the world
class BoundaryConstraint : public Constraint
{
public:
    BoundaryConstraint(int index, double val, bool xBoundary, bool greater, bool st = false);
    virtual ~BoundaryConstraint();

    void project(QList<Particle *> *estimates, int *counts);
    void draw(QList<Particle *> *particles);

    double evaluate(QList<Particle *> *estimates);
    glm::dvec2 gradient(QList<Particle *> *estimates, int respect);
    void updateCounts(int *counts);

private:
    int idx;
    double value;
    bool isX, isGreaterThan, stabile;
};

#endif // BOUNDARYCONSTRAINT_H
