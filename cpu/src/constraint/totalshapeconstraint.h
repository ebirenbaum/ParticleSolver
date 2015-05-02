#ifndef TOTALSHAPECONSTRAINT_H
#define TOTALSHAPECONSTRAINT_H

#include "particle.h"

class TotalShapeConstraint : public Constraint
{
public:
    TotalShapeConstraint(Body *bod, double stiff = 1.0);
    virtual ~TotalShapeConstraint();

    void project(QList<Particle *> *estimates, int *counts);
    void draw(QList<Particle *> *particles);

    double evaluate(QList<Particle *> *estimates);
    glm::dvec2 gradient(QList<Particle *> *estimates, int respect);
    void updateCounts(int *counts);

    glm::dvec2 guess(int idx);

private:
    Body *body;
};

#endif // TOTALSHAPECONSTRAINT_H
