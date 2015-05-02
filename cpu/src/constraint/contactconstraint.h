#ifndef CONTACTCONSTRAINT_H
#define CONTACTCONSTRAINT_H

#include "particle.h"

// Contact between two particles where AT LEAST ONE is not a solid
class ContactConstraint : public Constraint
{
public:
    ContactConstraint(int first, int second, bool st = false);
    virtual ~ContactConstraint();

    void project(QList<Particle *> *estimates);
    void draw(QList<Particle *> *particles);

    double evaluate(QList<Particle *> *estimates);
    glm::dvec2 gradient(QList<Particle *> *estimates, int respect);
    void updateCounts(int *counts);

private:
    int i1, i2;
    bool stabile;
};

#endif // CONTACTCONSTRAINT_H
