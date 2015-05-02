#ifndef DISTANCECONSTRAINT_H
#define DISTANCECONSTRAINT_H

#include "particle.h"

// Two particles must be exactly a certain distance away
class DistanceConstraint : public Constraint
{
public:
    DistanceConstraint(double distance, int first, int second, bool st = false);
    DistanceConstraint(int first, int second, QList<Particle *> *particles);
    virtual ~DistanceConstraint();

    void project(QList<Particle *> *estimates);
    void draw(QList<Particle *> *particles);

    double evaluate(QList<Particle *> *estimates);
    glm::dvec2 gradient(QList<Particle *> *estimates, int respect);
    void updateCounts(int *counts);

private:
    double d;
    int i1, i2;
    bool stabile;
};

#endif // DISTANCECONSTRAINT_H
