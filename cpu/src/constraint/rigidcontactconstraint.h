#ifndef RIGIDCONTACTCONSTRAINT_H
#define RIGIDCONTACTCONSTRAINT_H

#include "particle.h"

class RigidContactConstraint : public Constraint
{
public:
    RigidContactConstraint(int first, int second, QList<Body *> *bodies, bool st = false);
    virtual ~RigidContactConstraint();

    void initBoundary(Particle *p1, Particle *p2);

    void project(QList<Particle *> *estimates, int *counts);
    void draw(QList<Particle *> *particles);

    double evaluate(QList<Particle *> *estimates);
    glm::dvec2 gradient(QList<Particle *> *estimates, int respect);
    void updateCounts(int *counts);

private:
    QList<Body *> *bods;
    glm::dvec2 n;
    double d;
    int i1, i2;
    bool stabile;
};

#endif // RIGIDCONTACTCONSTRAINT_H
