#ifndef DISTANCECONSTRAINT_H
#define DISTANCECONSTRAINT_H

#include "constraint.h"

class DistanceConstraint : public Constraint
{
public:
    DistanceConstraint(float distance);
    virtual ~DistanceConstraint();

    float solveConstraint(glm::vec2 *particles, int numParticles);
    glm::vec2 getGradient(glm::vec2 *particles, int numParticles);

private:
    float m_distance;

};

#endif // DISTANCECONSTRAINT_H
