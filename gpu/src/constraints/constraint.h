#ifndef CONSTRAINT_H
#define CONSTRAINT_H

#include <glm.hpp>

class Constraint
{
public:
    Constraint() {}
    virtual ~Constraint() {}

    virtual float solveConstraint(glm::vec2 *particles, int numParticles) = 0;
    virtual glm::vec2 getGradient(glm::vec2 *particles, int numParticles) = 0;
};

#endif // CONSTRAINT_H
