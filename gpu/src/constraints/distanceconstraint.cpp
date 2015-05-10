#include "distanceconstraint.h"
#include <assert.h>

DistanceConstraint::DistanceConstraint(float distance)
{
    m_distance = distance;
}


DistanceConstraint::~DistanceConstraint()
{
}


float DistanceConstraint::solveConstraint(glm::vec2 *particles, int numParticles)
{
    assert(numParticles == 2);
    return glm::length(particles[0] - particles[1]) - m_distance;
}

glm::vec2 DistanceConstraint::getGradient(glm::vec2 *particles, int numParticles)
{
    assert(numParticles == 2);
    glm::vec2 diff = particles[0] - particles[1];
    float diffMag2 = glm::dot(diff, diff);
    if (diffMag2 > 0.00001f)
        return diff * float(1.f / sqrt(diffMag2));
    return glm::vec2();
}

