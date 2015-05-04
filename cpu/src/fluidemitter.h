#ifndef FLUIDEMITTER_H
#define FLUIDEMITTER_H

#include "includes.h"
#include "particle.h"
#include "totalfluidconstraint.h"

#define H 2.
#define H2 4.
#define H6 64.
#define H9 512.

class FluidEmitter
{
public:
    FluidEmitter(glm::dvec2 posn, double particlesPerSec, TotalFluidConstraint *fs);
    virtual ~FluidEmitter();
    void tick(QList<Particle *> *estimates, double secs);
    QList<Particle *> *getParticles();
    inline glm::dvec2 getPosn() { return m_posn; }

private:
    glm::dvec2 m_posn;
    double m_particlesPerSec;
    double timer;
    double totalTimer;
    TotalFluidConstraint *m_fs;
    QList<Particle *> grains;
};

#endif // FLUIDEMITTER_H
