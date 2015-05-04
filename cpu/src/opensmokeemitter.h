#ifndef SMOKEEMITTER_H
#define SMOKEEMITTER_H

#include "includes.h"
#include "particle.h"
#include "gasconstraint.h"

#define H 2.
#define H2 4.
#define H6 64.
#define H9 512.

class OpenSmokeEmitter
{
public:
    OpenSmokeEmitter(glm::dvec2 posn, double particlesPerSec, GasConstraint *gs);
    virtual ~OpenSmokeEmitter();
    void tick(QList<Particle *> *estimates, double secs);
    QList<Particle *> *getParticles();
    inline glm::dvec2 getPosn() { return m_posn; }

private:
    double poly6(double r2);
    glm::dvec2 spikyGrad(const glm::dvec2 &r, double rlen2);

    glm::dvec2 m_posn;
    double m_particlesPerSec;
    QList<Particle *> m_particles;
    double timer;
    GasConstraint *m_gs;
};

#endif // SMOKEEMITTER_H
