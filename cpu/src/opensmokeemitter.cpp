#include "opensmokeemitter.h"

OpenSmokeEmitter::OpenSmokeEmitter(glm::dvec2 posn, double particlesPerSec, GasConstraint *gs) :
    m_posn(posn), m_particlesPerSec(particlesPerSec), m_gs(gs)
{
    timer = 0;
}

OpenSmokeEmitter::~OpenSmokeEmitter() {
//    for(int i = m_particles.size()-1; i >= 0; i--) {
//        Particle *p = m_particles.at(i);
//        m_particles.removeAt(i);
//        delete(p);
//    }
}

void OpenSmokeEmitter::tick(QList<Particle *> *estimates, double secs) {
    timer += secs;
    while(timer >= 1./m_particlesPerSec) {
        timer -= 1./m_particlesPerSec;
        Particle *p = new Particle(m_posn, .1, GAS);
        m_particles.append(p);
        if(m_gs != NULL) {
            p = new Particle(m_posn, 1, GAS);
            m_gs->addParticle(p, estimates->size());
            estimates->append(p);
        }
    }
    for(Particle *p: m_particles) {
        if(p->ph == FLUID || p->ph == GAS) {
            p->v = glm::dvec2();
            double sum = 0;
            for(Particle *n: *estimates) {
                glm::dvec2 r = p->p - n->p;
                double p6 = poly6(glm::dot(r,r));
                p->v += n->v * p6;
                sum += p6;
            }

            if(sum > 0)
                p->p += p->v * secs / sum;
        }
    }
}

QList<Particle *> *OpenSmokeEmitter::getParticles()
{
    return &m_particles;
}

double OpenSmokeEmitter::poly6(double r2)
{
    if(r2 >= H2) return 0;
    double term2 = (H2 - r2);
    return (315. / (64. * M_PI * H9)) * (term2 * term2 * term2);
}

glm::dvec2 OpenSmokeEmitter::spikyGrad(const glm::dvec2 &r, double rlen2)
{
    if(rlen2 >= H) return glm::dvec2();
    if(rlen2 == 0) return glm::dvec2();
    return -glm::normalize(r) * (45. / (M_PI * H6)) * (H - rlen2) * (H - rlen2);
//    return -r / (H*H*rlen);
}
