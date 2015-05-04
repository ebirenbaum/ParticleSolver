#include "fluidemitter.h"

FluidEmitter::FluidEmitter(glm::dvec2 posn, double particlesPerSec, TotalFluidConstraint *fs) :
    m_posn(posn), m_particlesPerSec(particlesPerSec), m_fs(fs)
{
    timer = 0;
    totalTimer = 0;
}

FluidEmitter::~FluidEmitter() {
}

void FluidEmitter::tick(QList<Particle *> *estimates, double secs) {
    for(int i = m_fs->ps.size()-1; i >= 0; i--) {
        if(estimates->size() > m_fs->ps.at(i)) {
            Particle *p = estimates->at(m_fs->ps.at(i));
            std::cout << p << std::endl;
            double lambda = m_fs->lambdas[i];
            std::cout << lambda << std::endl;
            if(lambda >= 0 && glm::length(p->v) < .2) {
                Particle *newP = new Particle(p->p, 1, SOLID);
                newP->v = p->v;
                grains.append(newP);
                estimates->append(newP);
                estimates->removeAt(m_fs->ps.at(i));
                m_fs->removeParticle(i);
                delete p;
//                p->imass = 1;
//                p->ph = SOLID;
//                m_fs->ps.removeAt(i);
            }
        }
    }


//    for(int i=0; i<grains.size(); i++) {
//        Particle *p = grains.at(i);
//        if(glm::length(p->v) <= .02) {
//            grains.removeAt(i);
//            p->imass = 0;
//            p->v = glm::dvec2();
//            p->f = glm::dvec2();
//            p->ep = p->p;
//        }
//    }

    timer += secs;
    totalTimer += secs;
    while(totalTimer < 5 && timer >= 1./m_particlesPerSec) {
        timer -= 1./m_particlesPerSec;
        if(m_fs != NULL) {
            Particle *p = new Particle(m_posn, 1, FLUID);
            p->v = glm::dvec2(frand(),1);
            m_fs->addParticle(estimates->size());
            estimates->append(p);
        }
    }
}
