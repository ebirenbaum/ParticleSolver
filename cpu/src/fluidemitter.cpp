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
        Particle *p = estimates->at(m_fs->ps.at(i));
            //            std::cout << p << std::endl;
//            double lambda = m_fs->lambdas[i];
            //            std::cout << lambda << std::endl;
        //            if(lambda >= -.1 && glm::length(p->v) < .05 && glm::length(p->p - p->ep) < .05) {
        //            if(p->p.y >= 10 || fabs(p->p.x) >= 10 ) {
        if(glm::length(p->v) < .06 && p->p.y <= 5) {
            if(m_fs->lambdas[i] <= 0) {
                p->t -= 1;
                if(p->t <= 0) {
                    p->t = 0;

                    //                p->ph = SOLID;
                    //                Particle *newP = new Particle(p->p, 0, SOLID);
                    //                newP->v = p->v;
                    //            p->imass -= secs;
                    //            if(p->imass == 0)

                    p->imass = 0;
                    p->ph = SOLID;
                    p->ep = p->p;
                    p->v = glm::dvec2();
                    p->f = glm::dvec2();
                    //                grains.append(newP);
                    //                estimates->append(newP);
                    //                estimates->removeAt(m_fs->ps.at(i));
                    //                if(m_fs->ps.contains(i))
                    m_fs->removeParticle(i);
                    //                delete p;
                    //                p->imass = 1;
                    //                p->ph = SOLID;
                    //                m_fs->ps.removeAt(i);
                }
            } else {
                p->t += secs;
                if(p->t > 3) p->t = 3;
            }
            //        }
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
