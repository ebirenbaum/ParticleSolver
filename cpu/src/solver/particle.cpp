#include "particle.h"

SDFData Particle::getSDFData(QList<Body *> *bodies, int idx)
{
    if (ph != SOLID || bod < 0) {
        return SDFData();
    }

    Body *body = bodies->at(bod);
    SDFData out = body->sdf[idx];
    out.rotate(body->angle);
    return out;
}

void Body::updateCOM(QList<Particle *> *estimates, bool useEstimates)
{
    // Recompute center of mass
    glm::dvec2 total;
    for (int i = 0; i < particles.size(); i++) {
        Particle *p = estimates->at(particles[i]);
        total += (useEstimates ? p->ep : p->p) / p->imass;
    }
    center = total * imass;

    // Recompute angle delta guess
    angle = 0.0;
    double prev = 0.0;
    for (int i = 0; i < particles.size(); i++) {
        int index = particles[i];
        glm::dvec2 q = rs[index];
        if (glm::dot(q,q) == 0) {
            continue;
        }
        Particle *p = estimates->at(index);
        glm::dvec2 r = p->ep - center;

        double cos = r.x * q.x + r.y * q.y,
               sin = r.y * q.x - r.x * q.y,
               next = atan2(sin, cos) / p->imass;

        // Ensure all guesses are close to each other
        if (i > 0) {
            if (prev - next >= M_PI) {
                next += 2 * M_PI;
            }
        } else {
            if (next < 0) {
                next += 2 * M_PI;
            }
        }

        angle += next;
        prev = next;
    }
    angle *= imass;
}

void Body::computeRs(QList<Particle *> *estimates)
{
    imass = 0.0;
    for (int i = 0; i < particles.size(); i++) {
        int idx = particles[i];
        Particle *p = estimates->at(idx);
        glm::dvec2 r = p->p - center;
        rs[idx] = r;

        if (glm::dot(r,r) != 0) {
            imass += (1.0 / p->imass);
        }
    }
    imass = 1.0 / imass;
}
