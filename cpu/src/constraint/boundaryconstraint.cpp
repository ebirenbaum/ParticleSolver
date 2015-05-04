#include "boundaryconstraint.h"

BoundaryConstraint::BoundaryConstraint(int index, double val, bool xBoundary, bool greater, bool st)
    : Constraint(), idx(index), value(val), isX(xBoundary), isGreaterThan(greater), stabile(st)
{

}

BoundaryConstraint::~BoundaryConstraint()
{

}

void BoundaryConstraint::project(QList<Particle *> *estimates, int *counts)
{
    Particle *p = estimates->at(idx);

    // Add a little random jitter for fluids and gases so particles do not become trapped on boundaries
    double extra = p->ph == FLUID || p->ph == GAS ? frand() * .003 : 0;
    double d = (PARTICLE_RAD + extra);
    glm::dvec2 n = glm::dvec2();

    // Move the particle back into a valid spot (if necessary)
    if (isGreaterThan) {
        if (isX) {

            // Quit if no longer valid
            if (p->ep.x >= value + PARTICLE_RAD) {
                return;
            }
            p->ep.x = value + d;
            if (stabile) {
                p->p.x = value + d;
            }
            n = glm::dvec2(1,0);
        } else {

            // Quit if no longer valid
            if (p->ep.y >= value + PARTICLE_RAD) {
                return;
            }
            p->ep.y = value + d;
            if (stabile) {
                p->p.y = value + d;
            }
            n = glm::dvec2(0,1);
        }
    } else {
        if (isX) {

            // Quit if no longer valid
            if (p->ep.x <= value - PARTICLE_RAD) {
                return;
            }
            p->ep.x = value - d;
            if (stabile) {
                p->p.x = value - d;
            }
            n = glm::dvec2(-1,0);
        } else {

            // Quit if no longer valid
            if (p->ep.y <= value - PARTICLE_RAD) {
                return;
            }
            p->ep.y = value - d;
            if (stabile) {
                p->p.y = value - d;
            }
            n = glm::dvec2(0,-1);
        }
    }

    if (stabile) {
        return;
    }

    // Apply friction - boundaries have a coefficient of friction of 1
    glm::dvec2 dp = (p->ep - p->p) / (double)counts[idx],
               dpt = dp - glm::dot(dp, n) * n;
    double ldpt = glm::length(dpt);

    if (ldpt < EPSILON) {
        return;
    }

    // Choose between static and kinetic friction
    if (ldpt < sqrt(p->sFriction) * d) {
        p->ep -= dpt;
    } else {
        p->ep -= dpt * min(sqrt(p->kFriction )* d / ldpt, 1.);
    }
}

void BoundaryConstraint::draw(QList<Particle *> *particles)
{
}

double BoundaryConstraint::evaluate(QList<Particle *> *estimates)
{
    Particle *p = estimates->at(idx);
    if (isGreaterThan) {
        if (isX) {
            return (value + PARTICLE_RAD) - p->getP(stabile).x;
        } else {
            return (value + PARTICLE_RAD) - p->getP(stabile).y;
        }
    } else {
        if (isX) {
            return p->getP(stabile).x - (value - PARTICLE_RAD);
        } else {
            return p->getP(stabile).y - (value - PARTICLE_RAD);
        }
    }
}

glm::dvec2 BoundaryConstraint::gradient(QList<Particle *> *estimates, int respect)
{
    if (respect != idx) {
        return glm::dvec2();
    }

    if (isGreaterThan) {
        if (isX) {
            return glm::dvec2(-1,0);
        } else {
            return glm::dvec2(0,-1);
        }
    } else {
        if (isX) {
            return glm::dvec2(1,0);
        } else {
            return glm::dvec2(0,1);
        }
    }
}

void BoundaryConstraint::updateCounts(int *counts)
{
    counts[idx]++;
}
