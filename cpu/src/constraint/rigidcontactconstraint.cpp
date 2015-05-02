#include "rigidcontactconstraint.h"

RigidContactConstraint::RigidContactConstraint(int first, int second, QList<Body *> *bodies, bool st)
    : Constraint(), d(0.0), i1(first), i2(second), stabile(st), bods(bodies)
{

}

RigidContactConstraint::~RigidContactConstraint()
{
}

void RigidContactConstraint::initBoundary(Particle *p1, Particle *p2)
{
    glm::dvec2 x12 = p2->getP(stabile) - p1->getP(stabile);
    double len = glm::length(x12);
    d = len - PARTICLE_DIAM;
    if (d < EPSILON) return;
    x12 = len > EPSILON ? x12 / len : glm::dvec2(0,1);
    double dp = glm::dot(x12, n);
    if (dp < 0) {
        n = x12 - 2.0 * dp * n;
    } else {
        n = x12;
    }
}

void RigidContactConstraint::project(QList<Particle *> *estimates, int *counts)
{
    Particle *p1 = estimates->at(i1), *p2 = estimates->at(i2);
    SDFData dat1 = p1->getSDFData(bods, i1), dat2 = p2->getSDFData(bods, i2);

    if (dat1.distance < 0 || dat2.distance < 0) {
        glm::dvec2 x12 = p2->getP(stabile) - p1->getP(stabile);
        double len = glm::length(x12);
        d = PARTICLE_DIAM - len;
        if (d < EPSILON) return;
        n = -x12 / len;
    } else {
        if (dat1.distance < dat2.distance) {
            d = dat1.distance;
            n = dat1.gradient;
        } else {
            d = dat2.distance;
            n = -dat2.gradient;
        }

        if (d < PARTICLE_DIAM + EPSILON) {
            initBoundary(p1, p2);
        }
    }

    double wSum = p1->tmass + p2->tmass;
    glm::dvec2 dp = -(1.0 / wSum) * d * n,
              dp1 = -p1->tmass * dp  / (double)counts[i1],
              dp2 = p2->tmass * dp / (double)   counts[i2];

    p1->ep += dp1;
    p2->ep += dp2;

    if (stabile) {
        p1->p += dp1;
        p2->p += dp2;
    }

    // Apply friction
    glm::dvec2 nf = glm::normalize(n);
    glm::dvec2 dpf = (p1->ep - p1->p) - (p2->ep - p2->p),
               dpt = dpf - glm::dot(dpf, nf) * nf;
    double ldpt = glm::length(dpt);
    if (ldpt < EPSILON) {
        return;
    }
    double sFric = sqrt(p1->sFriction * p2->sFriction),
            kFric = sqrt(p1->kFriction * p2->kFriction);

    if (ldpt < sFric * d) {
        p1->ep += dpt * p1->tmass / wSum;
        p2->ep -= dpt * p2->tmass / wSum;
    } else {
        glm::dvec2 delta = dpt * min(kFric * d / ldpt, 1.);
        p1->ep += delta * p1->tmass / wSum;
        p2->ep -= delta * p2->tmass / wSum;
    }
}

void RigidContactConstraint::draw(QList<Particle *> *particles)
{

}

double RigidContactConstraint::evaluate(QList<Particle *> *estimates)
{
    Particle *p1 = estimates->at(i1), *p2 = estimates->at(i2);
    SDFData dat1 = p1->getSDFData(bods, i1), dat2 = p2->getSDFData(bods, i2);

    if (dat1.distance < 0 || dat2.distance < 0) {
        glm::dvec2 x12 = p2->getP(stabile) - p1->getP(stabile);
        double len = glm::length(x12);
        d = PARTICLE_DIAM - len;
        n = len > EPSILON ? -x12 / len : glm::dvec2(0,1);
    } else {
        if (dat1.distance < dat2.distance) {
            d = dat1.distance;
            n = dat1.gradient;
        } else {
            d = dat2.distance;
            n = -dat2.gradient;
        }

        if (d < PARTICLE_DIAM + EPSILON) {
            initBoundary(p1, p2);
        }
    }

    return d;
}

glm::dvec2 RigidContactConstraint::gradient(QList<Particle *> *estimates, int respect)
{
    if (respect == i1) {
        return -n;
    }

    if (respect == i2) {
        return n;
    }

    return glm::dvec2();
}

void RigidContactConstraint::updateCounts(int *counts)
{
    counts[i1]++;
    counts[i2]++;
}
