#include "contactconstraint.h"

ContactConstraint::ContactConstraint(int first, int second, bool st)
    : Constraint(), i1(first), i2(second), stabile(st)
{
}

ContactConstraint::~ContactConstraint()
{

}

void ContactConstraint::project(QList<Particle *> *estimates, int *counts)
{
    Particle *p1 = estimates->at(i1), *p2 = estimates->at(i2);
    if (p1->imass == 0.f && p2->imass == 0.f) {
        return;
    }

    glm::dvec2 diff = p1->getP(stabile) - p2->getP(stabile);
    double wSum = p1->imass + p2->imass,
            dist = glm::length(diff),
            mag = dist - PARTICLE_DIAM;

    // Previous iterations have moved particles out of collision
    if (mag > 0) {
        return;
    }

    double scale = mag / wSum;
    glm::dvec2 dp = (scale / dist) * diff,
            dp1 = -p1->imass * dp / (double)counts[i1],
              dp2 = p2->imass * dp / (double)counts[i2];

    p1->ep += dp1;
    p2->ep += dp2;

    if (stabile) {
        p1->p += dp1;
        p2->p += dp2;
    }
}

void ContactConstraint::draw(QList<Particle *> *particles)
{
    Particle *p1 = particles->at(i1), *p2 = particles->at(i2);

    glColor3f(1,1,0);
    glBegin(GL_LINES);

    glVertex2f(p1->p.x, p1->p.y);
    glVertex2f(p2->p.x, p2->p.y);

    glEnd();

    glPointSize(3);
    glBegin(GL_POINTS);

    glVertex2f(p1->p.x, p1->p.y);
    glVertex2f(p2->p.x, p2->p.y);

    glEnd();
}

double ContactConstraint::evaluate(QList<Particle *> *estimates)
{
    Particle *p1 = estimates->at(i1), *p2 = estimates->at(i2);
    double dist = glm::length(p1->getP(stabile) - p2->getP(stabile));
    return dist > PARTICLE_DIAM ? 0 : dist - PARTICLE_DIAM;
}

glm::dvec2 ContactConstraint::gradient(QList<Particle *> *estimates, int respect)
{
    if (!(respect == i1 || respect == i2)) {
        return glm::dvec2();
    }

    Particle *p1 = estimates->at(i1), *p2 = estimates->at(i2);
    glm::dvec2 diff = p1->getP(stabile) - p2->getP(stabile);
    double dist = glm::length(diff);

    if (dist > PARTICLE_DIAM) {
        return glm::dvec2();
    }

    glm::dvec2 n = diff / dist;
    if (respect == i1) {
        return n;
    } else {
        return -n;
    }
}

void ContactConstraint::updateCounts(int *counts)
{
    counts[i1]++;
    counts[i2]++;
}

