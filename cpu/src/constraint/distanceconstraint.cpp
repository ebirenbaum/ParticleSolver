#include "distanceconstraint.h"

DistanceConstraint::DistanceConstraint(double distance, int first, int second, bool st)
    : Constraint(), d(distance), i1(first), i2(second), stabile(st)
{

}

DistanceConstraint::DistanceConstraint(int first, int second, QList<Particle *> *particles)
    : Constraint(), d(0.0), i1(first), i2(second)
{
    d = glm::length(particles->at(i1)->p - particles->at(i2)->p);
}

DistanceConstraint::~DistanceConstraint()
{

}

void DistanceConstraint::project(QList<Particle *> *estimates)
{
    Particle *p1 = estimates->at(i1), *p2 = estimates->at(i2);

    if (p1->imass == 0.f && p2->imass == 0.f) {
        return;
    }

    glm::dvec2 diff = p1->ep - p2->ep;
    double wSum = p1->imass + p2->imass,
            dist = glm::length(diff),
            mag = dist - d,
            scale = mag / wSum;

    glm::dvec2 dp = (scale / dist) * diff,
              dp1 = -p1->imass * dp,
              dp2 = p2->imass * dp;

    p1->ep += dp1;
    p2->ep += dp2;
}

void DistanceConstraint::draw(QList<Particle *> *particles)
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

double DistanceConstraint::evaluate(QList<Particle *> *estimates)
{
    Particle *p1 = estimates->at(i1), *p2 = estimates->at(i2);
    return glm::length(p1->getP(stabile) - p2->getP(stabile)) - d;
}

glm::dvec2 DistanceConstraint::gradient(QList<Particle *> *estimates, int respect)
{
    if (!(respect == i1 || respect == i2)) {
        return glm::dvec2();
    }

    Particle *p1 = estimates->at(i1), *p2 = estimates->at(i2);
    glm::dvec2 n = glm::normalize(p1->getP(stabile) - p2->getP(stabile));
    if (respect == i1) {
        return n;
    } else {
        return -n;
    }
}

void DistanceConstraint::updateCounts(int *counts)
{
    counts[i1]++;
    counts[i2]++;
}
