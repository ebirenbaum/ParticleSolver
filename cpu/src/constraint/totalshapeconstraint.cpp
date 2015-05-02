#include "totalshapeconstraint.h"

TotalShapeConstraint::TotalShapeConstraint(Body *bod, double stiff)
    : Constraint(), body(bod)
{
    stiffness = stiff;
}

TotalShapeConstraint::~TotalShapeConstraint()
{

}

void TotalShapeConstraint::project(QList<Particle *> *estimates, int *counts)
{
    body->updateCOM(estimates);

    // implemented using http://labs.byhook.com/2010/06/29/particle-based-rigid-bodies-using-shape-matching/
    for (int i = 0; i < body->particles.size(); i++) {
        int idx = body->particles[i];
        Particle *p = estimates->at(idx);
        p->ep += (guess(idx) - p->ep) * stiffness;
    }
}

void TotalShapeConstraint::draw(QList<Particle *> *particles)
{
    glColor3f(0,1,0);
    glBegin(GL_LINES);

    for (int i = 0; i < body->particles.size(); i++) {
        int idx = body->particles[i];
        Particle *p = particles->at(idx);

        glVertex2f(p->p.x, p->p.y);
        glVertex2f(body->center.x, body->center.y);
    }

    glEnd();

    glPointSize(3);
    glBegin(GL_POINTS);
    glVertex2f(body->center.x, body->center.y);
    for (int i = 0; i < body->particles.size(); i++) {
        int idx = body->particles[i];
        Particle *p = particles->at(idx);
        glVertex2f(p->p.x, p->p.y);
    }
    glEnd();
}

double TotalShapeConstraint::evaluate(QList<Particle *> *estimates)
{
    (void) estimates;
    return 0;
}

glm::dvec2 TotalShapeConstraint::gradient(QList<Particle *> *estimates, int respect)
{
//    if (body->rs.contains(respect)) {
//        Particle *p = estimates->at(respect);
//        glm::dvec2 out = guess(respect) - p->ep;
//        if (out == glm::dvec2()) {
//            return glm::dvec2(0,0);
//        }
//        return -glm::normalize(out);
//    }
    (void) estimates;
    (void) respect;
    return glm::dvec2();
}

void TotalShapeConstraint::updateCounts(int *counts)
{
    for (int i = 0; i < body->particles.size(); i++) {
        counts[body->particles[i]]++;
    }
}

glm::dvec2 TotalShapeConstraint::guess(int idx)
{
    double c = cos(body->angle), s = sin(body->angle);

    glm::dvec2 q = body->rs[idx],
               d = glm::dvec2(c * q.x - s * q.y, s * q.x + c * q.y);
    return d + body->center;
}
