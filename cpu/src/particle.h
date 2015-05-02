#ifndef PARTICLE_H
#define PARTICLE_H

#include "includes.h"

#define PARTICLE_RAD .25
#define PARTICLE_DIAM .5

// Phase of mass for particles
enum Phase {
    SOLID,
    FLUID,
    GAS,
    NUM_PHASES
};

struct Body;
struct SDFData;

// Individual particle representation
struct Particle
{
    glm::dvec2 p, ep, v; // position, guess position, and velocity
    double imass, tmass, sFriction, kFriction; // inverse mass, temporary height-scaled mass, coeffs of friction
    int bod; // body (if any) this particle belongs to, for disabling collisions
    Phase ph; // phase of this particle

    Particle()
        : p(glm::dvec2()), v(glm::dvec2()), ph(NUM_PHASES) { init(0); }

    Particle(glm::dvec2 pos, double mass, Phase phase = SOLID)
        : p(pos), v(glm::dvec2()), ph(phase) { init(mass); }

    Particle(glm::dvec2 pos, glm::dvec2 vel, double mass, Phase phase)
        : p(pos), v(vel), ph(phase) { init(mass); }

    void init(double mass) {
        ep = glm::dvec2();
        bod = -1;

        if (mass <= 0) {
            imass = -mass;
        } else {
            imass = 1. / mass;
        }
        tmass = imass;

        sFriction = 0;
        kFriction = 0; // usually smaller the coefficient of static friction
    }

    inline void setStatic() { imass = 0.; }

    inline glm::dvec2 guess(double seconds) {
        return imass == 0. ? p : p + seconds * v;
    }

    inline void confirmGuess() {
        if (glm::length(ep - p) < EPSILON) {
            v = glm::dvec2(0,0); return;
        }
        p = ep;
    }

    void scaleMass() {
        if (imass != 0.0) {
            tmass = 1. / ((1. / imass) * exp(-p.y));
        } else {
            tmass = 0.0;
        }
    }

    // Used for stabilization-related constraints
    inline glm::dvec2 getP(bool stabile) { return stabile ? p : ep; }

    SDFData getSDFData(QList<Body *> *bodies, int idx);
};

// Signed distance field data for rigid-body collisions
struct SDFData {
    SDFData()
        : gradient(glm::dvec2()), distance(-1.0) {}

    SDFData(const glm::dvec2 grad, double dist)
        : gradient(grad), distance(dist) {}

    inline void rotate(double angle) { gradient = glm::rotate(gradient, angle); }

    glm::dvec2 gradient;
    double distance;
};

// Groups of constraints to be solved together
enum ConstraintGroup {
    STABILIZATION,
    CONTACT,
    STANDARD,
    SHAPE,
    NUM_CONSTRAINT_GROUPS
};

// Abstract superclass of all constraint types
class Constraint
{
public:
    Constraint() : stiffness(1) {}
    virtual ~Constraint() {}

    virtual void draw(QList<Particle *> *particles) = 0;

    // For iterative solving of constraints
    virtual void project(QList<Particle *> *estimates) = 0;

    // For matrix-oriented solving of constraints
    virtual double evaluate(QList<Particle *> *estimates) = 0;
    virtual glm::dvec2 gradient(QList<Particle *> *estimates, int respect) = 0;
    virtual void updateCounts(int *counts) = 0;

protected:
    double stiffness;
};

// A single rigid body
struct Body
{
    virtual ~Body() { delete shape; }
    QList<int> particles; // index into global particles list
    QHash<int, glm::dvec2> rs; // map from global particles index to r vector
    QHash<int, SDFData> sdf; // map from global particles index to SDF data
    Constraint *shape;
    glm::dvec2 center; // center of mass
    double imass, angle; // total inverse mass

    void draw(QList<Particle *> *parts);

    void updateCOM(QList<Particle *> *estimates, bool useEstimates = true);
    void computeRs(QList<Particle *> *estimates);
};

#endif // PARTICLE_H
