#ifndef SIMULATION_H
#define SIMULATION_H

#include "includes.h"
#include "particle.h"
#include "solver.h"

// Number of solver iterations per timestep
#define SOLVER_ITERATIONS 3

// Iterative or matrix solve
#define ITERATIVE

// Use stabilization pass or not, and if so how many iterations
//#define USE_STABILIZATION
#define STABILIZATION_ITERATIONS 2

// Gravity scaling factor for gases
#define ALPHA .3

// Built-in simulation scenes
enum SimulationType {
    FRICTION_TEST,
    GRANULAR_TEST,
    STACKS_TEST,
    WALL_TEST,
    PENDULUM_TEST,
    ROPE_TEST,
    FLUID_TEST,
    FLUID_SOLID_TEST,
    GAS_TEST,
    WATER_BALLOON_TEST,
    CRADLE_TEST,
    NUM_SIMULATION_TYPES
};

// The basic simulation, implementing the "main solve loop" from the paper.
class Simulation
{
public:
    Simulation();
    virtual ~Simulation();
    void init(SimulationType type);

    // Initializers for test scenes
    void initFriction();
    void initGranular();
    void initBoxes();
    void initPendulum();
    void initRope();
    void initFluid();
    void initFluidSolid();
    void initWall();
    void initGas();
    void initWaterBalloon();
    void initNewtonsCradle();

    // Basic interaction events
    void tick(double seconds);
    void draw();
    void resize(const glm::ivec2 &dim);
    void mousePressed(const glm::dvec2 &p);

    // Debug information and flags
    int getNumParticles();
    double getKineticEnergy();
    bool debug;

private:

    // Reset the simulation
    void clear();

    // Creation functions for different types of matter
    Body *createRigidBody(QList<Particle *> *verts, QList<SDFData> *sdfData);
    void createFluid(QList<Particle *> *particles, double density);
    void createGas(QList<Particle *> *particles, double density);

    // Simple drawing routines
    void drawGrid();
    void drawParticles();
    void drawBodies();
    void drawGlobals();
    void drawCircle();

    void setColor(int body, float alpha);

    // Counts for iterative particle solver
    int *m_counts;

    // Storage of global particles, rigid bodies, and general constraints
    QList<Particle *> m_particles;
    QList<Body *> m_bodies;
    QHash<ConstraintGroup, QList<Constraint *> > m_globalConstraints;

    // Solvers for regular and contact constraints
    Solver m_standardSolver;
    Solver m_contactSolver;

    // Drawing and boundary information
    glm::ivec2 m_dimensions;
    glm::dvec2 m_xBoundaries, m_yBoundaries;
    glm::dvec2 m_gravity;
    glm::dvec2 m_point;
};

#endif // SIMULATION_H
