#include "simulation.h"

#include "distanceconstraint.h"
#include "totalshapeconstraint.h"
#include "boundaryconstraint.h"
#include "contactconstraint.h"
#include "rigidcontactconstraint.h"
#include "totalfluidconstraint.h"
#include "gasconstraint.h"

Simulation::Simulation()
{
    m_counts = NULL;
    init(WRECKING_BALL);
    debug = true;
}

Simulation::~Simulation()
{
    clear();
}

void Simulation::clear() {
    for(int i = m_particles.size()-1; i >= 0; i--) {
        Particle *p = m_particles.at(i);
        m_particles.removeAt(i);
        delete(p);
    }
    for(int i = m_smokeEmitters.size()-1; i >= 0; i--) {
        OpenSmokeEmitter *p = m_smokeEmitters.at(i);
        m_smokeEmitters.removeAt(i);
        delete(p);
    }
    for(int i = m_fluidEmitters.size()-1; i >= 0; i--) {
        FluidEmitter *p = m_fluidEmitters.at(i);
        m_fluidEmitters.removeAt(i);
        delete(p);
    }
    for(int i = m_bodies.size()-1; i >= 0; i--) {
        Body *b = m_bodies.at(i);
        m_bodies.removeAt(i);
        delete(b);
    }
    for (int i = 0; i < NUM_CONSTRAINT_GROUPS; i++) {
        if(m_globalConstraints.contains((ConstraintGroup) i)) {
            QList<Constraint *> group = m_globalConstraints[(ConstraintGroup) i];
            for (int j = group.size()-1; j >=0; j--) {
                Constraint *c = group.at(j);
                for (int k = 0; k < NUM_CONSTRAINT_GROUPS; k++) {
                    if(m_globalConstraints.contains((ConstraintGroup) k)) {
                        m_globalConstraints[(ConstraintGroup) k].removeAll(c);
                    }
                }
                delete(c);
            }
        }
    }

    if (m_counts) {
        delete[] m_counts;
    }
}

void Simulation::init(SimulationType type)
{
    this->clear();

    // Default gravity value
    m_gravity = glm::dvec2(0,-9.8);

    switch (type) {
    case FRICTION_TEST:
        initFriction(); break;
    case SDF_TEST:
        initSdf(); break;
    case GRANULAR_TEST:
        initGranular(); break;
    case STACKS_TEST:
        initBoxes(); break;
    case WALL_TEST:
        initWall(); break;
    case PENDULUM_TEST:
        initPendulum(); break;
    case ROPE_TEST:
        initRope(); break;
    case FLUID_TEST:
        initFluid(); break;
    case FLUID_SOLID_TEST:
        initFluidSolid(); break;
    case GAS_ROPE_TEST:
        initRopeGas(); break;
    case WATER_BALLOON_TEST:
        initWaterBalloon(); break;
    case CRADLE_TEST:
        initNewtonsCradle(); break;
    case SMOKE_OPEN_TEST:
        initSmokeOpen(); break;
    case SMOKE_CLOSED_TEST:
        initSmokeClosed(); break;
    case VOLCANO_TEST:
        initVolcano(); break;
    case WRECKING_BALL:
        initWreckingBall(); break;
    default:
        initBoxes(); break;
    }

    // Set up the M^-1 matrix
    m_standardSolver.setupM(&m_particles);

    m_counts = new int[m_particles.size()];
}

// (#) in the main simulation loop refer to lines from the main loop in the paper
void Simulation::tick(double seconds)
{
    QHash<ConstraintGroup, QList<Constraint *> > constraints;

    // Add all rigid body shape constraints
    for (int i = 0; i < m_bodies.size(); i++) {
        Body *b = m_bodies[i];
        if (TotalShapeConstraint *c = dynamic_cast<TotalShapeConstraint *>(b->shape)) {
            constraints[SHAPE].append(c);
        } else {
            cout << "Rigid body's attached constraint was not a shape constraint." << endl;
            exit(1);
        }
    }

    // Add all other global constraints
    for (int i = 0; i < m_globalConstraints.size(); i++) {
        QList<Constraint *> group = m_globalConstraints[(ConstraintGroup) i];
        for (int j = 0; j < group.size(); j++) {
            constraints[(ConstraintGroup) i].append(group.at(j));
        }
    }

    // (1) For all particles
    for (int i = 0; i < m_particles.size(); i++) {
        Particle *p = m_particles[i];

        // (2) Apply forces
        glm::dvec2 myGravity = m_gravity;
        if(p->ph == GAS) myGravity *= ALPHA;
//        for(OpenSmokeEmitter *e: m_emitters) {
//            for(Particle *p: m_particles) {
//                if(glm::distance(p->p, e->getPosn()) < 1) {
//                    p->f += glm::dvec2(0,.03);
//                }
//            }
//        }
        p->v = p->v + seconds * myGravity + seconds * p->f;
        p->f = glm::dvec2();

        // (3) Predict positions, reset n
        p->ep = p->guess(seconds);
        m_counts[i] = 0;

        // (4) Apply mass scaling (used by certain constraints)
        p->scaleMass();
    }
    // (5) End for

    m_contactSolver.setupM(&m_particles, true);

    // (6) For all particles
    for (int i = 0; i < m_particles.size(); i++) {
        Particle *p = m_particles[i];

        // (7) Find neighboring particles and solid contacts, naive solution
        for (int j = i + 1; j < m_particles.size(); j++) {
            Particle *p2 = m_particles[j];

            // Skip collision between two immovables
            if (p->imass == 0 && p2->imass == 0) {
                continue;

            // Skip collisions betwee particles in the same rigid body
            } else if (p->ph == SOLID && p2->ph == SOLID && p->bod == p2->bod && p->bod != -1) {
                continue;
            } else {

                // Collision happens when circles overlap
                double dist = glm::distance(p->ep, p2->ep);
                if (dist < PARTICLE_DIAM - EPSILON) {

                    // Rigid contact constraints (which include friction) apply to solid-solid contact
                    if (p->ph == SOLID && p2->ph == SOLID) {
                        constraints[CONTACT].append(new RigidContactConstraint(i, j, &m_bodies));
#ifdef USE_STABILIZATION
                        constraints[STABILIZATION].append(new RigidContactConstraint(i, j, &m_bodies, true));
#endif
                    // Regular contact constraints (which have no friction) apply to other solid-other contact
                    } else if (p->ph == SOLID || p2->ph == SOLID) {
                        constraints[CONTACT].append(new ContactConstraint(i, j));
                    }
                }
            }
        }

        // (8) Find solid boundary contacts
        if (p->ep.x < m_xBoundaries.x + PARTICLE_RAD) {
            constraints[CONTACT].append(new BoundaryConstraint(i, m_xBoundaries.x, true, true));
#ifdef USE_STABILIZATION
            constraints[STABILIZATION].append(new BoundaryConstraint(i, m_xBoundaries.x, true, true, true));
#endif
        } else if (p->ep.x > m_xBoundaries.y - PARTICLE_RAD) {
            constraints[CONTACT].append(new BoundaryConstraint(i, m_xBoundaries.y, true, false));
#ifdef USE_STABILIZATION
            constraints[STABILIZATION].append(new BoundaryConstraint(i, m_xBoundaries.y, true, false, true));
#endif
        }

        if (p->ep.y < m_yBoundaries.x + PARTICLE_RAD) {
            constraints[CONTACT].append(new BoundaryConstraint(i, m_yBoundaries.x, false, true));
#ifdef USE_STABILIZATION
            constraints[STABILIZATION].append(new BoundaryConstraint(i, m_yBoundaries.x, false, true, true));
#endif
        } else if (p->ep.y > m_yBoundaries.y - PARTICLE_RAD) {
            constraints[CONTACT].append(new BoundaryConstraint(i, m_yBoundaries.y, false, false));
#ifdef USE_STABILIZATION
            constraints[STABILIZATION].append(new BoundaryConstraint(i, m_yBoundaries.y, false, false, true));
#endif
        }
    }
    // (9) End for

    m_contactSolver.setupSizes(m_particles.size(), &constraints[STABILIZATION]);

#ifdef ITERATIVE

    // (17) For constraint group
    for (int j = 0; j < (int) NUM_CONSTRAINT_GROUPS; j++) {
        ConstraintGroup g = (ConstraintGroup) j;

        // Skip the stabilization constraints
        if (g == STABILIZATION) {
            continue;
        }

        //  (18, 19, 20) Update n based on constraints in g
        for (int k = 0; k < constraints[g].size(); k++) {
            constraints[g].at(k)->updateCounts(m_counts);
        }
    }

#endif

#ifdef USE_STABILIZATION

    // (10) For stabilization iterations
    for (int i = 0; i < STABILIZATION_ITERATIONS; i++) {

#ifdef ITERATIVE
        // (11, 12, 13, 14) Solve contact constraints and update p, ep, and n
        for (int k = 0; k < constraints[STABILIZATION].size(); k++) {
            constraints[STABILIZATION].at(k)->project(&m_particles, m_counts);
        }
#else
        // (11, 12, 13, 14) Solve contact constraints and update p, ep, and n
        if (constraints[STABILIZATION].size() > 0) {
            m_contactSolver.solveAndUpdate(&m_particles, &constraints[STABILIZATION], true);
        } else {
            break;
        }
#endif

    }
    // (15) End for

#endif

#ifdef ITERATIVE

    // (16) For solver iterations
    for (int i = 0; i < SOLVER_ITERATIONS; i++) {

        // (17) For constraint group
        for (int j = 0; j < (int) NUM_CONSTRAINT_GROUPS; j++) {
            ConstraintGroup g = (ConstraintGroup) j;

            // Skip the stabilization constraints
            if (g == STABILIZATION) {
                continue;
            }

            //  (18, 19, 20) Solve constraints in g and update ep
            for (int k = 0; k < constraints[g].size(); k++) {
                constraints[g].at(k)->project(&m_particles, m_counts);
            }
        }
    }

#else

    m_standardSolver.setupSizes(m_particles.size(), &constraints[STANDARD]);
    m_contactSolver.setupSizes(m_particles.size(), &constraints[CONTACT]);

    // (16) For solver iterations
    for (int i = 0; i < SOLVER_ITERATIONS; i++) {

        // (17, 18, 19, 20) for constraint group, solve constraints and update ep
        if (constraints[CONTACT].size() > 0) {
            m_contactSolver.solveAndUpdate(&m_particles, &constraints[CONTACT]);
        }

        if (constraints[STANDARD].size() > 0) {
            m_standardSolver.solveAndUpdate(&m_particles, &constraints[STANDARD]);
        }

        if (constraints[SHAPE].size() > 0) {
            for (int j = 0; j < constraints[SHAPE].size(); j++) {
                constraints[SHAPE][j]->project(&m_particles, m_counts);
            }
        }
        // (21) End for
    }
    // (22) End for
#endif

    // (23) For all particles
    for (int i = 0; i < m_particles.size(); i++) {
        Particle *p = m_particles[i];

        // (24) Update velocities
        p->v = (p->ep - p->p) / seconds;

        // (25, 26) Advect diffuse particles, apply internal forces
        /// TODO

        // (27) Update positions or apply sleeping
        p->confirmGuess();
    }
    // (28) End for

    // Delete temporary conact constraints
    for(int i = constraints[CONTACT].size()-1; i >= 0; i--) {
        Constraint *c = constraints[CONTACT].at(i);
        constraints[CONTACT].removeAt(i);
        delete(c);
    }
    for(int i = constraints[STABILIZATION].size()-1; i >= 0; i--) {
        Constraint *c = constraints[STABILIZATION].at(i);
        constraints[STABILIZATION].removeAt(i);
        delete(c);
    }

    for(OpenSmokeEmitter *e: m_smokeEmitters) {
        e->tick(&m_particles, seconds);
        // (8) Find solid boundary contacts
        for(Particle *p: *(e->getParticles())) {
            if (p->p.x < m_xBoundaries.x) {
                p->p.x = m_xBoundaries.x;
            } else if (p->p.x > m_xBoundaries.y) {
                p->p.x = m_xBoundaries.y;
            }
            if (p->p.y < m_yBoundaries.x) {
                p->p.y = m_yBoundaries.x;
            } else if (p->p.y > m_yBoundaries.y) {
                p->p.y = m_yBoundaries.y;
            }
        }
    }
    for(FluidEmitter *e: m_fluidEmitters) {
        e->tick(&m_particles, seconds);
    }
    delete[] m_counts;
    m_counts = new int[m_particles.size()];
}

Body *Simulation::createRigidBody(QList<Particle *> *verts, QList<SDFData> *sdfData)
{
    if(verts->size() <= 1) {
        cout << "Rigid bodies must be at least 2 points." << endl;
        exit(1);
    }

    // Compute the total mass, add all the particles to the system and the body
    Body *body = new Body(); int offset = m_particles.size(), bodyIdx = m_bodies.size();
    double totalMass = 0.0;
    for (int i = 0; i < verts->size(); i++) {
        Particle *p = verts->at(i);
        p->bod = bodyIdx;
        p->ph = SOLID;

        if (p->imass == 0.0) {
            cout << "A rigid body cannot have a point of infinite mass." << endl;
            exit(1);
        }

        totalMass += (1.0 / p->imass);

        m_particles.append(p);
        body->particles.append(i + offset);
        body->sdf[i + offset] = sdfData->at(i);
    }

    // Update the body's global properties, including initial r_i vectors
    body->imass = 1.0 / totalMass;
    body->updateCOM(&m_particles, false);
    body->computeRs(&m_particles);
    body->shape = new TotalShapeConstraint(body);

    m_bodies.append(body);
    return body;
}

GasConstraint *Simulation::createGas(QList<Particle *> *verts, double density, bool open=false)
{
    int offset = m_particles.size();
    int bod = 100 * frand();
    QList<int> indices;
    for (int i = 0; i < verts->size(); i++) {
        Particle *p = verts->at(i);
        p->ph = GAS;
        p->bod = bod;

        if (p->imass == 0.0) {
            cout << "A fluid cannot have a point of infinite mass." << endl;
            exit(1);
        }

        m_particles.append(p);
        indices.append(offset + i);
    }
    GasConstraint *gs = new GasConstraint(density, &indices, open);
    m_globalConstraints[STANDARD].append(gs);
    return gs;
}

TotalFluidConstraint *Simulation::createFluid(QList<Particle *> *verts, double density)
{
    int offset = m_particles.size();
    int bod = 100 * frand();
    QList<int> indices;
    for (int i = 0; i < verts->size(); i++) {
        Particle *p = verts->at(i);
        p->ph = FLUID;
        p->bod = bod;

        if (p->imass == 0.0) {
            cout << "A fluid cannot have a point of infinite mass." << endl;
            exit(1);
        }

        m_particles.append(p);
        indices.append(offset + i);
    }
    TotalFluidConstraint *fs = new TotalFluidConstraint(density, &indices);
    m_globalConstraints[STANDARD].append(fs);
    return fs;
}

void Simulation::createSmokeEmitter(glm::dvec2 posn, double particlesPerSec, GasConstraint *gs)
{
    m_smokeEmitters.append(new OpenSmokeEmitter(posn, particlesPerSec, gs));
}

void Simulation::createFluidEmitter(glm::dvec2 posn, double particlesPerSec, TotalFluidConstraint *fs) {
    m_fluidEmitters.append(new FluidEmitter(posn, particlesPerSec, fs));
}

void Simulation::draw()
{
    drawGrid();
    if (debug) {
        drawParticles();
    }
    drawBodies();
    drawGlobals();
    drawSmoke();

    glColor3f(1,1,1);
    glPointSize(5);
    glBegin(GL_POINTS);
    glVertex2f(m_point.x, m_point.y);
    glEnd();
}

void Simulation::resize(const glm::ivec2 &dim)
{
    m_dimensions = dim;
}

void Simulation::drawGrid()
{
    glColor3f(.2,.2,.2);
    glBegin(GL_LINES);

    for (int x = -m_dimensions.x; x <= m_dimensions.x; x++) {
        glVertex2f(x, -m_dimensions.y);
        glVertex2f(x, m_dimensions.y);
    }
    for (int y = -m_dimensions.y; y <= m_dimensions.y; y++) {
        glVertex2f(-m_dimensions.y, y);
        glVertex2f(m_dimensions.y, y);
    }

    glColor3f(1,1,1);

    glVertex2f(-m_dimensions.x, 0);
    glVertex2f(m_dimensions.x, 0);
    glVertex2f(0, -m_dimensions.y);
    glVertex2f(0, m_dimensions.y);

    glEnd();

    glLineWidth(3);
    glBegin(GL_LINES);
    glVertex2f(m_xBoundaries.x, m_yBoundaries.x);
    glVertex2f(m_xBoundaries.x, m_yBoundaries.y);

    glVertex2f(m_xBoundaries.y, m_yBoundaries.x);
    glVertex2f(m_xBoundaries.y, m_yBoundaries.y);

    glVertex2f(m_xBoundaries.x, m_yBoundaries.x);
    glVertex2f(m_xBoundaries.y, m_yBoundaries.x);

    glVertex2f(m_xBoundaries.x, m_yBoundaries.y);
    glVertex2f(m_xBoundaries.y, m_yBoundaries.y);
    glEnd();
    glLineWidth(1);
}

void Simulation::drawParticles()
{
    for (int i = 0; i < m_particles.size(); i++) {
        const Particle *p = m_particles[i];

        if (p->imass == 0.f) {
            glColor3f(1,0,0);
        } else if (p->ph == FLUID || p->ph == GAS){
            glColor3f(0,p->bod / 100., 1-p->bod / 100.);
        } else if (p->ph == SOLID) {
            setColor(p->bod, 1);
        } else {
            glColor3f(0,0,1);
        }

        glPushMatrix();
        glTranslatef(p->p.x, p->p.y, 0);
        glScalef(PARTICLE_RAD, PARTICLE_RAD, 0);
        drawCircle();
        glPopMatrix();
    }

    glEnd();
}

void Simulation::drawBodies()
{
    for (int i = 0; i < m_bodies.size(); i++) {
        Body *b = m_bodies[i];
        if (debug) {
            b->shape->draw(&m_particles);
        } else {
            for (int i = 0; i < b->particles.size(); i++) {
                Particle *p = m_particles[(b->particles[i])];

                glPushMatrix();
                glTranslatef(p->p.x, p->p.y, 0);
                glPushMatrix();
                glRotatef(R2D(b->angle), 0, 0, 1);

                glEnable(GL_BLEND);
                setColor(p->bod, .6);
                glBegin(GL_QUADS);
                glVertex2f(-PARTICLE_RAD, -PARTICLE_RAD);
                glVertex2f(-PARTICLE_RAD, PARTICLE_RAD);
                glVertex2f(PARTICLE_RAD, PARTICLE_RAD);
                glVertex2f(PARTICLE_RAD, -PARTICLE_RAD);
                glEnd();
                glDisable(GL_BLEND);

                glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
                glBegin(GL_QUADS);
                glVertex2f(-PARTICLE_RAD, -PARTICLE_RAD);
                glVertex2f(-PARTICLE_RAD, PARTICLE_RAD);
                glVertex2f(PARTICLE_RAD, PARTICLE_RAD);
                glVertex2f(PARTICLE_RAD, -PARTICLE_RAD);
                glEnd();
                glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

//                glColor3f(0,0,0);
//                glPopMatrix();
//                glm::dvec2 s = b->sdf[b->particles[i]].gradient * b->sdf[b->particles[i]].distance;
//                s = glm::rotate(s, b->angle);
//                glBegin(GL_LINES);
//                glVertex2f(0,0);
//                glVertex2f(s.x, s.y);
//                glEnd();

                glPopMatrix();
            }
        }
    }
}

void Simulation::drawGlobals()
{
    for (int i = 0; i < m_globalConstraints.size(); i++) {
        for (int j = 0; j < m_globalConstraints[(ConstraintGroup) i].size(); j++) {
            m_globalConstraints[(ConstraintGroup)i ][j]->draw(&m_particles);
        }
    }
}

void Simulation::drawSmoke()
{
    glColor3f(1,1,1);
    glBegin(GL_QUADS);
    double rad = PARTICLE_RAD/7.;
    for (int i = 0; i < m_smokeEmitters.size(); i++) {
        QList<Particle *> *particles = m_smokeEmitters.at(i)->getParticles();
        for(int j = 0; j < particles->size(); j++) {
            Particle *p = particles->at(j);
            glVertex2d(p->p.x-rad, p->p.y-rad);
            glVertex2d(p->p.x+rad, p->p.y-rad);
            glVertex2d(p->p.x+rad, p->p.y+rad);
            glVertex2d(p->p.x-rad, p->p.y+rad);
//            glPushMatrix();
//            glTranslatef(p->p.x, p->p.y, 0);
//            glScalef(PARTICLE_RAD/7., PARTICLE_RAD/7., 0);
//            drawCircle();
//            glPopMatrix();
        }
    }
}

void Simulation::setColor(int body, float alpha)
{
    int choice = abs(body) % 5;
    if (choice == 0) {
        glColor4f(1,.7,0,alpha);
    } else if (choice == 1) {
        glColor4f(.35,.75,.95,alpha);
    } else if (choice == 2) {
        glColor4f(1,.55,.8,alpha);
    } else if (choice == 3) {
        glColor4f(.1,.85,.9,alpha);
    } else {
        glColor4f(.1,.9,.6,alpha);
    }
}

void Simulation::drawCircle()
{
    glBegin(GL_TRIANGLE_FAN);

    glVertex2f(0,0);
    for (int f = 0; f <= 32; f++) {
        double a = f * M_PI / 16.f;
        glVertex2f(sin(a), cos(a));
    }

    glEnd();
}

void Simulation::initFriction()
{
    m_xBoundaries = glm::dvec2(-20,20);
    m_yBoundaries = glm::dvec2(0,1000000);

    double root2 = sqrt(2);
    QList<Particle *> vertices;
    QList<SDFData> data;
    data.append(SDFData(glm::normalize(glm::dvec2(-1,-1)), PARTICLE_RAD * root2));
    data.append(SDFData(glm::normalize(glm::dvec2(-1,1)), PARTICLE_RAD * root2));
    data.append(SDFData(glm::normalize(glm::dvec2(0,-1)), PARTICLE_RAD));
    data.append(SDFData(glm::normalize(glm::dvec2(0,1)), PARTICLE_RAD));
    data.append(SDFData(glm::normalize(glm::dvec2(1,-1)), PARTICLE_RAD * root2));
    data.append(SDFData(glm::normalize(glm::dvec2(1,1)), PARTICLE_RAD * root2));

    glm::ivec2 dim = glm::ivec2(3,2);
    for (int x = 0; x < dim.x; x++) {
        double xVal = PARTICLE_DIAM * ((x % dim.x) - dim.x / 2);
        for (int y = 0; y < dim.y; y++) {
            double yVal = (dim.y + (y % dim.y) + 1) * PARTICLE_DIAM;
            Particle *part =new Particle(glm::dvec2(xVal, yVal), (x == 0 && y == 0 ? 1 : 1.));
            part->v.x = 5;
            part->kFriction = .01;
            part->sFriction = .1;
            vertices.append(part);
        }
    }
    Body *body = createRigidBody(&vertices, &data);
}

void Simulation::initGranular()
{
    m_xBoundaries = glm::dvec2(-100,100);
    m_yBoundaries = glm::dvec2(-5, 1000);
    m_gravity = glm::dvec2(0,-9.8);

    for (int i = -15; i <= 15; i++) {
        for (int j = 0; j < 30; j++) {
            glm::dvec2 pos = glm::dvec2(i * (PARTICLE_DIAM + EPSILON), pow(j,1) * (PARTICLE_DIAM) + PARTICLE_RAD + m_yBoundaries.x);
            Particle *part= new Particle(pos, 1, SOLID);
            part->sFriction = .35;
            part->kFriction = .3;
            m_particles.append(part);
        }
    }

    Particle *jerk = new Particle(glm::dvec2(-25.55, 40), 100.f, SOLID);
    jerk->v.x = 8.5;
    m_particles.append(jerk);
}

void Simulation::initSdf()
{
    m_xBoundaries = glm::dvec2(-20,20);
    m_yBoundaries = glm::dvec2(0,1000000);

    int numBoxes = 2;
    double root2 = sqrt(2);
    QList<Particle *> vertices;
    QList<SDFData> data;
    data.append(SDFData(glm::normalize(glm::dvec2(-1,-1)), PARTICLE_RAD * root2));
    data.append(SDFData(glm::normalize(glm::dvec2(-1,0)), PARTICLE_RAD));
    data.append(SDFData(glm::normalize(glm::dvec2(-1,1)), PARTICLE_RAD * root2));
    data.append(SDFData(glm::normalize(glm::dvec2(1,-1)), PARTICLE_RAD * root2));
    data.append(SDFData(glm::normalize(glm::dvec2(1,0)), PARTICLE_RAD));
    data.append(SDFData(glm::normalize(glm::dvec2(1,1)), PARTICLE_RAD * root2));

    glm::ivec2 dim = glm::ivec2(2,3);
    for (int i = numBoxes - 1; i >= 0; i--) {
        for (int x = 0; x < dim.x; x++) {
            double xVal = PARTICLE_DIAM * ((x % dim.x) - dim.x / 2) + i * PARTICLE_RAD;
            for (int y = 0; y < dim.y; y++) {
                double yVal = ((40 * i) * dim.y + (y % dim.y) + 1) * PARTICLE_DIAM;
                Particle *part = new Particle(glm::dvec2(xVal, yVal), 4.);
                if (i > 0) part->v.y = -120;
                vertices.append(part);
            }
        }
        Body *body = createRigidBody(&vertices, &data);
        vertices.clear();
    }
}

void Simulation::initBoxes()
{
    m_xBoundaries = glm::dvec2(-20,20);
    m_yBoundaries = glm::dvec2(0,1000000);

    int numBoxes = 7, numColumns = 2;
    double root2 = sqrt(2);
    QList<Particle *> vertices;
    QList<SDFData> data;
    data.append(SDFData(glm::normalize(glm::dvec2(-1,-1)), PARTICLE_RAD * root2));
    data.append(SDFData(glm::normalize(glm::dvec2(-1,1)), PARTICLE_RAD * root2));
    data.append(SDFData(glm::normalize(glm::dvec2(0,-1)), PARTICLE_RAD));
    data.append(SDFData(glm::normalize(glm::dvec2(0,1)), PARTICLE_RAD));
    data.append(SDFData(glm::normalize(glm::dvec2(1,-1)), PARTICLE_RAD * root2));
    data.append(SDFData(glm::normalize(glm::dvec2(1,1)), PARTICLE_RAD * root2));

    for (int j = -numColumns; j <= numColumns; j++) {
        glm::ivec2 dim = glm::ivec2(3,2);
        for (int i = numBoxes - 1; i >= 0; i--) {
            for (int x = 0; x < dim.x; x++) {
                double xVal = j * 4 + PARTICLE_DIAM * ((x % dim.x) - dim.x / 2);
                for (int y = 0; y < dim.y; y++) {
                    double yVal = ((2 * i + 1) * dim.y + (y % dim.y) + 1) * PARTICLE_DIAM;
                    Particle *part = new Particle(glm::dvec2(xVal, yVal), 4.);
                    part->sFriction = 1.;
                    part->kFriction = 1.;
                    vertices.append(part);
                }
            }
            Body *body = createRigidBody(&vertices, &data);
            vertices.clear();
        }
    }
}

void Simulation::initWall()
{
    m_xBoundaries = glm::dvec2(-50,50);
    m_yBoundaries = glm::dvec2(0,1000000);

    glm::dvec2 dim = glm::dvec2(6,2);
    int height = 9, width = 5;
    double root2 = sqrt(2);
    QList<Particle *> vertices;
    QList<SDFData> data;
    data.append(SDFData(glm::normalize(glm::dvec2(-1,-1)), PARTICLE_RAD * root2));
    data.append(SDFData(glm::normalize(glm::dvec2(-1,1)), PARTICLE_RAD * root2));

    for (int i = 0; i < dim.x - 2; i++) {
        data.append(SDFData(glm::normalize(glm::dvec2(0,-1)), PARTICLE_RAD));
        data.append(SDFData(glm::normalize(glm::dvec2(0,1)), PARTICLE_RAD));
    }

    data.append(SDFData(glm::normalize(glm::dvec2(1,-1)), PARTICLE_RAD * root2));
    data.append(SDFData(glm::normalize(glm::dvec2(1,1)), PARTICLE_RAD * root2));

    for (int j = -width; j <= width; j++) {
        for (int i = height - 1; i >= 0; i--) {
            for (int x = 0; x < dim.x; x++) {
                double num = (i % 2 == 0 ? 3 : -1);
                double xVal = j * (EPSILON + dim.x / 2.) + PARTICLE_DIAM * (x % (int)dim.x) - num * PARTICLE_RAD;
                for (int y = 0; y < dim.y; y++) {
                    double yVal = (i * dim.y + (y % (int)dim.y) + EPSILON) * PARTICLE_DIAM + PARTICLE_RAD;
                    Particle *part = new Particle(glm::dvec2(xVal, yVal), 1.);
                    part->sFriction = 1;
                    part->kFriction = 0;
                    vertices.append(part);
                }
            }
            Body *body = createRigidBody(&vertices, &data);
            vertices.clear();
        }
    }
}

void Simulation::initPendulum()
{
    m_xBoundaries = glm::dvec2(-10,10);
    m_yBoundaries = glm::dvec2(0,1000000);

    int chainLength = 3;
    m_particles.append(new Particle(glm::dvec2(0, chainLength * 3 + 6) * PARTICLE_DIAM + glm::dvec2(0,2), 0, SOLID));

    QList<SDFData> data;
    data.append(SDFData(glm::normalize(glm::dvec2(-1,-1)), PARTICLE_RAD));
    data.append(SDFData(glm::normalize(glm::dvec2(-1,1)), PARTICLE_RAD));
    data.append(SDFData(glm::normalize(glm::dvec2(0,-1)), PARTICLE_RAD));
    data.append(SDFData(glm::normalize(glm::dvec2(0,1)), PARTICLE_RAD));
    data.append(SDFData(glm::normalize(glm::dvec2(1,-1)), PARTICLE_RAD));
    data.append(SDFData(glm::normalize(glm::dvec2(1,1)), PARTICLE_RAD));

    QList<Particle *> vertices;
    double xs[6] = {-1,-1,0,0,1,1};

    for (int i = chainLength; i >= 0; i--) {
        for (int j = 0; j < 6; j++) {
            double y = ((i + 1) * 3 + (j % 2)) * PARTICLE_DIAM + 2;
            vertices.append(new Particle(glm::dvec2(xs[j] * PARTICLE_DIAM, y), 1.));
        }
        Body *body = createRigidBody(&vertices, &data);
        vertices.clear();

        if (i < chainLength) {
            int basePrev = 1 + (chainLength - i - 1) * 6, baseCur = basePrev + 6;
            m_globalConstraints[STANDARD].append(new DistanceConstraint(baseCur + 1, basePrev, &m_particles));
            m_globalConstraints[STANDARD].append(new DistanceConstraint(baseCur + 5, basePrev + 4, &m_particles));
        }
    }

    m_globalConstraints[STANDARD].append(new DistanceConstraint(0, 4, &m_particles));
}

void Simulation::initRope()
{
    double scale = 5.;
    m_xBoundaries = glm::dvec2(-scale,scale);
    m_yBoundaries = glm::dvec2(0,1000000);

    double top = 6, dist = PARTICLE_RAD;

    Particle *e1 = new Particle(glm::dvec2(m_xBoundaries.x, top), 0, SOLID);
    e1->bod = -2;
    m_particles.append(e1);

    for (double i = m_xBoundaries.x + dist; i < m_xBoundaries.y - dist; i += dist) {
        Particle *part = new Particle(glm::dvec2(i, top), 1., SOLID);
        part->bod = -2;
        m_particles.append(part);
        m_globalConstraints[STANDARD].append(
                    new DistanceConstraint(dist, m_particles.size() - 2, m_particles.size() - 1));
    }

    Particle *e2 = new Particle(glm::dvec2(m_xBoundaries.y, top), 0, SOLID);
    e2->bod = -2;
    m_particles.append(e2);

    m_globalConstraints[STANDARD].append(
                new DistanceConstraint(dist, m_particles.size() - 2, m_particles.size() - 1));

    double delta = .7;
    QList<Particle *> particles;

    for(double x = -scale; x < scale; x += delta) {
        for(double y = 10; y < 10 + scale; y += delta) {
            particles.append(new Particle(glm::dvec2(x,y) + .2 * glm::dvec2(frand() - .5, frand() - .5), 1));
        }
    }
    createGas(&particles, 1.75);
}

void Simulation::initFluid()
{
    double scale = 4., delta = .7;
    m_gravity = glm::dvec2(0,-9.8);
    m_xBoundaries = glm::dvec2(-2 * scale,2 * scale);
    m_yBoundaries = glm::dvec2(-2 * scale, 10 * scale);
    QList<Particle *> particles;

    double num = 2.;
    for (int d = 0; d < num; d++) {
        double start = -2 * scale + 4 * scale * (d / num);
        for(double x = start; x < start + (4 * scale / num); x += delta) {
            for(double y = -2 * scale; y < scale; y += delta) {
                particles.append(new Particle(glm::dvec2(x,y) + .2 * glm::dvec2(frand() - .5, frand() - .5), 1));
            }
        }
        createFluid(&particles, 1 + .75 * d);
        particles.clear();
    }
}

void Simulation::initFluidSolid()
{
    double scale = 3., delta = .7;
    m_gravity = glm::dvec2(0, -9.8);
    m_xBoundaries = glm::dvec2(-2 * scale,2 * scale);
    m_yBoundaries = glm::dvec2(-2 * scale, 100 * scale);
    QList<Particle *> particles;

    double num = 1.;
    for (int d = 0; d < num; d++) {
        double start = -2 * scale + 4 * scale * (d / num);
        for(double x = start; x < start + (4 * scale / num); x += delta) {
            for(double y = -2 * scale; y < 2 * scale; y += delta) {
                particles.append(new Particle(glm::dvec2(x,y + 3) + .2 * glm::dvec2(frand() - .5, frand() - .5), 1));
            }
        }
        createFluid(&particles, 1. + 1.25 * (d + 1));
        particles.clear();
    }

    if(true) {
        particles.clear();
        QList<SDFData> data;
        double root2 = sqrt(2);
        glm::ivec2 dim = glm::ivec2(5,2);
        data.append(SDFData(glm::normalize(glm::dvec2(-1,-1)), PARTICLE_RAD * root2));
        data.append(SDFData(glm::normalize(glm::dvec2(-1,1)), PARTICLE_RAD * root2));
        for (int i = 0; i < dim.x - 2; i++) {
            data.append(SDFData(glm::normalize(glm::dvec2(0,-1)), PARTICLE_RAD));
            data.append(SDFData(glm::normalize(glm::dvec2(0,1)), PARTICLE_RAD));
        }
        data.append(SDFData(glm::normalize(glm::dvec2(1,-1)), PARTICLE_RAD * root2));
        data.append(SDFData(glm::normalize(glm::dvec2(1,1)), PARTICLE_RAD * root2));

        for (int x = 0; x < dim.x; x++) {
            double xVal = PARTICLE_DIAM * ((x % dim.x) - dim.x / 2);
            for (int y = 0; y < dim.y; y++) {
                double yVal = (dim.y + (y % dim.y) + 1) * PARTICLE_DIAM;
                particles.append(new Particle(glm::dvec2(xVal-3, yVal + 10), 2));
            }
        }
        Body *body = createRigidBody(&particles, &data);
    }

    if(true) {
        particles.clear();
        QList<SDFData> data;
        double root2 = sqrt(2);
        glm::ivec2 dim = glm::ivec2(5,2);
        data.append(SDFData(glm::normalize(glm::dvec2(-1,-1)), PARTICLE_RAD * root2));
        data.append(SDFData(glm::normalize(glm::dvec2(-1,1)), PARTICLE_RAD * root2));
        for (int i = 0; i < dim.x - 2; i++) {
            data.append(SDFData(glm::normalize(glm::dvec2(0,-1)), PARTICLE_RAD));
            data.append(SDFData(glm::normalize(glm::dvec2(0,1)), PARTICLE_RAD));
        }
        data.append(SDFData(glm::normalize(glm::dvec2(1,-1)), PARTICLE_RAD * root2));
        data.append(SDFData(glm::normalize(glm::dvec2(1,1)), PARTICLE_RAD * root2));

        for (int x = 0; x < dim.x; x++) {
            double xVal = PARTICLE_DIAM * ((x % dim.x) - dim.x / 2);
            for (int y = 0; y < dim.y; y++) {
                double yVal = (dim.y + (y % dim.y) + 1) * PARTICLE_DIAM;
                particles.append(new Particle(glm::dvec2(xVal+3, yVal + 10), .2));
            }
        }
        Body *body = createRigidBody(&particles, &data);
    }
}


void Simulation::initGas()
{
    double scale = 2., delta = .7;
    m_gravity = glm::dvec2(0, -9.8);
    m_xBoundaries = glm::dvec2(-2  * scale,2 * scale);
    m_yBoundaries = glm::dvec2(-2  * scale, 10 * scale);
    QList<Particle *> particles;

    double num = 2.;
    for (int d = 0; d < num; d++) {
        double start = -2 * scale + 4 * scale * (d / num);
        for(double x = start; x < start + (4 * scale / num); x += delta) {
            for(double y = -2 * scale; y < 2 * scale; y += delta) {
                particles.append(new Particle(glm::dvec2(x,y) + .2 * glm::dvec2(frand() - .5, frand() - .5), 1));
            }
        }
        createGas(&particles, .75 + 3*(d));
        particles.clear();
    }

    scale = 3;
    for (int d = 0; d < num; d++) {
        double start = -2 * scale + 4 * scale * (d / num);
        for(double x = start; x < start + (4 * scale / num); x += delta) {
            for(double y = -2 * scale; y < 2 * scale; y += delta) {
                particles.append(new Particle(glm::dvec2(x,y+10) + .2 * glm::dvec2(frand() - .5, frand() - .5), 1));
            }
        }
        createFluid(&particles, 4. + .75 * (d + 1));
        particles.clear();
    }
}

void Simulation::initWaterBalloon()
{
    double scale = 10.;
    m_xBoundaries = glm::dvec2(-scale,scale);
    m_yBoundaries = glm::dvec2(-10,1000000);

    double samples = 60, da = 360. / samples;

    for (int i = 0; i < samples; i++) {
        double angle = D2R(i * da);
        Particle *part = new Particle(glm::dvec2(sin(angle), cos(angle)) * 3., 1);
        part->bod = -2;
        int idx = m_particles.size();
        m_particles.append(part);

        if (i > 0) {
            m_globalConstraints[STANDARD].append(
                        new DistanceConstraint(idx, idx - 1, &m_particles));
        }
    }
    m_globalConstraints[STANDARD].append(
                new DistanceConstraint(0, m_particles.size() - 1, &m_particles));
    int idk = m_particles.size();

    for (int i = 0; i < samples; i++) {
        double angle = D2R(i * da);
        Particle *part = new Particle(glm::dvec2(sin(angle), cos(angle) + 3) * 3., 1);
        part->bod = -3;
        int idx = m_particles.size();
        m_particles.append(part);

        if (i > 0) {
            m_globalConstraints[STANDARD].append(
                        new DistanceConstraint(idx, idx - 1, &m_particles));
        }
    }
    m_globalConstraints[STANDARD].append(
                new DistanceConstraint(idk, m_particles.size() - 1, &m_particles));

    double delta = 1.5 * PARTICLE_RAD;
    QList<Particle *> particles;

    for(double x = -2; x <= 2; x += delta) {
        for(double y = -2; y <= 2; y += delta) {
            particles.append(new Particle(glm::dvec2(x,y) + .2 * glm::dvec2(frand() - .5, frand() - .5), 1));
        }
    }
    createFluid(&particles, 1.75);

    particles.clear();
    for(double x = -2; x <= 2; x += delta) {
        for(double y = -2; y <= 2; y += delta) {
            particles.append(new Particle(glm::dvec2(x,y + 9) + .2 * glm::dvec2(frand() - .5, frand() - .5), 1));
        }
    }
    createFluid(&particles, 1.75);
}

void Simulation::initNewtonsCradle()
{
    m_xBoundaries = glm::dvec2(-10,10);
    m_yBoundaries = glm::dvec2(-5,1000000);

    int n = 2;

    for (int i = -n; i <= n; i++) {
        int idx = m_particles.size();
        m_particles.append(new Particle(glm::dvec2(i * PARTICLE_DIAM, 0), 0.f));
        if (i != -n) {
            m_particles.append(new Particle(glm::dvec2(i * PARTICLE_DIAM, -3), 1.f));
        } else {
            Particle *part = new Particle(glm::dvec2(i * PARTICLE_DIAM - 3, 0), 1.f);
            m_particles.append(part);
        }
        m_globalConstraints[STANDARD].append(new DistanceConstraint(idx, idx+1, &m_particles));
    }
}

void Simulation::initSmokeOpen()
{
    double scale = 2., delta = .63;
    m_gravity = glm::dvec2(0, -9.8);
    m_xBoundaries = glm::dvec2(-3  * scale,3 * scale);
    m_yBoundaries = glm::dvec2(-2  * scale,100 * scale);
    QList<Particle *> particles;

    double start = -2 * scale;
    for(double x = start; x < start + (4 * scale); x += delta) {
        for(double y = -2 * scale; y < 2 * scale; y += delta) {
            particles.append(new Particle(glm::dvec2(x,y) + .2 * glm::dvec2(frand() - .5, frand() - .5), 1));
        }
    }
    GasConstraint *gs = createGas(&particles, 1.5, true);
    particles.clear();

    createSmokeEmitter(glm::dvec2(0,-2*scale+1), 15, gs);
}

void Simulation::initSmokeClosed()
{
    double scale = 2., delta = .63;
    m_gravity = glm::dvec2(0, -9.8);
    m_xBoundaries = glm::dvec2(-2  * scale,2 * scale);
    m_yBoundaries = glm::dvec2(-2  * scale,2 * scale);
    QList<Particle *> particles;

    double start = -2 * scale;
    for(double x = start; x < start + (4 * scale); x += delta) {
        for(double y = -2 * scale; y < 2 * scale; y += delta) {
            particles.append(new Particle(glm::dvec2(x,y) + .2 * glm::dvec2(frand() - .5, frand() - .5), 1));
        }
    }
    GasConstraint *gs = createGas(&particles, 1.5, false);
    particles.clear();

    createSmokeEmitter(glm::dvec2(0,-2*scale+1), 15, NULL);
}

void Simulation::initRopeGas()
{
    double scale = 2., delta = .63;
    m_gravity = glm::dvec2(0, -9.8);
    m_xBoundaries = glm::dvec2(-4  * scale,4 * scale);
    m_yBoundaries = glm::dvec2(-2  * scale,100 * scale);

    double top = 6, dist = PARTICLE_RAD;

    Particle *e1 = new Particle(glm::dvec2(-2.*scale, top), 0, SOLID);
    e1->bod = -2;
    m_particles.append(e1);

    for (double i = -2*scale + dist; i < 2*scale - dist; i += dist) {
        Particle *part = new Particle(glm::dvec2(i, top), .2, SOLID);
        part->bod = -2;
        m_particles.append(part);
        m_globalConstraints[STANDARD].append(
                    new DistanceConstraint(dist, m_particles.size() - 2, m_particles.size() - 1));
    }

    Particle *e2 = new Particle(glm::dvec2(2*scale, top), 0, SOLID);
    e2->bod = -2;
    m_particles.append(e2);

    m_globalConstraints[STANDARD].append(
                new DistanceConstraint(dist, m_particles.size() - 2, m_particles.size() - 1));

    QList<Particle *> particles;

    double start = -.5 * scale;
    for(double x = start; x < start + (1 * scale); x += delta) {
        for(double y = -.5 * scale; y < .5 * scale; y += delta) {
            particles.append(new Particle(glm::dvec2(x,y) + .2 * glm::dvec2(frand() - .5, frand() - .5), 1));
        }
    }
    GasConstraint *gs = createGas(&particles, 1.5, true);

    createSmokeEmitter(glm::dvec2(0,0), 15, gs);
    particles.clear();
}

void Simulation::initVolcano()
{
    double scale = 10., delta = .2;

    for(double x = 1.; x <= scale; x+=delta) {
        m_particles.append(new Particle(glm::dvec2(-x,scale-x), 0));
        m_particles.append(new Particle(glm::dvec2(x,scale-x), 0));
    }

    m_gravity = glm::dvec2(0,-9.8);
    m_xBoundaries = glm::dvec2(-2 * scale,2 * scale);
    m_yBoundaries = glm::dvec2(0, 10 * scale);
    QList<Particle *> particles;

    delta = .8;
    for(double y = 0.; y < scale-1.; y+=delta) {
        for(double x = 0.; x < scale-y-1; x += delta) {
            particles.append(new Particle(glm::dvec2(x,y) + .2 * glm::dvec2(frand() - .5, frand() - .5), 1));
            particles.append(new Particle(glm::dvec2(-x,y) + .2 * glm::dvec2(frand() - .5, frand() - .5), 1));
        }
    }
    TotalFluidConstraint *fs = createFluid(&particles, 1);
    particles.clear();

    createFluidEmitter(glm::dvec2(0,0), scale*4, fs);

//    double top = scale-.5, dist = PARTICLE_RAD;

//    Particle *e1 = new Particle(glm::dvec2(-1-dist, top), 0, SOLID);
//    e1->bod = -2;
//    m_particles.append(e1);

//    for (double i = -1; i <= 2; i += dist) {
//        Particle *part = new Particle(glm::dvec2(i, top), 1, SOLID);
//        part->bod = -2;
//        m_particles.append(part);
//        m_globalConstraints[STANDARD].append(
//                    new DistanceConstraint(dist, m_particles.size() - 2, m_particles.size() - 1));
//    }
}

void Simulation::initWreckingBall()
{
    m_xBoundaries = glm::dvec2(-15,100);
    m_yBoundaries = glm::dvec2(0,1000000);

    glm::dvec2 dim = glm::dvec2(6,2);
    int height = 8, width = 2;
    double root2 = sqrt(2);
    QList<Particle *> vertices;
    QList<SDFData> data;
    data.append(SDFData(glm::normalize(glm::dvec2(-1,-1)), PARTICLE_RAD * root2));
    data.append(SDFData(glm::normalize(glm::dvec2(-1,1)), PARTICLE_RAD * root2));

    for (int i = 0; i < dim.x - 2; i++) {
        data.append(SDFData(glm::normalize(glm::dvec2(0,-1)), PARTICLE_RAD));
        data.append(SDFData(glm::normalize(glm::dvec2(0,1)), PARTICLE_RAD));
    }

    data.append(SDFData(glm::normalize(glm::dvec2(1,-1)), PARTICLE_RAD * root2));
    data.append(SDFData(glm::normalize(glm::dvec2(1,1)), PARTICLE_RAD * root2));

    for (int j = -width; j <= width; j++) {
        for (int i = height - 1; i >= 0; i--) {
            for (int x = 0; x < dim.x; x++) {
                double num = (i % 2 == 0 ? 3 : -1);
                double xVal = j * (EPSILON + dim.x / 2.) + PARTICLE_DIAM * (x % (int)dim.x) - num * PARTICLE_RAD;
                for (int y = 0; y < dim.y; y++) {
                    double yVal = (i * dim.y + (y % (int)dim.y) + EPSILON) * PARTICLE_DIAM + PARTICLE_RAD;
                    Particle *part = new Particle(glm::dvec2(xVal, yVal), 30.);
                    part->sFriction = 1;
                    part->kFriction = 1;
                    vertices.append(part);
                }
            }
            Body *body = createRigidBody(&vertices, &data);
            vertices.clear();
        }
    }

    double scale = 6., delta = .4;
    QList<Particle *> particles;

    double num = 1.;
    double start = m_xBoundaries.x + 1;
    for(double x = start; x < start + (scale / num); x += delta) {
        for(double y = 0; y < 1.2 * scale; y += delta) {
            particles.append(new Particle(glm::dvec2(x,y) + .2 * glm::dvec2(frand() - .5, frand() - .5), 1));
        }
    }
    createFluid(&particles, 2.5);
    particles.clear();

    int idx = m_particles.size();
    m_particles.append(new Particle(glm::dvec2(10, 50), 0));
    data.clear();

    glm::dvec2 base = glm::dvec2(57, 50);
    particles.append(new Particle(base, 1000));
    for (double a = 0; a <= 360; a+=30) {
        glm::dvec2 vec = glm::dvec2(cos(D2R(a)), sin(D2R(a)));
        particles.append(new Particle(vec * PARTICLE_RAD + base, 1000));
        data.append(SDFData(vec, PARTICLE_RAD * 1.5));
    }
    data.append(SDFData());
    createRigidBody(&particles, &data);

    m_globalConstraints[STANDARD].append(new DistanceConstraint(idx, idx + 1, &m_particles));
}

int Simulation::getNumParticles()
{
    return m_particles.size();
}

double Simulation::getKineticEnergy()
{
    double energy = 0;
    for (int i = 0; i < m_particles.size(); i++) {
        Particle *p = m_particles[i];
        if (p->imass != 0.) {
            energy += .5 * glm::dot(p->v, p->v) / p->imass;
        }
    }
    return energy;
}

void Simulation::mousePressed(const glm::dvec2 &p)
{
    for (int i = 0; i < m_particles.size(); i++) {
        Particle *part = m_particles.at(i);

        glm::dvec2 to = glm::normalize(p - part->p);
        part->v += 7. * to;
    }
    m_point = p;
}
