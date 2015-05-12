/*
 * This class controls helper classes for rendering and
 * updating the particle system.
 *
 * Mouse events, key events, window sizing, timing updates,
 * and render calls are all received by this application then
 * passed to the corresponding render or particlesystem class.
 */

#include <cuda_runtime.h>
#include <QMouseEvent>
#include <QWheelEvent>
#include <QKeyEvent>
#include <random>
#include <unistd.h>

#include "particleapp.h"
#include "particlesystem.h"
#include "renderer.h"
#include "helper_math.h"
#include "util.cuh"

#define MAX_PARTICLES 15000 // (vbo size)
#define PARTICLE_RADIUS 0.25f
#define GRID_SIZE make_uint3(64, 64, 64) // 3D


ParticleApp::ParticleApp()
    : m_particleSystem(NULL),
      m_renderer(NULL),
      m_mouseDownL(false),
      m_mouseDownR(false),
      m_fluidEmmiterOn(false),
      m_timer(-1.f)
{
    cudaInit();

    m_particleSystem = new ParticleSystem(PARTICLE_RADIUS, GRID_SIZE, MAX_PARTICLES, make_int3(-50, 0, -50), make_int3(50, 200, 50), 5);
    m_renderer = new Renderer(m_particleSystem->getMinBounds(), m_particleSystem->getMaxBounds());
    m_renderer->createVAO(m_particleSystem->getCurrentReadBuffer(),
                          m_particleSystem->getParticleRadius());
    makeInitScene();
}


ParticleApp::~ParticleApp()
{
    if (m_particleSystem)
        delete m_particleSystem;
    if (m_renderer)
        delete m_renderer;

    // cudaDeviceReset causes the driver to clean up all state. While
    // not mandatory in normal operation, it is good practice.  It is also
    // needed to ensure correct operation when the application is being
    // profiled. Calling cudaDeviceReset causes all profile data to be
    // flushed before the application exits
    cudaDeviceReset();
}

inline float frand()
{
    return rand() / (float) RAND_MAX;
}

void ParticleApp::makeInitScene()
{
    m_particleSystem->addRope(make_float3(0, 20, 0), make_float3(0, -.5, 0), .4f, 32, 1.f, true);
}


void ParticleApp::tick(float secs)
{
    if (m_fluidEmmiterOn && m_timer <= 0.f)
    {
        m_particleSystem->addFluid(make_int3(-1,0,-1), make_int3(1,1,1), 1.f, 1.f, make_float3(0,0,1));
        m_timer = 0.1f;
    }
    m_timer -= secs;

    m_particleSystem->update(secs);
    m_renderer->update(secs);
}


void ParticleApp::render()
{
    m_renderer->render(m_particleSystem->getColorIndex(), m_particleSystem->getColors());
}


void ParticleApp::mousePressed(QMouseEvent *e, float x, float y)
{
    // shoot a particle into the sceen on left mouse click
    if (e->button() == Qt::LeftButton)
    {
        m_particleSystem->setParticleToAdd(m_renderer->getEye(), m_renderer->getDir(x, y) * 30.f, 2.f);
        m_mouseDownL = true;
    }
    else if (e->button() == Qt::RightButton)
    {
        m_mouseDownR = true;
    }

}

void ParticleApp::mouseReleased(QMouseEvent *e, float, float)
{
    if (e->button() == Qt::LeftButton)
    {
        m_mouseDownL = false;
    }
    if (e->button() == Qt::RightButton)
    {
        m_mouseDownR = false;
    }
}


void ParticleApp::mouseMoved(QMouseEvent *e, float deltaX, float deltaY)
{
    m_renderer->mouseMoved(e, deltaX, deltaY);
}

void ParticleApp::mouseScrolled(QWheelEvent *) {}


void ParticleApp::keyPressed(QKeyEvent *e)
{
    m_renderer->keyPressed(e);
}

void ParticleApp::keyReleased(QKeyEvent *e)
{
    bool resetVbo = true;
    float3 h, vec;
    float angle;


    // numbers 0-9 toggle different scenes
    switch (e->key())
    {
    case Qt::Key_1: // single rope
        delete m_particleSystem;
        m_particleSystem = new ParticleSystem(PARTICLE_RADIUS, GRID_SIZE, MAX_PARTICLES, make_int3(-50, 0, -50), make_int3(50, 200, 50), 5);
        makeInitScene();
        break;
    case Qt::Key_2: // single cloth
        delete m_particleSystem;
        m_particleSystem = new ParticleSystem(PARTICLE_RADIUS, GRID_SIZE, MAX_PARTICLES, make_int3(-50, 0, -50), make_int3(50, 200, 50), 5);
        m_particleSystem->addHorizCloth(make_int2(0, -3), make_int2(6,3), make_float3(.5f,7.f,.5f), make_float2(.3f, .3f), 3.f, false);
        break;
    case Qt::Key_3: // two fluids, different densities
        delete m_particleSystem;
        m_particleSystem = new ParticleSystem(PARTICLE_RADIUS, GRID_SIZE, MAX_PARTICLES, make_int3(-7, 0, -5), make_int3(7, 20, 5), 5);
        m_particleSystem->addFluid(make_int3(-7, 0, -5), make_int3(7, 5, 5), 1.f, 2.f, colors[rand() % numColors]);
        m_particleSystem->addFluid(make_int3(-7, 5, -5), make_int3(7, 10, 5), 1.f, 3.f, colors[rand() % numColors]);
        break;
    case Qt::Key_4: // one solid particle stack
        delete m_particleSystem;
        m_particleSystem = new ParticleSystem(PARTICLE_RADIUS, GRID_SIZE, MAX_PARTICLES, make_int3(-50, 0, -50), make_int3(50, 200, 50), 5);
        m_particleSystem->addParticleGrid(make_int3(-3, 0, -3), make_int3(3, 20, 3), 1.f, false);
        break;
    case Qt::Key_5: // three solid particle stacks
        delete m_particleSystem;
        m_particleSystem = new ParticleSystem(PARTICLE_RADIUS, GRID_SIZE, MAX_PARTICLES, make_int3(-50, 0, -50), make_int3(50, 200, 50), 5);
        m_particleSystem->addParticleGrid(make_int3(-10, 0, -3), make_int3(-7, 10, 3), 1.f, false);
        m_particleSystem->addParticleGrid(make_int3(-3, 0, -3), make_int3(3, 10, 3), 1.f, false);
        m_particleSystem->addParticleGrid(make_int3(7, 0, -3), make_int3(10, 10, 3), 1.f, false);
        break;
    case Qt::Key_6: // particles on cloth
        delete m_particleSystem;
        m_particleSystem = new ParticleSystem(PARTICLE_RADIUS, GRID_SIZE, MAX_PARTICLES, make_int3(-50, 0, -50), make_int3(50, 200, 50), 5);
        m_particleSystem->addHorizCloth(make_int2(-10, -10), make_int2(10, 10), make_float3(.3f, 5.5f, .3f), make_float2(.1f, .1f), 10.f, true);
        m_particleSystem->addParticleGrid(make_int3(-3, 6, -3), make_int3(3, 15, 3), 1.f, false);
        break;
    case Qt::Key_7: // fluid blob
        delete m_particleSystem;
        m_particleSystem = new ParticleSystem(PARTICLE_RADIUS, GRID_SIZE, MAX_PARTICLES, make_int3(-50, 0, -50), make_int3(50, 200, 50), 5);
        m_particleSystem->addFluid(make_int3(-7, 6, -7), make_int3(7, 13, 7), 1.f, 1.5f, colors[rand() % numColors]);
        break;
    case Qt::Key_8: // combo scene
        delete m_particleSystem;
        m_particleSystem = new ParticleSystem(PARTICLE_RADIUS, GRID_SIZE, MAX_PARTICLES, make_int3(-50, 0, -50), make_int3(50, 200, 50), 5);
        m_particleSystem->addHorizCloth(make_int2(14, -4), make_int2(24, 6), make_float3(.3f, 2.5f, .3f), make_float2(.25f, .25f), 10.f, true);
        m_particleSystem->addHorizCloth(make_int2(10, -10), make_int2(25, -5), make_float3(.3f, 15.5f, .3f), make_float2(.25f, .25f), 3.f, false);
        m_particleSystem->addRope(make_float3(-17, 20, -17), make_float3(0, -.5, 0.001f), .4f, 30, 1.f, true);
        m_particleSystem->addRope(make_float3(-16, 20, -17), make_float3(0, 0, .5f), .4f, 50, 1.f, true);
        m_particleSystem->addRope(make_float3(-17, 20, -16), make_float3(0, -.5, 0.001f), .4f, 40, 1.f, true);
        m_particleSystem->addParticleGrid(make_int3(17, 6, 0), make_int3(21, 11, 4), 1.f, false);
        m_particleSystem->addParticleGrid(make_int3(-12, 0, -20), make_int3(0, 12, -17), 1.f, false);
        m_particleSystem->addParticleGrid(make_int3(-18, 0, -15), make_int3(-16, 9, -12), 1.f, false);
        m_particleSystem->addStaticSphere(make_int3(5, 5, -10), make_int3(10, 10, -5), .5f);
        break;
    case Qt::Key_9: // ropes on immovable sphere
        delete m_particleSystem;
        m_particleSystem = new ParticleSystem(PARTICLE_RADIUS, GRID_SIZE, MAX_PARTICLES, make_int3(-50, 0, -50), make_int3(50, 200, 50), 5);

        h = make_float3(0, 10, 0);
        for(int i = 0; i < 50; i++)
        {
            angle = M_PI * i * 0.02f;
            vec = make_float3(cos(angle), sin(angle), 0.f);
            m_particleSystem->addRope(vec*5.f + h, vec*.5f, .35f, 30, 1.f, true);
        }
        m_particleSystem->addStaticSphere(make_int3(-4, 7, -4), make_int3(4, 16, 4), .5f);

        break;
    case Qt::Key_0: // empty scene
        delete m_particleSystem;
        m_particleSystem = new ParticleSystem(PARTICLE_RADIUS, GRID_SIZE, MAX_PARTICLES, make_int3(-50, 0, -50), make_int3(50, 200, 50), 5);
        break;
    case Qt::Key_Space: // toggle fluids at origin
        m_fluidEmmiterOn = !m_fluidEmmiterOn;
        break;
    default:
        resetVbo = false;
        m_renderer->keyReleased(e);
        break;
    }
    if (resetVbo)
    {
        m_renderer->createVAO(m_particleSystem->getCurrentReadBuffer(),
                              m_particleSystem->getParticleRadius());
    }
}


void ParticleApp::resize(int w, int h)
{
    m_renderer->resize(w, h);
}


