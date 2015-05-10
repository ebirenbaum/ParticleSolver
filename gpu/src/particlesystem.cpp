#include "particlesystem.h"
#include "GL/glew.h"
#include <string.h>
#include <assert.h>
#include <math.h>


#include "wrappers.cuh"
#include "kernel.cuh"
#include "util.cuh"
#include "shared_variables.cuh"
#include "helper_math.h"
//#include "kernel.cuh"

//#include "debugprinting.h"

ParticleSystem::ParticleSystem(uint numParticles, float particleRadius, uint3 gridSize, uint maxParticles, int3 minBounds, int3 maxBounds, int iterations)
    : m_initialized(false),
      m_particleRadius(particleRadius),
      m_maxParticles(maxParticles),
      m_numParticles(numParticles),
//      m_dPos(0),
      m_posVbo(0),
      m_cuda_posvbo_resource(0),
      m_gridSize(gridSize),
      m_rigidIndex(0),
      m_minBounds(minBounds),
      m_maxBounds(maxBounds),
      m_solverIterations(iterations)
{
    m_numGridCells = m_gridSize.x * m_gridSize.y * m_gridSize.z;

    m_gridSortBits = 18;

    // set simulation parameters
    m_params.gridSize = m_gridSize;
    m_params.numCells = m_numGridCells;
    m_params.numBodies = m_numParticles;

    m_params.particleRadius = m_particleRadius;

    m_params.worldOrigin = make_float3(0.f, 0.f, 0.f);
    float cellSize = m_params.particleRadius * 2.0f;  // cell size equal to particle diameter
    m_params.cellSize = make_float3(cellSize);

//    m_params.spring = 0.f;
//    m_params.damping = 0.f;
//    m_params.shear = 0.f;
//    m_params.attraction = 0.f;
//    m_params.boundaryDamping = 0.f;
    m_params.spring = 0.5f;
    m_params.damping = 0.02f;
    m_params.shear = 0.01f;
    m_params.attraction = 0.0f;
    m_params.boundaryDamping = -0.5f;

//    m_params.gravity = make_float3(0.0f, -.0f, 0.0f);
    m_params.gravity = make_float3(0.0f, -9.8f, 0.0f);
    m_params.globalDamping = 1.0f;

    _init(numParticles, maxParticles);

    mousePos = make_float4(-1);
}


ParticleSystem::~ParticleSystem()
{
    _finalize();
}


inline float frand()
{
    return rand() / (float) RAND_MAX;
}

// step the simulation
void ParticleSystem::update(float deltaTime)
{
    assert(m_initialized);

//    std::cout << deltaTime << std::endl;
    deltaTime = std::min(deltaTime, .05f);
//    deltaTime = .01f;

    if (m_numParticles == 0)
    {
        addNewStuff();
        return;
    }

    float *dPos = (float *) mapGLBufferObject(&m_cuda_posvbo_resource);

    // update constants
    setParameters(&m_params);

    /* ALGORITHM 1 SIMULATION LOOP */

    /**********************************************************************************
     * 1    for all particles i do:
     * 2 (x)    apply forces vi = vi + deltaTime * forcesExt(xi)
     * 3 (x)    apply forces xi* = xi + deltaTime * vi
     * 4        apply mass scaling mi* = mi * exp(-k*h(xi*))
     * 5    end for
     */

    integrateSystem(dPos,
                    deltaTime,
                    m_numParticles);

    /**********************************************************************************
     * 6    for all particles i do:
     * 7        find neighboring particles Ni(xi*)
     * 8        find solid contacts
     * 9    end for
     */

//    // calculate grid hash
//    calcHash(   m_dGridParticleHash,
//                m_dGridParticleIndex,
//                dPos,
//                m_numParticles);

//    // sort particles based on hash
//    sortParticles(m_dGridParticleHash, m_dGridParticleIndex, m_numParticles);

//    // reorder particle arrays into sorted order and
//    // find start and end of each cell
//    reorderDataAndFindCellStart(
//                m_dCellStart,
//                m_dCellEnd,
//                m_dSortedPos,
//                m_dGridParticleHash,
//                m_dGridParticleIndex,
//                dPos,
//                m_numParticles,
//                m_numGridCells);


    //    /**********************************************************************************
    //     * 10   while iter < stabalizationIterations do:
    //     * 11       deltaX = 0, n = 0
    //     * 12       solve contact constraints for deltaX, n
    //     * 13       update x  =  x + deltaX / n
    //     * 14       update x* = x* + deltaX / n
    //     * 15   end while
    //     */

//    collide(    dPos,
//                m_dSortedPos,
//                m_dGridParticleIndex,
//                m_dCellStart,
//                m_dCellEnd,
//                m_numParticles,
//                m_numGridCells);

//    collideWorld(dPos,
//                 m_numParticles);

    /**********************************************************************************
     * 16   while iter < solverIterations do:
     * 17       for each constraint group G do:
     * 18           deltaX = 0, n = 0
     * 19           solve all constraints in G for deltaX, n
     * 24           update x* = x* + deltaX / n
     * 21       end for
     * 22   end while
     */

    sortByType(dPos,
               m_numParticles);

    for (uint i = 0; i < m_solverIterations; i++)
    {

        // calculate grid hash
        calcHash(   m_dGridParticleHash,
                    m_dGridParticleIndex,
                    dPos,
                    m_numParticles);

        // sort particles based on hash
        sortParticles(m_dGridParticleHash, m_dGridParticleIndex, m_numParticles);

        // reorder particle arrays into sorted order and
        // find start and end of each cell
        reorderDataAndFindCellStart(
                    m_dCellStart,
                    m_dCellEnd,
                    m_dSortedPos,
                    m_dSortedW,
                    m_dSortedPhase,
                    m_dGridParticleHash,
                    m_dGridParticleIndex,
                    dPos,
                    m_numParticles,
                    m_numGridCells);

        collide(    dPos,
                    m_dSortedPos,
                    m_dSortedW,
                    m_dSortedPhase,
                    m_dGridParticleIndex,
                    m_dCellStart,
                    m_dCellEnd,
                    m_numParticles,
                    m_numGridCells);

        solveFluids(    m_dSortedPos,
                        m_dSortedPhase,
                        m_dGridParticleIndex,
                        m_dCellStart,
                        m_dCellEnd,
                        dPos,
                        m_numParticles,
                        m_numGridCells,
                        mousePos);

        collideWorld(dPos,
                     m_dSortedPos,
                     m_numParticles,
                     m_minBounds,
                     m_maxBounds);


        solveDistanceConstraints(dPos);

        solvePointConstraints(dPos);

    }


    /**********************************************************************************
     * 23   for all particles i do:
     * 24       update velocity v = (xi* - xi) / deltaTime:
     * 25       advect diffuse particles
     * 26       apply internal forces fDrag, fVort
     * 27       update positions xi = xi* or apply sleeping
     * 28   end for
     */

    calcVelocity(dPos,
                 deltaTime,
                 m_numParticles);

    // note: do unmap at end here to avoid unnecessary graphics/CUDA context switch
    unmapGLBufferObject(m_cuda_posvbo_resource);

    addNewStuff();
}


void ParticleSystem::addNewStuff()
{
    addParticles();
    addFluids();
}


void ParticleSystem::addParticles()
{
    if (m_particlesToAdd.empty())
        return;

    uint size = m_particlesToAdd.size() / 2;
    float4 pos, vel;
    for (uint i = 0; i < size; i++)
    {
        pos = m_particlesToAdd.front();
        m_particlesToAdd.pop_front();
        vel = m_particlesToAdd.front();
        m_particlesToAdd.pop_front();
        addParticle(pos, make_float4(vel.x, vel.y, vel.z, 0), vel.w, 1.5f, SOLID);
    }
}

void ParticleSystem::addFluids()
{
    if (m_fluidsToAdd.empty())
        return;

    uint start = m_numParticles;

    uint size = m_fluidsToAdd.size() / 2;
    float4 pos, color;
    for (uint i = 0; i < size; i++)
    {
        pos = m_fluidsToAdd.front();
        m_fluidsToAdd.pop_front();
        color = m_fluidsToAdd.front();
        m_fluidsToAdd.pop_front();
        addParticle(make_float4(make_float3(pos), 1), make_float4(0,-1,0,0), pos.w, color.w, FLUID);
    }

    m_colorIndex.push_back(make_int2(start, m_numParticles));
    m_colors.push_back(make_float4(make_float3(color), 1.f));
}
void ParticleSystem::addFluidBlock()
{
//    if (m_particlesToAdd.empty())
//        return;

//    uint size = m_particlesToAdd.size() / 2;
//    int3 ll, ur;
//    for (uint i = 0; i < size; i++)
//    {
//        ll = m_fluidsToAdd.front();
//        m_fluidsToAdd.pop_front();
//        ur = m_fluidsToAdd.front();
//        m_fluidsToAdd.pop_front();
//        addFluid(ll, ur, 1.f, 1.5f);
//    }

}


void ParticleSystem::setParticleToAdd(float3 pos, float3 vel, float mass)
{
    float jitter = m_particleRadius * 0.01f;
    pos.x += (frand()*2.0f-1.0f) * jitter;
    pos.y += (frand()*2.0f-1.0f) * jitter;
    m_particlesToAdd.push_back(make_float4(pos, 1.f));
    m_particlesToAdd.push_back(make_float4(vel, mass));

    m_colorIndex.push_back(make_int2(m_numParticles, m_numParticles+1));
    m_colors.push_back(make_float4(frand(), frand(), frand(), 1.f));
}


void ParticleSystem::addParticle(float4 pos, float4 vel, float mass, float ro, int phase)
{
    if (m_numParticles == m_maxParticles)
        return;

    float *data = (float*)&pos;
    unregisterGLBufferObject(m_cuda_posvbo_resource);
    glBindBuffer(GL_ARRAY_BUFFER, m_posVbo);
    glBufferSubData(GL_ARRAY_BUFFER, m_numParticles*4*sizeof(float), 4*sizeof(float), data);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    registerGLBufferObject(m_posVbo, &m_cuda_posvbo_resource);

    appendIntegrationParticle(vel, ro, 1);
    appendPhaseAndMass(phase, mass, 1);
    appendSolverParticle();
    m_numParticles++;
}

void ParticleSystem::setFluidToAdd(float3 pos, float3 color, float mass, float density)
{
    m_fluidsToAdd.push_back(make_float4(pos, mass));
    m_fluidsToAdd.push_back(make_float4(color, density));
}


void ParticleSystem::addFluid(int3 ll, int3 ur, float mass, float density, float3 color)
{
    int start = m_numParticles;
    float jitter = m_particleRadius * 0.01f;
    float distance = m_particleRadius * 2.5f/*1.667 * density*/;

    int3 count = make_int3((int)ceil(ur.x - ll.x) / distance, (int)ceil(ur.y - ll.y) / distance, (int)ceil(ur.z - ll.z) / distance);

//    std::cout << count.x << ", " << count.y << ", " << count.x << ", " << std::endl;
    float4 pos;

#ifndef TWOD
    for (int z = 0; z < count.z; z++)
    {
#endif
        for (int y = 0; y < count.y; y++)
        {
            for (int x = 0; x < count.x; x++)
            {
//                std::cout << x << ", " << y << ", " << z << std::endl;

                pos = make_float4(ll.x + x * distance + (frand()*2.0f-1.0f)*jitter,
                                  ll.y + y * distance + (frand()*2.0f-1.0f)*jitter,
            #ifdef TWOD
                                  ZPOS,
            #else
                                  ll.z + z * distance + (frand()*2.0f-1.0f)*jitter,
            #endif
                                  1.f);

                addParticle(pos, make_float4(0.f), mass, density, FLUID);
            }
        }
#ifndef TWOD
    }
#endif

    m_colorIndex.push_back(make_int2(start, m_numParticles));
    m_colors.push_back(make_float4(color, 1.f));
}


void ParticleSystem::addParticleGrid(int3 ll, int3 ur, float mass, bool addJitter)
{
    int start = m_numParticles;
    float jitter = 0.f;
    if (addJitter)
        jitter = m_particleRadius * 0.01f;
    float distance = m_particleRadius * 2.002f;

    int3 count = make_int3((int)ceil(ur.x - ll.x) / distance, (int)ceil(ur.y - ll.y) / distance, (int)ceil(ur.z - ll.z) / distance);

    std::cout << count.x << ", " << count.y << ", " << count.z << std::endl;
    float4 pos;

#ifndef TWOD
    for (int z = 0; z < count.z; z++)
    {
#endif
        for (int y = 0; y < count.y; y++)
        {
            for (int x = 0; x < count.x; x++)
            {
                pos = make_float4(ll.x + x * distance + (frand()*2.0f-1.0f)*jitter,
                                  ll.y + y * distance + (frand()*2.0f-1.0f)*jitter,
            #ifdef TWOD
                                  ZPOS,
            #else
                                  ll.z + z * distance + (frand()*2.0f-1.0f)*jitter,
            #endif
                                  1.f);

                addParticle(pos, make_float4(0.f), mass, 1.f, SOLID);
            }
        }
#ifndef TWOD
    }
#endif

    m_colorIndex.push_back(make_int2(start, m_numParticles));
    m_colors.push_back(make_float4(frand(), frand(), frand(), 1.f));
}


void ParticleSystem::addHorizCloth(int2 ll, int2 ur, float3 spacing, float2 dist, float mass, bool holdEdges)
{
    int start = m_numParticles;

    int2 count = make_int2((int)ceil(ur.x - ll.x) / spacing.x, (int)ceil(ur.y - ll.y) / spacing.z);

    std::cout << count.x << ", " << count.y << std::endl;
    float4 pos;

#ifndef TWOD
    for (int z = 0; z < count.y; z++)
    {
#endif
        for (int x = 0; x < count.x; x++)
        {
            pos = make_float4(ll.x + x * spacing.x,
                              spacing.y,
        #ifdef TWOD
                              ZPOS,
        #else
                              ll.y + z * spacing.z,
        #endif
                              1.f);

            addParticle(pos, make_float4(0.f), mass, 1.f, RIGID + m_rigidIndex);

            uint index = start + z * count.x + x;
            if (x > 0)
                addDistanceConstraint(make_uint2(index - 1, index), dist.x);
            else
                addPointConstraint(index, make_float3(pos));
            if (z > 0)
                addDistanceConstraint(make_uint2(index - count.x, index), dist.y);
            else if (holdEdges)
                addPointConstraint(index, make_float3(pos));
            if (x == count.x - 1 && holdEdges)
                addPointConstraint(index, make_float3(pos));
            if (z == count.y - 1 && holdEdges)
                addPointConstraint(index, make_float3(pos));


        }
#ifndef TWOD
    }
#endif

    m_colorIndex.push_back(make_int2(start, m_numParticles));
    m_colors.push_back(make_float4(frand(), frand(), frand(), 1.f));
    m_rigidIndex++;
}

void ParticleSystem::addRope(float3 start, float3 spacing, float dist, int numLinks, float mass, bool constrainStart)
{
    int startI = m_numParticles;

    addParticle(make_float4(start, 1.f), make_float4(0.f), mass, 1.f, RIGID + m_rigidIndex);
    if (constrainStart)
        addPointConstraint(startI, start);

    float4 pos;
    for (int i = 1; i < numLinks; i++)
    {
        pos = make_float4(start + i * spacing, 1.f);
        addParticle(pos, make_float4(0.f), mass, 1.f, RIGID + m_rigidIndex);
        addDistanceConstraint(make_uint2(startI + i-1, startI + i), dist);
    }

    m_colorIndex.push_back(make_int2(startI, m_numParticles));
    m_colors.push_back(make_float4(.5f, .4f, .2f, 1.f));
    m_rigidIndex++;
}

void ParticleSystem::addStaticSphere(int3 ll, int3 ur, float spacing)
{
    uint startI = m_numParticles;
    int3 count = make_int3((int)ceil(ur.x - ll.x) / spacing, (int)ceil(ur.y - ll.y) / spacing, (int)ceil(ur.z - ll.z) / spacing);
    float4 pos;
    float radius = (ur.x - ll.x) * .5f;
    float3 center = make_float3(ll) + make_float3(radius);

#ifndef TWOD
    for (int z = 0; z < count.z; z++)
    {
#endif
        for (int y = 0; y < count.y; y++)
        {
            for (int x = 0; x < count.x; x++)
            {

                pos = make_float4(ll.x + x * spacing,
                                  ll.y + y * spacing,
            #ifdef TWOD
                                  ZPOS,
            #else
                                  ll.z + z * spacing,
            #endif
                                  1.f);
                if (length(make_float3(pos) - center) < radius)
                {
                    addParticle(pos, make_float4(0.f), 10.f, 1.f, RIGID + m_rigidIndex);
                    addPointConstraint(m_numParticles-1, make_float3(pos));
                }
            }
        }
#ifndef TWOD
    }
#endif

    m_colorIndex.push_back(make_int2(startI, m_numParticles));
    m_colors.push_back(make_float4(frand(),frand(),frand(), 1.f));
    m_rigidIndex++;


}


void ParticleSystem::makePointConstraint(uint index, float3 point)
{
    addPointConstraint(index, point);
}

void ParticleSystem::makeDistanceConstraint(uint2 index, float distance)
{
    addDistanceConstraint(index, distance);
}


void ParticleSystem::setArray(bool isPosArray, const float *data, int start, int count)
{
    assert(m_initialized);

    if (isPosArray)
    {
        unregisterGLBufferObject(m_cuda_posvbo_resource);
        glBindBuffer(GL_ARRAY_BUFFER, m_posVbo);
        glBufferSubData(GL_ARRAY_BUFFER, start*4*sizeof(float), count*4*sizeof(float), data);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        registerGLBufferObject(m_posVbo, &m_cuda_posvbo_resource);
    }
}


GLuint ParticleSystem::createVBO(uint size)
{
    GLuint vbo;
    glGenBuffers(1, &vbo);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, size, 0, GL_DYNAMIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    return vbo;
}


void ParticleSystem::_init(uint numParticles, uint maxParticles)
{
    m_maxParticles = maxParticles;
    m_numParticles = numParticles;
    cudaInit();
    initIntegration();
    initHandles();

    /*
     *  allocate GPU data
     */
    uint memSize = sizeof(GLfloat) * 4 * m_maxParticles;

    m_posVbo = createVBO(sizeof(GLfloat) * 4 * m_maxParticles);
    registerGLBufferObject(m_posVbo, &m_cuda_posvbo_resource);

    // grid and collisions
    allocateArray((void **)&m_dSortedPos, memSize);
    allocateArray((void **)&m_dSortedW, m_maxParticles*sizeof(float));
    allocateArray((void **)&m_dSortedPhase, m_maxParticles*sizeof(int));

    allocateArray((void **)&m_dGridParticleHash, m_maxParticles*sizeof(uint));
    allocateArray((void **)&m_dGridParticleIndex, m_maxParticles*sizeof(uint));

    allocateArray((void **)&m_dCellStart, m_numGridCells*sizeof(uint));
    allocateArray((void **)&m_dCellEnd, m_numGridCells*sizeof(uint));


    setParameters(&m_params);

    m_initialized = true;
}


void ParticleSystem::_finalize()
{
    assert(m_initialized);

    freeArray(m_dSortedPos);
    freeArray(m_dSortedW);
    freeArray(m_dSortedPhase);

    freeArray(m_dGridParticleHash);
    freeArray(m_dGridParticleIndex);
    freeArray(m_dCellStart);
    freeArray(m_dCellEnd);

    unregisterGLBufferObject(m_cuda_posvbo_resource);
    glDeleteBuffers(1, (const GLuint *)&m_posVbo);

    freeIntegrationVectors();
    freeSolverVectors();
    freeSharedVectors();
    destroyHandles();
}

