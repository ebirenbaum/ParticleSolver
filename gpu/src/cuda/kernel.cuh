#ifndef PARTICLES_KERNEL_H
#define PARTICLES_KERNEL_H

//#define TWOD
#ifdef TWOD
#define ZPOS .5f
#endif

#define FETCH(t, i) tex1Dfetch(t##Tex, i)

#include "vector_types.h"

// simulation parameters
struct SimParams
{
//    float3 colliderPos;
//    float  colliderRadius;

    float3 gravity;
    float globalDamping;
    float particleRadius;

    uint3 gridSize;
    unsigned int numCells;
    float3 worldOrigin;
    float3 cellSize;

    unsigned int numBodies;
    unsigned int maxParticlesPerCell;


    float spring;
    float damping;
    float shear;
    float attraction;
    float boundaryDamping;

};

#endif //PARTICLES_KERNEL_H
