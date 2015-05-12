#ifndef PARTICLES_KERNEL_H
#define PARTICLES_KERNEL_H

#define FETCH(t, i) tex1Dfetch(t##Tex, i)

#include "vector_types.h"

// simulation parameters
struct SimParams
{
    float3 gravity;
    float globalDamping;
    float particleRadius;

    uint3 gridSize;
    unsigned int numCells;
    float3 worldOrigin;
    float3 cellSize;

    unsigned int numBodies;
    unsigned int maxParticlesPerCell;
};

#endif //PARTICLES_KERNEL_H
