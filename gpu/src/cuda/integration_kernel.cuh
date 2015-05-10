#ifndef INTEGRATION_KERNEL_H
#define INTEGRATION_KERNEL_H

#include <stdio.h>
#include <math.h>
#include <curand.h>
#include <thrust/sort.h>

#include "helper_math.h"
#include "math_constants.h"
#include "kernel.cuh"
#include "shared_variables.cuh"

//#define X_BOUNDARY 7.f
//#define X_BOUNDARY 50.f
//#define Z_BOUNDARY 50.f

//#define

#define EPS 0.001f

////////////// fluid constants /////////////
#define MAX_FLUID_NEIGHBORS 500

#define H 2.f       // kernel radius
#define H2 4.f      // H^2
#define H6 64.f     // H^6
#define H9 512.f    // H^9
#define POLY6_COEFF 0.00305992474f // 315 / (64 * pi * H9)
#define SPIKEY_COEFF 0.22381163872f // 45 / (pi * H6)

#define FLUID_RELAXATION .01f // epsilon used when calculating lambda
#define K_P .1f              // scales artificial pressure
#define E_P 4.f              // exponent to art. pressure
#define DQ_P .2f             // between .1 and .3 (for art pressure)


/////////////////// friction ///////////////
#define S_FRICTION .005f
#define K_FRICTION .0002f
//#define S_FRICTION .15f
//#define K_FRICTION .003f

// textures for particle position and velocity
texture<float4, 1, cudaReadModeElementType> oldPosTex;
texture<float, 1, cudaReadModeElementType> invMassTex;
texture<int, 1, cudaReadModeElementType> oldPhaseTex;

texture<uint, 1, cudaReadModeElementType> gridParticleHashTex;
texture<uint, 1, cudaReadModeElementType> cellStartTex;
texture<uint, 1, cudaReadModeElementType> cellEndTex;


// simulation parameters in constant memory
__constant__ SimParams params;

struct collide_world_functor
{
    float *rands;
    int3 minBounds;
    int3 maxBounds;

    __host__ __device__
    collide_world_functor(float *_rands, int3 _minBounds, int3 _maxBounds)
        : rands(_rands), minBounds(_minBounds), maxBounds(_maxBounds) {}

    template <typename Tuple>
    __device__
    void operator()(Tuple t)
    {
        float4 posData = thrust::get<0>(t);
        float4 Xstar = thrust::get<1>(t);
        int phase = thrust::get<2>(t);

        float3 epos = make_float3(posData.x, posData.y, posData.z);
        float3 pos = make_float3(Xstar.x, Xstar.y, Xstar.z);

        float3 n = make_float3(0.f);

        float d = params.particleRadius;
        float eps = d * 0.f;
        if (phase < SOLID)
            eps = d * 0.01f;

        if (epos.y < minBounds.y + params.particleRadius)
        {
            epos.y = minBounds.y + params.particleRadius + rands[5] * eps;
            n += make_float3(0,1,0);
        }

        eps = d * 0.01f;

        if (epos.x > maxBounds.x - params.particleRadius)
        {
            epos.x = maxBounds.x - (params.particleRadius + rands[0] * eps);
            n += make_float3(-1,0,0);
        }

        if (epos.x < minBounds.x + params.particleRadius)
        {
            epos.x = minBounds.x + (params.particleRadius + rands[1] * eps);
            n += make_float3(1,0,0);
        }

        if (epos.y > maxBounds.y - params.particleRadius)
        {
            epos.y = maxBounds.y - (params.particleRadius + rands[2] * eps);
            n += make_float3(0,-1,0);
        }

#ifndef TWOD
        if (epos.z > maxBounds.z - params.particleRadius)
        {
            epos.z = maxBounds.z - (params.particleRadius + rands[3] * eps);
            n += make_float3(0,0,-1);
        }

        if (epos.z < minBounds.z + params.particleRadius)
        {
            epos.z = minBounds.z + (params.particleRadius + rands[4] * eps);
            n += make_float3(0,0,1);
        }
#endif


#ifdef TWOD
        epos.z = ZPOS; // 2D
        pos.z = ZPOS;
#endif

        if (length(n) < EPS || phase < SOLID)
        {
            thrust::get<0>(t) = make_float4(epos, posData.w);
            return;
        }

        float3 dp = (epos - pos);
        float3 dpt = dp - dot(dp, n) * n;
        float ldpt = length(dpt);

        if (ldpt < EPS)
        {
            thrust::get<0>(t) = make_float4(epos, posData.w);
            return;
        }


        if (ldpt < sqrt(S_FRICTION) * d)
            epos -= dpt;
        else
            epos -= dpt * min(sqrt(K_FRICTION) * d / ldpt, 1.f);

        // store new position and velocity

        thrust::get<0>(t) = make_float4(epos, posData.w);
    }
};

struct integrate_functor
{
    float deltaTime;

    __host__ __device__
    integrate_functor(float delta_time)
        : deltaTime(delta_time) {}

    template <typename Tuple>
    __device__
    void operator()(Tuple t)
    {
        volatile float4 posData = thrust::get<0>(t);
        volatile float4 velData = thrust::get<1>(t);
        float3 pos = make_float3(posData.x, posData.y, posData.z);
        float3 vel = make_float3(velData.x, velData.y, velData.z);

        vel += params.gravity * deltaTime;

        // new position = old position + velocity * deltaTime
        pos += vel * deltaTime;

        // store new position and velocity
        thrust::get<0>(t) = make_float4(pos, posData.w);
    }
};

// calculate position in uniform grid
__device__ int3 calcGridPos(float3 p)
{
    int3 gridPos;
    gridPos.x = floor((p.x - params.worldOrigin.x) / params.cellSize.x);
    gridPos.y = floor((p.y - params.worldOrigin.y) / params.cellSize.y);
    gridPos.z = floor((p.z - params.worldOrigin.z) / params.cellSize.z);
    return gridPos;
}

// calculate address in grid from position (clamping to edges)
__device__ uint calcGridHash(int3 gridPos)
{
    gridPos.x = gridPos.x & (params.gridSize.x-1);  // wrap grid, assumes size is power of 2
    gridPos.y = gridPos.y & (params.gridSize.y-1);
    gridPos.z = gridPos.z & (params.gridSize.z-1);
    return __umul24(__umul24(gridPos.z, params.gridSize.y), params.gridSize.x) + __umul24(gridPos.y, params.gridSize.x) + gridPos.x;
}

// calculate grid hash value for each particle
__global__
void calcHashD(uint   *gridParticleHash,  // output
               uint   *gridParticleIndex, // output
               float4 *pos,               // input: positions
               uint    numParticles)
{
    uint index = __umul24(blockIdx.x, blockDim.x) + threadIdx.x;

    if (index >= numParticles) return;

    volatile float4 p = pos[index];

    // get address in grid
    int3 gridPos = calcGridPos(make_float3(p.x, p.y, p.z));
    uint hash = calcGridHash(gridPos);

    // store grid hash and particle index
    gridParticleHash[index] = hash;
    gridParticleIndex[index] = index;
}

// rearrange particle data into sorted order, and find the start of each cell
// in the sorted hash array
__global__
void reorderDataAndFindCellStartD(uint   *cellStart,        // output: cell start index
                                  uint   *cellEnd,          // output: cell end index
                                  float4 *sortedPos,        // output: sorted positions
                                  float  *sortedW,          // output: sorted inverse masses
                                  int    *sortedPhase,      // output: sorted phase values
                                  uint   *gridParticleHash, // input: sorted grid hashes
                                  uint   *gridParticleIndex,// input: sorted particle indices
                                  float4 *oldPos,           // input: position array
                                  float  *W,
                                  int    *phase,
                                  uint    numParticles)
{
    extern __shared__ uint sharedHash[];    // blockSize + 1 elements
    uint index = __umul24(blockIdx.x,blockDim.x) + threadIdx.x;

    uint hash;

    // handle case when no. of particles not multiple of block size
    if (index < numParticles)
    {
        hash = gridParticleHash[index];

        // Load hash data into shared memory so that we can look
        // at neighboring particle's hash value without loading
        // two hash values per thread
        sharedHash[threadIdx.x+1] = hash;

        if (index > 0 && threadIdx.x == 0)
        {
            // first thread in block must load neighbor particle hash
            sharedHash[0] = gridParticleHash[index-1];
        }
    }

    __syncthreads();

    if (index < numParticles)
    {
        // If this particle has a different cell index to the previous
        // particle then it must be the first particle in the cell,
        // so store the index of this particle in the cell.
        // As it isn't the first particle, it must also be the cell end of
        // the previous particle's cell

        if (index == 0 || hash != sharedHash[threadIdx.x])
        {
            cellStart[hash] = index;

            if (index > 0)
                cellEnd[sharedHash[threadIdx.x]] = index;
        }

        if (index == numParticles - 1)
        {
            cellEnd[hash] = index + 1;
        }

        // Now use the sorted index to reorder the pos and vel data
        uint sortedIndex = gridParticleIndex[index];
        float4 pos = FETCH(oldPos, sortedIndex);       // macro does either global read or texture fetch
        float w = FETCH(invMass, sortedIndex);       // macro does either global read or texture fetch
        int phase = FETCH(oldPhase, sortedIndex);       // macro does either global read or texture fetch

        sortedPos[index] = pos;
        sortedW[index] = w;
        sortedPhase[index] = phase;
    }


}


// collide a particle against all other particles in a given cell
__device__
void collideCell(int3    gridPos,
                 uint    index,
                 float3  pos,
                 int     phase,
                 float4 *oldPos,
                 uint   *cellStart,
                 uint   *cellEnd,
                 uint   *neighbors,
                 uint   *numNeighbors)
{
    uint gridHash = calcGridHash(gridPos);

    // get start of bucket for this cell
    uint startIndex = FETCH(cellStart, gridHash);

    float collideDist = params.particleRadius * 2.001f; // slightly bigger radius
    float collideDist2 = collideDist * collideDist;

//    float3 delta = make_float3(0.0f);

    if (startIndex != 0xffffffff)          // cell is not empty
    {
        // iterate over particles in this cell
        uint endIndex = FETCH(cellEnd, gridHash);

        for (uint j=startIndex; j<endIndex; j++)
        {
            if (j != index)                // check not colliding with self
            {
                float3 pos2 = make_float3(FETCH(oldPos, j));
                int phase2 = FETCH(oldPhase, j);

                if (phase > SOLID && phase == phase2)
                    continue;

                // collide two spheres
                float3 diff = pos - pos2;

                float mag2 = dot(diff, diff);

                if (mag2 < collideDist2 && numNeighbors[index] < MAX_FLUID_NEIGHBORS)
                {
                    // neighbor stuff
                    neighbors[index * MAX_FLUID_NEIGHBORS + numNeighbors[index]] = j;
                    numNeighbors[index] += 1;

//                    delta += diff * (sqrt(mag2) - collideDist) * -.5f;
                }
            }
        }
    }
}


__global__
void collideD(float4 *newPos,               // output: new pos
              float4 *prevPositions,
              float4 *sortedPos,               // input: sorted positions
              float  *sortedW,
              int    *sortedPhase,
              uint   *gridParticleIndex,    // input: sorted particle indices
              uint   *cellStart,
              uint   *cellEnd,
              uint    numParticles,
              uint   *neighbors,
              uint   *numNeighbors)
{
    uint index = __mul24(blockIdx.x,blockDim.x) + threadIdx.x;

    if (index >= numParticles) return;

    int phase = FETCH(oldPhase, index);
    if (phase < SOLID) return;

    // read particle data from sorted arrays
    float3 pos = make_float3(FETCH(oldPos, index));

    // get address in grid
    int3 gridPos = calcGridPos(pos);

    // examine neighbouring cells
    float3 delta = make_float3(0.f);

    numNeighbors[index] = 0;
    for (int z=-1; z<=1; z++)
    {
        for (int y=-1; y<=1; y++)
        {
            for (int x=-1; x<=1; x++)
            {
                int3 neighbourPos = gridPos + make_int3(x, y, z);
                collideCell(neighbourPos, index, pos, phase, sortedPos, cellStart, cellEnd, neighbors, numNeighbors);
            }
        }
    }

    float collideDist = params.particleRadius * 2.001f;

    float w = FETCH(invMass, index);
    float sW = (w != 0.f ? (1.f / ((1.f / w) * exp(-pos.y))) : w);

    uint originalIndex = gridParticleIndex[index];
//    float3 currPos = make_float3(newPos[originalIndex]);
    float3 prevPos = make_float3(prevPositions[originalIndex]);

    for (uint i = 0; i < numNeighbors[index]; i++)
    {
        float3 pos2 =  make_float3(FETCH(oldPos, neighbors[index * MAX_FLUID_NEIGHBORS + i]));
        float w2 =  FETCH(invMass, neighbors[index * MAX_FLUID_NEIGHBORS + i]);
        int phase2 =  FETCH(oldPhase, neighbors[index * MAX_FLUID_NEIGHBORS + i]);

        float3 diff = pos - pos2;
        float dist = length(diff);
        float mag = dist - collideDist;

        float colW = w;
        float colW2 = w2;

        if (phase >= SOLID && phase2 >= SOLID)
        {
            colW = sW;
            colW2 = (w2 != 0.f ? (1.f / ((1.f / w2) * exp(-pos.y))) : w2);
        }

//        colWsum = colW + colW1);
        float scale = mag / (colW + colW2);
        float3 dp = diff * (scale / dist);
        float3 dp1 = -colW * dp / numNeighbors[index];
        float3 dp2 = colW2 * dp / numNeighbors[index];

        delta += dp1;



        ////////////////////// friction //////////////////
        if (phase < SOLID || phase2 < SOLID)
            continue;

        uint neighborIndex = gridParticleIndex[neighbors[index * MAX_FLUID_NEIGHBORS + i]];
        float3 prevPos2 = make_float3(prevPositions[neighborIndex]);
//        float3 currPos2 = make_float3(newPos[neighbors[index * MAX_FLUID_NEIGHBORS + i]]);

        float3 nf = normalize(diff);
        float3 dpRel = (pos + dp1 - prevPos) - (prevPos + dp2 - prevPos2);
        float3 dpt = dpRel - dot(dpRel, nf) * nf;
        float ldpt = length(dpt);

        if (ldpt < EPS)
            continue;

        if (ldpt < (S_FRICTION) * dist)
            delta -= dpt * colW / (colW + colW2);
        else
            delta -= dpt * min((K_FRICTION) * dist / ldpt, 1.f);
    }

    // write new velocity back to original unsorted location
    newPos[originalIndex] = make_float4(pos + delta, 1.0f);
}


struct subtract_functor
{
    const float time;

    subtract_functor(float _time) : time(_time) {}

    __device__
    float4 operator()(const float4& orig, const float4& solved) const {
        return (solved - orig) / -time;
    }
};




// collide a particle against all other particles in a given cell
__device__
void collideCellRadius(int3    gridPos,
                         uint    index,
                         float3  pos,
                         float4 *oldPos,
                         uint   *cellStart,
                         uint   *cellEnd,
                         uint   *neighbors,
                         uint   *numNeighbors)
{
    uint gridHash = calcGridHash(gridPos);

    // get start of bucket for this cell
    uint startIndex = FETCH(cellStart, gridHash);

    if (startIndex != 0xffffffff)          // cell is not empty
    {
        // iterate over particles in this cell
        uint endIndex = FETCH(cellEnd, gridHash);

        for (uint j=startIndex; j<endIndex; j++)
        {
            if (j != index)                // check not colliding with self
            {
                float3 pos2 = make_float3(FETCH(oldPos, j));

                float3 relPos = pos - pos2;
                float dist2 = dot(relPos, relPos);
                if (dist2 < H2 && numNeighbors[index] < MAX_FLUID_NEIGHBORS)
                {
                    // neighbor stuff
                    neighbors[index * MAX_FLUID_NEIGHBORS + numNeighbors[index]] = j;
                    numNeighbors[index] += 1;
                }
            }
        }
    }

}


__global__
void findLambdasD(float  *lambda,
                  float4 *oldPos,               // input: sorted positions
                  uint   *gridParticleIndex,    // input: sorted particle indices
                  uint   *cellStart,
                  uint   *cellEnd,
                  uint    numParticles,
                  uint   *neighbors,
                  uint   *numNeighbors,
                  float  *ros)
{
    uint index = __mul24(blockIdx.x,blockDim.x) + threadIdx.x;

    if (index >= numParticles) return;

    int phase = FETCH(oldPhase, index);
    if (phase != FLUID) return;

    // read particle data from sorted arrays
    float3 pos = make_float3(FETCH(oldPos, index));

    // get address in grid
    int3 gridPos = calcGridPos(pos);

    // examine neighbouring cells

    int rad = (int)ceil(H / params.cellSize.x);

    numNeighbors[index] = 0;
    for (int z=-rad; z<=rad; z++)
    {
        for (int y=-rad; y<=rad; y++)
        {
            for (int x=-rad; x<=rad; x++)
            {
                int3 neighbourPos = gridPos + make_int3(x, y, z);
                collideCellRadius(neighbourPos, index, pos, oldPos, cellStart, cellEnd, neighbors, numNeighbors);
            }
        }
    }

    float ro = 0.f;
    float denom = 0.f;
    float3 grad = make_float3(0.f);
    for (uint i = 0; i < numNeighbors[index]; i++)
    {
        float3 pos2 =  make_float3(FETCH(oldPos, neighbors[index * MAX_FLUID_NEIGHBORS + i]));
        float3 r = pos - pos2;
        float rlen2 = dot(r, r);
        float rlen = sqrt(rlen2);
        float hMinus2 = H2 - rlen2;
        float hMinus = H - rlen;

        // do fluid solid scaling hurr
        ro += (POLY6_COEFF * hMinus2*hMinus2*hMinus2 ) / 1.f; // <-- should be mass

        float3 spikeyGrad;
        if (rlen < 0.0001f)
            spikeyGrad = make_float3(0.f); // randomize a little
        else
            spikeyGrad = (r / rlen) * -SPIKEY_COEFF * hMinus*hMinus;
        spikeyGrad /= ros[gridParticleIndex[index]];

        grad += -spikeyGrad;
        denom += dot(spikeyGrad, spikeyGrad);
    }
    ro += (POLY6_COEFF * H6 ) / 1.f; // <-- should be mass too
    denom += dot(grad, grad);

    lambda[index] = - ((ro / ros[gridParticleIndex[index]]) - 1) / (denom + FLUID_RELAXATION);

//    float3 relPosM = make_float3(mousePos) - pos;
//    float distM = length(relPosM);

//    if (distM < params.particleRadius)
//    {
//        for (uint i = 0; i < numNeighbors[index]; i++)
//        {
//            neighbors[index * MAX_FLUID_NEIGHBORS + i] = gridParticleIndex[neighbors[index * MAX_FLUID_NEIGHBORS + i]];
//        }
//    }
//    else
//    {
//        numNeighbors[index] = 0;
//    }
}


__global__
void solveFluidsD(float  *lambda,
                  float4 *oldPos,               // input: sorted positions
                  uint   *gridParticleIndex,    // input: sorted particle indices
                  float4 *particles,
                  uint    numParticles,
                  uint   *neighbors,
                  uint   *numNeighbors,
                  float  *ros)
{
    uint index = __mul24(blockIdx.x,blockDim.x) + threadIdx.x;

    if (index >= numParticles) return;

    int phase = FETCH(oldPhase, index);
    if (phase != FLUID) return;

    float4 pos = FETCH(oldPos, index);

    float4 delta = make_float4(0.f);
    for (uint i = 0; i < numNeighbors[index]; i++)
    {
        float4 pos2 =  FETCH(oldPos, neighbors[index * MAX_FLUID_NEIGHBORS + i]);
        float4 r = pos - pos2;
        float rlen2 = dot(r, r);
        float rlen = sqrt(rlen2);
        float hMinus2 = H2 - rlen2;
        float hMinus = H - rlen;

        float4 spikeyGrad;
        if (rlen < 0.0001f)
            spikeyGrad = make_float4(0,EPS,0,0) * -SPIKEY_COEFF * hMinus*hMinus;
        else
            spikeyGrad = (r / rlen) * -SPIKEY_COEFF * hMinus*hMinus;

        float term2 = H2 - (DQ_P * DQ_P * H2);

        float numer = (POLY6_COEFF * hMinus2*hMinus2*hMinus2 ) ;
        float denom = (POLY6_COEFF * term2*term2*term2 );
        float lambdaCorr = -K_P * pow(numer / denom, E_P);

        delta += (lambda[index] + lambda[neighbors[index * MAX_FLUID_NEIGHBORS + i]] + lambdaCorr) * spikeyGrad;
    }

    uint origIndex = gridParticleIndex[index];
    particles[origIndex] += delta / (ros[gridParticleIndex[index]] + numNeighbors[index]);

}

#endif // INTEGRATION_KERNEL_H

