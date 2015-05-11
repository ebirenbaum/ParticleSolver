#ifndef SOLVER_KERNEL_H
#define SOLVER_KERNEL_H

#include <stdio.h>
#include <math.h>
#include "helper_math.h"
#include "math_constants.h"
#include "thrust/tuple.h"

struct point_constraint_functor
{
    float4 *particles;

    __host__ __device__
    point_constraint_functor(float4 *particles_) : particles(particles_) {}

    template <typename Tuple>
    __device__
    void operator()(Tuple t)
    {
        uint index = thrust::get<0>(t);
        float4 pos = particles[index];
        particles[index] = make_float4(thrust::get<1>(t), pos.w);
    }
};

struct delta_computing_functor
{
    float4 *particles;

    __host__ __device__
    delta_computing_functor(float4 *particles_) : particles(particles_) {}

    template <typename Tuple>
    __device__
    void operator()(Tuple t)
    {
        /*
         * 0: index
         * 1: dists
         * 2: sortedI1
         * 3: sortedI2
         * 4: delta1
         * 5: delta2
         */
        uint2 index = thrust::get<0>(t);
        float4 p1 = particles[index.x];
        float4 p2 = particles[index.y];

        thrust::get<2>(t) = index.x;
        thrust::get<3>(t) = index.y;

        float4 relPos = p1 - p2;
        relPos.w = 0.f; // inverse masses not needed

        float dist = length(relPos);
        if (dist > 0.0001f)
        {
            float4 grad = relPos / dist;
            float mag = (thrust::get<1>(t) - dist) * .5f;
            float4 delta = grad * mag;

            thrust::get<4>(t) = delta;
            thrust::get<5>(t) = -delta;
        }
        else
        {
            thrust::get<4>(t) = make_float4(0);
            thrust::get<5>(t) = make_float4(0);
        }
    }
};

struct gradient_functor
{
    float4 *particles;
    float *val;
    float *JT;
    uint numParticles;
    uint numConstraints;

    __host__ __device__
    gradient_functor(float4 *particles_, float *val_, float *JT_, uint numParticles_, uint numConstraints_)
        : particles(particles_), val(val_), JT(JT_), numParticles(numParticles_), numConstraints(numConstraints_) {}

    template <typename Tuple>
    __device__
    void operator()(Tuple t)
    {
        uint constraintIndex = thrust::get<0>(t);
        uint p1i = thrust::get<1>(t);
        uint p2i = thrust::get<2>(t);
        float4 p14 = particles[p1i];
        float4 p24 = particles[p2i];
        float3 p1 = make_float3(p14.x, p14.y, p14.z);
        float3 p2 = make_float3(p24.x, p24.y, p24.z);

        float3 relPos = p1 - p2;

        // 2D
        relPos.z = 0.f;

        float dist = length(relPos);

        if (dist > 0.0001f)
        {
            float3 grad = relPos / dist;

            // set sparse matrix values for J
            if (p1i < p2i)
            {
                val[constraintIndex * 4] = grad.x;
                val[constraintIndex*4+1] = grad.y;
                val[constraintIndex*4+2] = -grad.x;
                val[constraintIndex*4+3] = -grad.y;
            }
            else
            {
                val[constraintIndex * 4] = -grad.x;
                val[constraintIndex*4+1] = -grad.y;
                val[constraintIndex*4+2] = grad.x;
                val[constraintIndex*4+3] = grad.y;
            }

            // set dense matrix values for JT
            uint col = constraintIndex * numParticles * 2;
            JT[col + p1i * 2] = grad.x;
            JT[col + p1i*2+1] = grad.y;
            JT[col + p2i * 2] = -grad.x;
            JT[col + p2i*2+1] = -grad.y;
        }
        else
        {
            // reset to zero if gradient can't be computed

            // sparse J
            val[constraintIndex * 4] = 0.f;
            val[constraintIndex*4+1] = 0.f;
            val[constraintIndex*4+2] = 0.f;
            val[constraintIndex*4+3] = 0.f;

            // dense JT
            uint col = constraintIndex * numParticles * 2;
            JT[col + p1i * 2] = 0.f;
            JT[col + p1i*2+1] = 0.f;
            JT[col + p2i * 2] = 0.f;
            JT[col + p2i*2+1] = 0.f;
        }

        // negative distance for vec b
        thrust::get<4>(t) = (thrust::get<3>(t) - dist);
    }
};

struct distance_constraint_functor
{
    template <typename Tuple>
    __device__
    void operator()(Tuple t)
    {
        thrust::get<0>(t) += thrust::get<1>(t) / thrust::get<2>(t);
    }
};

struct identity_functor
{
    const uint w;

    identity_functor(uint _w) : w(_w) {}

    __device__
    float operator()(const uint& x) const
    {
        if (x / w == x % w)
            return 1;
        return 0;
    }
};

struct diag_extraction_functor
{
    float *A;
    uint width;

    diag_extraction_functor(float *_A, uint _width)
        : A(_A), width(_width) {}

    template <typename Tuple>
    __device__
    void operator()(Tuple t)
    {
        uint index = thrust::get<0>(t);
        index = index * width + index;
        float elmt = A[index];
        thrust::get<1>(t) = (abs(elmt) < 0.0001f ? 0.f : 1.f / elmt);
        A[index] = 0.f;
    }
};

struct col_reduction
{
    float4 *particles;
    float *b;
    float *J;
    const uint rowsJ;
    const uint colsJ;

    __host__ __device__
    col_reduction(float4 *_particles, float *_b, float *_J, const uint _rowsJ, const uint _colsJ)
        : particles(_particles), b(_b), J(_J), rowsJ(_rowsJ), colsJ(_colsJ) {}


    template <typename Tuple>
    __device__
    void operator()(Tuple t)
    {
        uint row = thrust::get<0>(t);
        uint index = row*2;
        float4 temp = make_float4(0.f);
        for (uint i = 0; i < colsJ; i++)
        {
            temp.x += J[index] * b[i];
            temp.y += J[index+1] * b[i];
            index += rowsJ;
        }

        particles[row] += temp / (thrust::get<1>(t) * thrust::get<2>(t));
    }
};


__global__
void constraintCentricSolveD(float      *A,
                             float      *B,
                             float      *JT,
                             float      *lambdas,
                             float      *particles,
                             uint       *indices,
                             uint       *occurences,
                             uint       numParticles,
                             uint       numConstraints,
                             float      minErr,
                             uint       maxIter,
                             float      omega)
{
    uint index = __umul24(blockIdx.x,blockDim.x) + threadIdx.x;

    if (index < numConstraints)
    {
        int vars = numConstraints - 1;
        float bi = B[index];

        float err = 1.f;
        int iter = 0;
        float sum, lambda, old_lambda = 0.f;
        uint row;

        __syncthreads();
        while ((err > minErr /*|| iter < 10*/) && iter < maxIter)
        {
            sum = bi;
            row = index;
            for (int i = 0; i < vars; i++)
            {
                if (i != index)
                    sum -= lambdas[i] * A[row];
                row += numConstraints;              // next column of matrix A (stored in column-major order)
            }
            lambda = sum / A[index * numConstraints + index]; // divide by corresponding diagonal coefficient
            lambdas[index] = lambda;

            err = lambda - old_lambda;
            old_lambda = lambda;
            iter++;
        }

        uint pi = indices[index];
        atomicAdd(particles + pi * 4, omega * lambda * JT[index * numParticles * 2 + pi * 2] / occurences[pi]);
        atomicAdd(particles + pi*4+1, omega * lambda * JT[index * numParticles * 2 + pi*2+1] / occurences[pi]);

        pi = indices[index + numConstraints];
        atomicAdd(particles + pi * 4, omega * lambda * JT[index * numParticles * 2 + pi * 2] / occurences[pi]);
        atomicAdd(particles + pi*4+1, omega * lambda * JT[index * numParticles * 2 + pi*2+1] / occurences[pi]);
    }

    __syncthreads();
}

struct rigid_body_functor
{
    float2 *rVecs;
    float4 *particles;

    __host__ __device__
    rigid_body_functor(float2 *_rVecs, float4 *_particles)
        : rVecs(_rVecs), particles(_particles) {}

    template <typename Tuple>
    __device__
    void operator()(Tuple t)
    {
        /* Tuple:
         * 0: angle
         * 1: center mass
         * 2: particle index start
         * 3: rVec index start
         * 4: size
         */


     }
};


struct occurence_functor
{
    uint *occ;

    __host__ __device__
    occurence_functor(uint *_occ)
        : occ(_occ) {}

    template <typename Tuple>
    __device__
    void operator()(Tuple t)
    {
        occ[thrust::get<0>(t)] += thrust::get<1>(t);
    }
};

#endif // SOLVER_KERNEL_H

