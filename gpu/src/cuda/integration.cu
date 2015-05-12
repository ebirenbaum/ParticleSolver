
#include <cuda_runtime.h>
#include <cuda_gl_interop.h>

#include <curand.h>
#include <stdio.h>

#include <thrust/device_ptr.h>
#include <thrust/device_vector.h>
#include <thrust/for_each.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/transform.h>

#include "helper_cuda.h"
#include "integration_kernel.cuh"
#include "util.cuh"
#include "shared_variables.cuh"

//#define PRINT

curandGenerator_t gen(0);

thrust::device_vector<float> V; // particle velocities
thrust::device_vector<float> lambda;
thrust::device_vector<float> denom;

thrust::device_vector<float> ros;

thrust::device_vector<uint> neighbors;
thrust::device_vector<uint> numNeighbors;

thrust::device_vector<float> textureVec;

float *rands;

extern "C"
{
    /*****************************************************************************
     *                              INITIALIZATION
     *****************************************************************************/


    void initIntegration()
    {
        allocateArray((void **)&rands, 6 * sizeof(float));
        checkCudaErrors(curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_DEFAULT));
        checkCudaErrors(curandSetPseudoRandomGeneratorSeed(gen, 1234ULL));
    }

    void appendIntegrationParticle(float *v, float *ro, uint numParticles)
    {
        int sizeV = V.size();
        int sizeRo = ros.size();

        // resize the vectors
        V.resize(sizeV + 4 * numParticles);
        ros.resize(sizeRo + numParticles);

        // get raw pointers to the data
        float *dV = thrust::raw_pointer_cast(V.data());
        float *dRos = thrust::raw_pointer_cast(ros.data());

        // copy the new data over to the gpu
        copyArrayToDevice(dV + sizeV, v, 0, 4 * numParticles * sizeof(float));
        copyArrayToDevice(dRos + sizeRo, ro, 0, numParticles * sizeof(float));

        // resize but don't need to fill
        lambda.resize(ros.size());
        numNeighbors.resize(ros.size());
        neighbors.resize(V.size() * MAX_FLUID_NEIGHBORS);
        textureVec.resize(V.size());
    }

    void freeIntegrationVectors()
    {
         V.clear();
         lambda.clear();
         denom.clear();
         ros.clear();
         neighbors.clear();
         numNeighbors.clear();
         textureVec.clear();

         V.shrink_to_fit();
         lambda.shrink_to_fit();
         denom.shrink_to_fit();
         ros.shrink_to_fit();
         neighbors.shrink_to_fit();
         numNeighbors.shrink_to_fit();
         textureVec.shrink_to_fit();

         checkCudaErrors(curandDestroyGenerator(gen));
         freeArray(rands);
    }

    void setParameters(SimParams *hostParams)
    {
        // copy parameters to constant memory
        checkCudaErrors(cudaMemcpyToSymbol(params, hostParams, 1 * sizeof(SimParams)));
    }

















    /*****************************************************************************
     *                              UPDATE POSITIONS
     *****************************************************************************/

    void integrateSystem(float *pos, float deltaTime, uint numParticles)
    {
        thrust::device_ptr<float4> d_pos4((float4 *)pos);
        thrust::device_ptr<float4> d_vel4((float4 *)thrust::raw_pointer_cast(V.data()));

        // copy current positions for reference later
        copyToXstar(pos, numParticles);

        // guess new positions based on forces
        thrust::for_each(
            thrust::make_zip_iterator(thrust::make_tuple(d_pos4, d_vel4)),
            thrust::make_zip_iterator(thrust::make_tuple(d_pos4+numParticles, d_vel4+numParticles)),
            integrate_functor(deltaTime));
    }





















    /*****************************************************************************
     *                              BUILD GRID
     *****************************************************************************/

    void calcHash(uint *gridParticleHash, uint *gridParticleIndex, float *pos, int numParticles)
    {
        uint numThreads, numBlocks;
        computeGridSize(numParticles, 256, numBlocks, numThreads);

        // execute the kernel
        calcHashD<<< numBlocks, numThreads >>>(gridParticleHash, gridParticleIndex, (float4 *) pos, numParticles);

#ifdef PRINT
        printf("HASHES:\n");
        thrust::device_ptr<uint> dGPH(gridParticleHash);
        for (uint i = 0; i < numParticles; i++)
        {
            printf("particle: %u: hash: %u\n", i, (uint)*(dGPH + i));
        }
        printf("\n");
#endif

        // check if kernel invocation generated an error
        getLastCudaError("Kernel execution failed");
    }


    void reorderDataAndFindCellStart(uint  *cellStart,
                                     uint  *cellEnd,
                                     float *sortedPos,
                                     float *sortedW,
                                     int   *sortedPhase,
                                     uint  *gridParticleHash,
                                     uint  *gridParticleIndex,
                                     float *oldPos,
                                     uint   numParticles,
                                     uint   numCells)
    {
        uint numThreads, numBlocks;
        computeGridSize(numParticles, 256, numBlocks, numThreads);

        // set all cells to empty
        checkCudaErrors(cudaMemset(cellStart, 0xffffffff, numCells*sizeof(uint)));

        float *dW = getWRawPtr();
        int *dPhase = getPhaseRawPtr();

        thrust::device_ptr<float> d_pos(oldPos);
        thrust::device_ptr<float> d_W(dW);
        thrust::device_ptr<int> d_Phase(dPhase);

//        printf("particles\n");
//        for (int i = 0; i < numParticles; i++)
//        {
//            printf("%.2f %.2f %.2f %.2f, ", (float)*(d_pos+i*4), (float)*(d_pos+i*4+1), (float)*(d_pos+i*4+2), (float)*(d_pos+i*4+3));
//            printf("w: %.2f, phase: %d\n", (float)*(d_W+i), (int)*(d_Phase+i));
//        }
//        printf("\n");
//        float4 *pos4 = (float4*)oldPos;
        float *pos;
        checkCudaErrors(cudaMalloc((void**)&pos, numParticles*4*sizeof(float)));
//        textureVec.resize(4 * numParticles);
//        float *pos = thrust::raw_pointer_cast(textureVec.data());
        checkCudaErrors(cudaMemcpy(pos, oldPos, numParticles*4*sizeof(float), cudaMemcpyDeviceToDevice));

        checkCudaErrors(cudaBindTexture(0, oldPosTex, pos, numParticles*sizeof(float4)));
        checkCudaErrors(cudaBindTexture(0, invMassTex, dW, numParticles*sizeof(float)));
        checkCudaErrors(cudaBindTexture(0, oldPhaseTex, dPhase, numParticles*sizeof(int)));

        uint smemSize = sizeof(uint)*(numThreads+1);
        reorderDataAndFindCellStartD<<< numBlocks, numThreads, smemSize>>>(cellStart,
                                                                           cellEnd,
                                                                           (float4 *) sortedPos,
                                                                           sortedW,
                                                                           sortedPhase,
                                                                           gridParticleHash,
                                                                           gridParticleIndex,
                                                                           (float4 *) oldPos,
                                                                           dW,
                                                                           dPhase,
                                                                           numParticles);
        getLastCudaError("Kernel execution failed: reorderDataAndFindCellStartD");

#ifdef PRINT
        printf("Sorted:\n");
        thrust::device_ptr<uint> dGPH(gridParticleHash);
        thrust::device_ptr<uint> dGPI(gridParticleIndex);
        for (uint i = 0; i < numParticles; i++)
        {
            printf("i: %u: hash: %u\n", i, (uint)*(dGPH + i));
            printf("i: %u: part: %u\n", i, (uint)*(dGPI + i));
        }
        printf("\n");


        printf("Sorted:\n");
        thrust::device_ptr<uint> dstart(cellStart);
        thrust::device_ptr<uint> dend(cellEnd);
        for (uint i = 0; i < 16; i++)
        {
            printf("i: %u: start: %u\n", i, (uint)*(dstart + i));
            printf("i: %u: end: %u\n", i, (uint)*(dend + i));
        }
        printf("\n");
#endif

        checkCudaErrors(cudaUnbindTexture(oldPosTex));
        checkCudaErrors(cudaUnbindTexture(invMassTex));
        checkCudaErrors(cudaUnbindTexture(oldPhaseTex));

        checkCudaErrors(cudaFree(pos));
    }

    void sortParticles(uint *dGridParticleHash, uint *dGridParticleIndex, uint numParticles)
    {
        thrust::sort_by_key(thrust::device_ptr<uint>(dGridParticleHash),
                            thrust::device_ptr<uint>(dGridParticleHash + numParticles),
                            thrust::device_ptr<uint>(dGridParticleIndex));
    }


































    /*****************************************************************************
     *                              PROCESS COLLISIONS
     *****************************************************************************/

    void sortByType(float *dPos, uint numParticles)
    {

    }

    void collideWorld(float *pos, float *sortedPos, uint numParticles, int3 minBounds, int3 maxBounds)
    {
        thrust::device_ptr<float4> d_pos4((float4 *)pos);
        thrust::device_ptr<float4> d_Xstar((float4*)getXstarRawPtr());
        thrust::device_ptr<int> d_phase(getPhaseRawPtr());

        // create random vars for boundary collisions
        checkCudaErrors(curandGenerateUniform(gen, rands, 6));

        // check for boundary collisions and move particles
//        thrust::for_each
//        thrust::transform(d_pos4, d_pos4 + numParticles, d_Xstar, d_pos4, collide_world_functor(rands, minBounds, maxBounds));

        thrust::for_each(
            thrust::make_zip_iterator(thrust::make_tuple(d_pos4, d_Xstar, d_phase)),
            thrust::make_zip_iterator(thrust::make_tuple(d_pos4+numParticles, d_Xstar+numParticles, d_phase+numParticles)),
            collide_world_functor(rands, minBounds, maxBounds));
    }

    void collide(float *particles,
                 float *sortedPos,
                 float *sortedW,
                 int   *sortedPhase,
                 uint  *gridParticleIndex,
                 uint  *cellStart,
                 uint  *cellEnd,
                 uint   numParticles,
                 uint   numCells)
    {
        checkCudaErrors(cudaBindTexture(0, oldPosTex, sortedPos, numParticles*sizeof(float4)));
        checkCudaErrors(cudaBindTexture(0, invMassTex, sortedW, numParticles*sizeof(float)));
        checkCudaErrors(cudaBindTexture(0, oldPhaseTex, sortedPhase, numParticles*sizeof(int)));

        checkCudaErrors(cudaBindTexture(0, cellStartTex, cellStart, numCells*sizeof(uint)));
        checkCudaErrors(cudaBindTexture(0, cellEndTex, cellEnd, numCells*sizeof(uint)));

        // store neighbors
        uint *dNeighbors = thrust::raw_pointer_cast(neighbors.data());
        uint *dNumNeighbors = thrust::raw_pointer_cast(numNeighbors.data());
        float *dXstar = getXstarRawPtr();

        // thread per particle
        uint numThreads, numBlocks;
        computeGridSize(numParticles, 64, numBlocks, numThreads);

        // execute the kernel
        collideD<<< numBlocks, numThreads >>>((float4 *)particles,
                                              (float4 *)dXstar,
                                              (float4 *)sortedPos,
                                              sortedW,
                                              sortedPhase,
                                              gridParticleIndex,
                                              cellStart,
                                              cellEnd,
                                              numParticles,
                                              dNeighbors,
                                              dNumNeighbors);

        // check if kernel invocation generated an error
        getLastCudaError("Kernel execution failed");

        checkCudaErrors(cudaUnbindTexture(oldPosTex));
        checkCudaErrors(cudaUnbindTexture(invMassTex));
        checkCudaErrors(cudaUnbindTexture(oldPhaseTex));

        checkCudaErrors(cudaUnbindTexture(cellStartTex));
        checkCudaErrors(cudaUnbindTexture(cellEndTex));
    }

























    /*****************************************************************************
     *                              UPDATE VELOCITIES
     *****************************************************************************/

    void calcVelocity(float *dpos, float deltaTime, uint numParticles)
    {
        float *dXstar = getXstarRawPtr();
        thrust::device_ptr<float4> d_Xstar((float4*)dXstar);
        thrust::device_ptr<float4> d_pos((float4*)dpos);
        thrust::device_ptr<float4> d_vel((float4*)thrust::raw_pointer_cast(V.data()));


        thrust::transform(d_pos, d_pos + numParticles, d_Xstar, d_vel, subtract_functor(deltaTime));

    }























    /*****************************************************************************
     *                              SOLVE FLUIDS
     *****************************************************************************/
    void solveFluids(float *sortedPos,
                     float *sortedW,
                     int   *sortedPhase,
                     uint  *gridParticleIndex,
                     uint  *cellStart,
                     uint  *cellEnd,
                     float *particles,
                     uint   numParticles,
                     uint   numCells)
    {
        checkCudaErrors(cudaBindTexture(0, oldPosTex, sortedPos, numParticles*sizeof(float4)));
        checkCudaErrors(cudaBindTexture(0, invMassTex, sortedW, numParticles*sizeof(float)));
        checkCudaErrors(cudaBindTexture(0, oldPhaseTex, sortedPhase, numParticles*sizeof(float4)));
        checkCudaErrors(cudaBindTexture(0, cellStartTex, cellStart, numCells*sizeof(uint)));
        checkCudaErrors(cudaBindTexture(0, cellEndTex, cellEnd, numCells*sizeof(uint)));

        // thread per particle
        uint numThreads, numBlocks;
        computeGridSize(numParticles, 256, numBlocks, numThreads);

        float *dLambda = thrust::raw_pointer_cast(lambda.data());
//        float *dDenom = thrust::raw_pointer_cast(denom.data());
        uint *dNeighbors = thrust::raw_pointer_cast(neighbors.data());
        uint *dNumNeighbors = thrust::raw_pointer_cast(numNeighbors.data());
        float *dRos = thrust::raw_pointer_cast(ros.data());

//        printf("ros: %u, numParts: %u\n", (uint)ros.size(), numParticles);

        // execute the kernel
        findLambdasD<<< numBlocks, numThreads >>>(dLambda,
                                                  gridParticleIndex,
                                                  cellStart,
                                                  cellEnd,
                                                  numParticles,
                                                  dNeighbors,
                                                  dNumNeighbors,
                                                  dRos);

        // execute the kernel
        solveFluidsD<<< numBlocks, numThreads >>>(dLambda,
                                                  gridParticleIndex,
                                                  (float4 *) particles,
                                                  numParticles,
                                                  dNeighbors,
                                                  dNumNeighbors,
                                                  dRos);

        // check if kernel invocation generated an error
        getLastCudaError("Kernel execution failed");

        checkCudaErrors(cudaUnbindTexture(oldPosTex));
        checkCudaErrors(cudaUnbindTexture(invMassTex));
        checkCudaErrors(cudaUnbindTexture(oldPhaseTex));
        checkCudaErrors(cudaUnbindTexture(cellStartTex));
        checkCudaErrors(cudaUnbindTexture(cellEndTex));
    }
}
