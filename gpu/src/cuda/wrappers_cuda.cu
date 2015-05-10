
//#define PRINT
//#define JACOBI
//#define GAUSS
//#define SOLVER


#include <cuda_runtime.h>
#include <cuda_gl_interop.h>


#include <cublas_v2.h>
#include <cusparse.h>
#ifdef SOLVER
#include <cusolverDn.h>
#endif

#include "helper_cuda.h"
#include "kernel_impl.cuh"

#include "thrust/device_ptr.h"
#include "thrust/device_vector.h"
#include "thrust/for_each.h"
#include "thrust/iterator/zip_iterator.h"
//#include "thrust/sort.h"
//#include "thrust/reduce.h"
#include "thrust/transform.h"

//#include <stdlib.h>
//#include <string.h>
#include <stdio.h>

#ifdef PRINT
#include "cusp/csr_matrix.h"
#include "cusp/print.h"
//#include "cusp/multiply.h"
//#include "cusp/transpose.h"
#include "cusp/array1d.h"
#endif

typedef typename cusp::array1d_view< thrust::device_vector<uint>::iterator > IndexArrayView;
typedef typename cusp::array1d_view< thrust::device_vector<float>::iterator > ValueArrayView;
typedef cusp::csr_matrix_view<IndexArrayView, IndexArrayView, ValueArrayView> CSRView;

cublasHandle_t cublasHandle;
cusparseHandle_t cusparseHandle;
cusparseMatDescr_t matDescr;

#ifdef SOLVER
cusolverDnHandle_t cusolverHandle;
#endif

thrust::device_vector<uint> Jr; // row offsets for J
thrust::device_vector<uint> Jc; // col indices for J
thrust::device_vector<float> Jv; // values for J

// D contains Diag of A
thrust::device_vector<uint> Dr; // row offsets for D
thrust::device_vector<uint> Dc; // col indices for D
thrust::device_vector<float> Dv; // values for D

thrust::device_vector<uint> occurences;     // number of constraints affecting a particle
thrust::device_vector<float> W;     // vector of inverse masses
thrust::device_vector<float> JT;    // dense matrix of J transpose
thrust::device_vector<float> A;     // Solution to J * W * JT // doubles as T in jacobi iteration
thrust::device_vector<float> B;
thrust::device_vector<float> B2;
thrust::device_vector<float> C;


//thrust::device_vector<float> deltas;
float *d_oldPos;

int workSize;
float *d_work;
int *d_devIpiv;
int *d_devInfo;

////////////////////////// RIGID BODIES ///////////////////////

thrust::device_vector<float> d_rbAngle;
thrust::device_vector<float> d_rbCenterMass;    // dim2
thrust::device_vector<uint> d_rbPartStart;
thrust::device_vector<uint> d_rbRStart;
thrust::device_vector<uint> d_rbSize;

thrust::device_vector<float> d_Rs;              // dim 2


extern "C"
{

    void allocateArray(void **devPtr, size_t size)
    {
        checkCudaErrors(cudaMalloc(devPtr, size));
    }

    void freeArray(void *devPtr)
    {
        checkCudaErrors(cudaFree(devPtr));
    }

    void cudaInit()
    {
        int devID;

        // use command-line specified CUDA device, otherwise use device with highest Gflops/s
        devID = findCudaDevice();

        if (devID < 0)
        {
            printf("No CUDA Capable devices found, exiting...\n");
            exit(EXIT_SUCCESS);
        }

        checkCudaErrors(cublasCreate(&cublasHandle));
        checkCudaErrors(cusparseCreate(&cusparseHandle));

        checkCudaErrors(cusparseCreateMatDescr(&matDescr));
        cusparseSetMatType(matDescr, CUSPARSE_MATRIX_TYPE_GENERAL);
        cusparseSetMatIndexBase(matDescr,CUSPARSE_INDEX_BASE_ZERO);
#ifdef SOLVER
        checkCudaErrors(cusolverDnCreate(&cusolverHandle));
#endif
    }

    void cudaClean()
    {
        // Destroy the handles
        checkCudaErrors(cublasDestroy(cublasHandle));
        checkCudaErrors(cusparseDestroy(cusparseHandle));
#ifdef SOLVER
        checkCudaErrors(cusolverDnDestroy(cusolverHandle));
#endif
    }

    void initDeviceVectors(uint *indices, uint numParticles, uint numConstraints)
    {
        thrust::device_ptr<uint> d_i(indices);

        occurences.resize(numParticles);
        thrust::device_vector<uint> dv_is(numConstraints * 2);
        thrust::copy(d_i, d_i + numConstraints * 2, dv_is.begin());
        thrust::device_vector<uint> d_occurs(dv_is.size(), 1.f);
        thrust::sort(dv_is.begin(), dv_is.end());
        thrust::reduce_by_key(dv_is.begin(), dv_is.end(), d_occurs.begin(), dv_is.begin(), occurences.begin());

        uint p1i, p2i;
        // row offsets and column indices for matrix J
        for (uint i = 0; i < numConstraints; i++)
        {
            Dr.push_back(i);
            Dc.push_back(i);

            Jr.push_back(i*4);

            p1i =  d_i[i];
            p2i =  d_i[numConstraints + i];
            if (p1i < p2i)
            {
                Jc.push_back(p1i * 2);
                Jc.push_back(p1i*2+1);
                Jc.push_back(p2i * 2);
                Jc.push_back(p2i*2+1);
            }
            else
            {
                Jc.push_back(p2i * 2);
                Jc.push_back(p2i*2+1);
                Jc.push_back(p1i * 2);
                Jc.push_back(p1i*2+1);
            }
        }
        Dr.push_back(numConstraints);
        Jr.push_back(numConstraints*4);

        W.resize(numParticles);
        thrust::fill(W.begin(), W.end(), 1.f);

#ifdef PRINT
        thrust::device_ptr<uint> O(occurences.data());
        printf("Occurs:\n");
        for (uint i = 0; i < occurences.size(); i++)
        {
            printf("%u ", (uint)*(O + i));
        }
        printf("\n");

        thrust::device_ptr<float> d_W(W.data());
        printf("W:\n");
        for (uint i = 0; i < W.size(); i++)
        {
            printf("%.2f ", (float)*(d_W + i));
        }
        printf("\n");
#endif

        uint numParticles2 = numParticles * 2;

        Dv.resize(numConstraints);
        Jv.resize(numConstraints * 4);
        JT.resize(numConstraints * numParticles2);
        thrust::fill(JT.begin(), JT.end(), 0.f);
        A.resize(numConstraints * numConstraints);
        B.resize(numConstraints);
        B2.resize(numConstraints);
        C.resize(numConstraints);
//        deltas.resize(numParticles2);

        allocateArray((void **)&d_oldPos, numParticles * 4 * sizeof(float));

#ifdef SOLVER
        float *pA = thrust::raw_pointer_cast(A.data());
        cusolverDnSpotrf_bufferSize(cusolverHandle, CUBLAS_FILL_MODE_LOWER, numConstraints, pA, numConstraints, &workSize);
        allocateArray((void**)&d_work, workSize * sizeof(float));
        allocateArray((void **)&d_devIpiv, numConstraints * sizeof(int));
        allocateArray((void**)&d_devInfo, sizeof(int));
#endif
    }

    void freeDeviceVectors()
    {
        // set size to zero
        Jr.clear();
        Jc.clear();
        Jv.clear();

        Dr.clear();
        Dc.clear();
        Dv.clear();

        occurences.clear();
        W.clear();
        JT.clear();
        A.clear();
        B.clear();
        B2.clear();
        C.clear();
//        deltas.clear();

        // free memory
        Jr.shrink_to_fit();
        Jc.shrink_to_fit();
        Jv.shrink_to_fit();

        Dr.shrink_to_fit();
        Dc.shrink_to_fit();
        Dv.shrink_to_fit();

        occurences.shrink_to_fit();
        W.shrink_to_fit();
        JT.shrink_to_fit();
        A.shrink_to_fit();
        B.shrink_to_fit();
        B2.shrink_to_fit();
        C.shrink_to_fit();
//        deltas.shrink_to_fit();
        freeArray(d_oldPos);

#ifdef SOLVER
        freeArray(d_work);
        freeArray(d_devIpiv);
        freeArray(d_devInfo);
#endif
        d_rbAngle.clear();
        d_rbCenterMass.clear();
        d_rbPartStart.clear();
        d_rbRStart.clear();
        d_rbSize.clear();

        d_Rs.clear();

        d_rbAngle.shrink_to_fit();
        d_rbCenterMass.shrink_to_fit();
        d_rbPartStart.shrink_to_fit();
        d_rbRStart.shrink_to_fit();
        d_rbSize.shrink_to_fit();

        d_Rs.shrink_to_fit();

    }

    void copyArrayToDevice(void *device, const void *host, int offset, int size)
    {
        checkCudaErrors(cudaMemcpy((char *) device + offset, host, size, cudaMemcpyHostToDevice));
    }

    void registerGLBufferObject(unsigned int vbo, struct cudaGraphicsResource **cuda_vbo_resource)
    {
        checkCudaErrors(cudaGraphicsGLRegisterBuffer(cuda_vbo_resource, vbo,
                                                     cudaGraphicsMapFlagsNone));
    }

    void unregisterGLBufferObject(struct cudaGraphicsResource *cuda_vbo_resource)
    {
        checkCudaErrors(cudaGraphicsUnregisterResource(cuda_vbo_resource));
    }

    void *mapGLBufferObject(struct cudaGraphicsResource **cuda_vbo_resource)
    {
        void *ptr;
        checkCudaErrors(cudaGraphicsMapResources(1, cuda_vbo_resource, 0));
        size_t num_bytes;
        checkCudaErrors(cudaGraphicsResourceGetMappedPointer((void **)&ptr, &num_bytes,
                                                             *cuda_vbo_resource));
        return ptr;
    }

    void unmapGLBufferObject(struct cudaGraphicsResource *cuda_vbo_resource)
    {
        checkCudaErrors(cudaGraphicsUnmapResources(1, &cuda_vbo_resource, 0));
    }

    void copyArrayFromDevice(void *host, const void *device, int size)
    {
        checkCudaErrors(cudaMemcpy(host, device, size, cudaMemcpyDeviceToHost));
    }

    void setParameters(SimParams *hostParams)
    {
        // copy parameters to constant memory
        checkCudaErrors(cudaMemcpyToSymbol(params, hostParams, sizeof(SimParams)));
    }

    //Round a / b to nearest higher integer value
    uint iDivUp(uint a, uint b)
    {
        return (a % b != 0) ? (a / b + 1) : (a / b);
    }

    // compute grid and thread block size for a given number of elements
    void computeGridSize(uint n, uint blockSize, uint &numBlocks, uint &numThreads)
    {
        numThreads = min(blockSize, n);
        numBlocks = iDivUp(n, numThreads);
    }

    void integrateSystem(float *pos, float *vel, float deltaTime, uint numParticles)
    {
        thrust::device_ptr<float4> d_pos4((float4 *)pos);
        thrust::device_ptr<float4> d_vel4((float4 *)vel);

        thrust::device_ptr<float4> d_old4((float4 *)d_oldPos);
        thrust::copy(d_pos4, d_pos4 + numParticles, d_old4);

        thrust::for_each(
            thrust::make_zip_iterator(thrust::make_tuple(d_pos4, d_vel4)),
            thrust::make_zip_iterator(thrust::make_tuple(d_pos4+numParticles, d_vel4+numParticles)),
            integrate_functor(deltaTime));

#ifdef PRINT
        for (uint i = 0; i < numParticles; i++)
        {
            printf("Integrated Pos : (%.2f, %.2f, %.2f, %.2f)\n",
                   ((float4)*(d_pos4+i)).x,((float4)*(d_pos4+i)).y,
                   ((float4)*(d_pos4+i)).z,((float4)*(d_pos4+i)).w);
        }
#endif
    }

    void calcHash(uint *gridParticleHash, uint *gridParticleIndex, float *pos, int numParticles)
    {
        uint numThreads, numBlocks;
        computeGridSize(numParticles, 256, numBlocks, numThreads);

        // execute the kernel
        calcHashD<<< numBlocks, numThreads >>>(gridParticleHash, gridParticleIndex, (float4 *) pos, numParticles);

        // check if kernel invocation generated an error
        getLastCudaError("Kernel execution failed");
    }


    void reorderDataAndFindCellStart(uint  *cellStart,
                                     uint  *cellEnd,
                                     float *sortedPos,
                                     float *sortedVel,
                                     uint  *gridParticleHash,
                                     uint  *gridParticleIndex,
                                     float *oldPos,
                                     float *oldVel,
                                     uint   numParticles,
                                     uint   numCells)
    {
        uint numThreads, numBlocks;
        computeGridSize(numParticles, 256, numBlocks, numThreads);

        // set all cells to empty
        checkCudaErrors(cudaMemset(cellStart, 0xffffffff, numCells*sizeof(uint)));

        checkCudaErrors(cudaBindTexture(0, oldPosTex, oldPos, numParticles*sizeof(float4)));
        checkCudaErrors(cudaBindTexture(0, oldVelTex, oldVel, numParticles*sizeof(float4)));

        uint smemSize = sizeof(uint)*(numThreads+1);
        reorderDataAndFindCellStartD<<< numBlocks, numThreads, smemSize>>>(
            cellStart,
            cellEnd,
            (float4 *) sortedPos,
            (float4 *) sortedVel,
            gridParticleHash,
            gridParticleIndex,
            (float4 *) oldPos,
            (float4 *) oldVel,
            numParticles);
        getLastCudaError("Kernel execution failed: reorderDataAndFindCellStartD");

        checkCudaErrors(cudaUnbindTexture(oldPosTex));
        checkCudaErrors(cudaUnbindTexture(oldVelTex));
    }

    void collide(float *newVel,
                 float *sortedPos,
                 float *sortedVel,
                 uint  *gridParticleIndex,
                 uint  *cellStart,
                 uint  *cellEnd,
                 uint   numParticles,
                 uint   numCells)
    {
        checkCudaErrors(cudaBindTexture(0, oldPosTex, sortedPos, numParticles*sizeof(float4)));
        checkCudaErrors(cudaBindTexture(0, oldVelTex, sortedVel, numParticles*sizeof(float4)));
        checkCudaErrors(cudaBindTexture(0, cellStartTex, cellStart, numCells*sizeof(uint)));
        checkCudaErrors(cudaBindTexture(0, cellEndTex, cellEnd, numCells*sizeof(uint)));

        // thread per particle
        uint numThreads, numBlocks;
        computeGridSize(numParticles, 64, numBlocks, numThreads);

        // execute the kernel
        collideD<<< numBlocks, numThreads >>>((float4 *)newVel,
                                              (float4 *)sortedPos,
                                              (float4 *)sortedVel,
                                              gridParticleIndex,
                                              cellStart,
                                              cellEnd,
                                              numParticles);

        // check if kernel invocation generated an error
        getLastCudaError("Kernel execution failed");

        checkCudaErrors(cudaUnbindTexture(oldPosTex));
        checkCudaErrors(cudaUnbindTexture(oldVelTex));
        checkCudaErrors(cudaUnbindTexture(cellStartTex));
        checkCudaErrors(cudaUnbindTexture(cellEndTex));
    }

    void sortParticles(uint *dGridParticleHash, uint *dGridParticleIndex, uint numParticles)
    {
        thrust::sort_by_key(thrust::device_ptr<uint>(dGridParticleHash),
                            thrust::device_ptr<uint>(dGridParticleHash + numParticles),
                            thrust::device_ptr<uint>(dGridParticleIndex));
    }


//    void solvePointConstraints(uint *indices, float *points, float *particles, uint numConstraints)
//    {
//        thrust::device_ptr<uint> d_indices(indices);
//        thrust::device_ptr<float4> d_points((float4 *) points);

//        thrust::for_each(
//            thrust::make_zip_iterator(thrust::make_tuple(d_indices, d_points)),
//            thrust::make_zip_iterator(thrust::make_tuple(d_indices+numConstraints, d_points+numConstraints)),
//            point_constraint_functor((float4 *)particles));
//    }

//    void iterativeSolveDistanceConstraints(uint *indices, float *distance, float *particles, uint numConstraints)
//    {
//        // copy indices to device vector (this gets manipulated)
//        thrust::device_vector<uint> d_indices(indices, indices + numConstraints*2);

//        // pointers to the locations of the first and second indices
//        thrust::device_ptr<uint> d_index1 = d_indices.data();
//        thrust::device_ptr<uint> d_index2 = d_indices.data() + numConstraints;

//        // pointer to the distance values
//        thrust::device_ptr<float> d_distances(distance);

//        // device vector of delta values (to accumulate then sum)
//        thrust::device_vector<float4> d_deltas_vec(numConstraints*2);
//        thrust::device_ptr<float4> d_deltas_ptr = d_deltas_vec.data();

//        // calculate the delta values for each point
//        thrust::for_each(
//            thrust::make_zip_iterator(thrust::make_tuple(d_index1, d_index2, d_distances, d_deltas_ptr, d_deltas_ptr+ numConstraints)),
//            thrust::make_zip_iterator(thrust::make_tuple(d_index1+numConstraints, d_index2+numConstraints, d_distances+numConstraints, d_deltas_ptr+numConstraints, d_deltas_ptr+numConstraints*2)),
//            delta_computing_functor((float4 *) particles));

//        thrust::sort_by_key(d_index1, d_index1+numConstraints*2, d_deltas_ptr);

//        thrust::pair<thrust::device_ptr<uint>,thrust::device_ptr<float4>> new_end;
//        new_end = thrust::reduce_by_key(d_index1, d_index1+numConstraints*2, d_deltas_ptr, d_index1, d_deltas_ptr);

//        thrust::for_each(
//            thrust::make_zip_iterator(thrust::make_tuple(d_index1, d_deltas_ptr)),
//            thrust::make_zip_iterator(thrust::make_tuple(new_end.first, new_end.second)),
//            distance_constraint_functor((float4 *)particles));
//    }

    void calcVelocity(float *hpos, float *dpos, float *vel, float deltaTime, uint numParticles)
    {
//        thrust::host_vector<float> h_origPos(hpos, hpos + numParticles*4);
//        thrust::device_vector<float> d_origPos = h_origPos;
        thrust::device_ptr<float4> d_origPos((float4*)d_oldPos);
        thrust::device_ptr<float4> d_solvedPos((float4*)dpos);
        thrust::device_ptr<float4> d_vel((float4*)vel);

        thrust::transform(d_origPos, d_origPos + numParticles * 4, d_solvedPos, d_vel, subtract_functor(deltaTime));
    }


//    void buildJAndB(uint *indices, float *distance, float *particles, float *jay_, float *b, uint numParticles, uint numConstraints)
//    {
////        uint2 sizeJ = make_uint2(numParticles * 2, numConstraints); // rows, cols
////        // raw pointer to device memory
//////        float *_J, *_b;
//////        allocateArray((void **)&_J, sizeJ.x * sizeJ.y * sizeof(float));
//////        allocateArray((void **)&_b, numConstraints * sizeof(float));

////        // wrap raw pointer with a device_ptr
////        thrust::device_ptr<float> d_J(J);
//        thrust::device_ptr<float> d_b(b);

////        // use device_ptr in thrust algorithms
////        thrust::fill(d_J, d_J + sizeJ.x * sizeJ.y, 0.f);
////        thrust::fill(d_b, d_b + numConstraints, 0.f);

//        uint numPartsExp = numParticles * 2; // * 3 for 3D
////        cusp::print(W);

////        J.column_indices(1, 1) = 0;




//#ifdef PRINT
//        thrust::device_ptr<float4> d_particles((float4*)particles);
//        for (uint i = 0; i < numParticles; i++)
//        {
//            printf("particles before : (%.2f, %.2f, %.2f, %.2f)\n",
//                   ((float4)*(d_particles+i)).x,((float4)*(d_particles+i)).y,
//                   ((float4)*(d_particles+i)).z,((float4)*(d_particles+i)).w);
//        }

//        printf("MATRIX J (before):\n");
//        for (uint r = 0; r < sizeJ.x; r++)
//        {
//            for (uint c = 0; c < sizeJ.y; c++)
//            {
//                uint index = c * sizeJ.x + r;
//                printf("%.2f ", (float)*(d_J+index));
//            }
//            printf("\n");
//        }

//        printf("B (before):\n");
//        for (uint i = 0; i < numConstraints * 4; i++)
//        {
//            printf("%.2f, ", (float)*(d_b+i));
//        }
//        printf("\b \n");
//#endif

////        // copy indices to device vector (this gets manipulated)
////        thrust::device_vector<uint> d_indices(indices, indices + numConstraints*2);

////        // pointers to the locations of the first and second indices
////        thrust::device_ptr<uint> d_index1 = d_indices.data();
////        thrust::device_ptr<uint> d_index2 = d_indices.data() + numConstraints;

//        // pointer to the indices
//        thrust::device_ptr<uint> d_index1(indices);
//        thrust::device_ptr<uint> d_index2(indices + numConstraints);
//        // pointer to the distance values
//        thrust::device_ptr<float> d_distances(distance);

//        // keeps track of the index
//        thrust::counting_iterator<uint> first(0);
//        thrust::counting_iterator<uint> last = first + numConstraints;

////        J.row_offsets[0] = 0;
////        J.row_offsets[1] = 4;
////        J.row_offsets[2] = 7;

////        J.column_indices[0] = 1; J.values[0] = 1.f;
////        J.column_indices[1] = 3; J.values[1] = 2.f;
////        J.column_indices[2] = 4; J.values[2] = 3.f;
////        J.column_indices[3] = 5; J.values[3] = 4.f;
////        J.column_indices[4] = 0; J.values[4] = 5.f;
////        J.column_indices[5] = 1; J.values[5] = 6.f;
////        J.column_indices[6] = 2; J.values[6] = 7.f;
////        J.column_indices[7] = 4; J.values[7] = 8.f;

//        cusp::csr_matrix<uint, float, cusp::device_memory> J(numConstraints, numPartsExp, numConstraints * 4);

//        uint *row_off = thrust::raw_pointer_cast(&J.row_offsets[0]);
//        uint *col_ind = thrust::raw_pointer_cast(&J.column_indices[0]);
//        float *val = thrust::raw_pointer_cast(&J.values[0]);

////        CSRView j(2, 6, 8, );



//        // fill J and b
////        thrust::for_each(
////                    thrust::make_zip_iterator(thrust::make_tuple(first, d_index1, d_index2, d_distances, d_b)),
////                    thrust::make_zip_iterator(thrust::make_tuple(last, d_index1+numConstraints, d_index2+numConstraints, d_distances+numConstraints, d_b+numConstraints)),
////                    gradient_functor((float4 *) particles, row_off, col_ind, val, numConstraints));

//        J.row_offsets[numConstraints] = numConstraints * 4 - 1;

//#ifdef PRINT
//        printf("MATRIX J (after):\n");
//        for (uint r = 0; r < sizeJ.x; r++)
//        {
//            for (uint c = 0; c < sizeJ.y; c++)
//            {
//                uint index = c * sizeJ.x + r;
//                printf("%.2f ", (float)*(d_J+index));
//            }
//            printf("\n");
//        }

//        printf("B (after):\n");
//        for (uint i = 0; i < numConstraints; i++)
//        {
//            printf("%.2f, ", (float)*(d_b+i));
//        }
//        printf("\b \n");
//#endif

//        cusp::csr_matrix<uint, float, cusp::device_memory> JW(numConstraints, numPartsExp, numConstraints * 4);
//        cusp::csr_matrix<uint, float, cusp::device_memory> JT(numPartsExp, numConstraints, numConstraints * 4);
//        cusp::csr_matrix<uint, float, cusp::device_memory> A(numConstraints, numConstraints, numConstraints * numConstraints);

////        cusp::multiply(J, W, JW);
////        cusp::transpose(J, JT);
////        cusp::multiply(JW,JT,A);

////        printf("J:\n");
////        cusp::print(J);
////        printf("JW:\n");
////        cusp::print(JW);
////        printf("JT:\n");
////        cusp::print(J);
////        printf("A:\n");
////        cusp::print(A);

////        freeArray(_J);
////        freeArray(_b);
//    }


//    void makeA(float *J, float *A, uint numParticles, uint numConstraints)
//    {
//        float alpha = 1.f, beta = 0.f;
//        checkCudaErrors(cublasSgemm(cublasHandle, CUBLAS_OP_T, CUBLAS_OP_N, numConstraints, numConstraints, numParticles * 2, &alpha, J, numParticles * 2, J, numParticles * 2, &beta, A, numConstraints));
//#ifdef PRINT
//        thrust::device_ptr<float> d_A(A);
//        printf("MATRIX A:\n");
//        for (uint r = 0; r < numConstraints; r++)
//        {
//            for (uint c = 0; c < numConstraints; c++)
//            {
//                uint index = c * numConstraints + r;
//                printf("%.2f ", (float)*(d_A+index));
//            }
//            printf("\n");
//        }
//#endif
//    }

    void solveAxb(uint *indices, float *distances, float *particles, uint numParticles, uint numConstraints, int gaussIters, float omega)
    {
        uint numParticles2 = numParticles * 2;

        float *val = thrust::raw_pointer_cast(Jv.data());

        // pointer to the indices
        thrust::device_ptr<uint> d_index1(indices);
        thrust::device_ptr<uint> d_index2(indices + numConstraints);
        // pointer to the distance values
        thrust::device_ptr<float> d_distances(distances);

#ifdef PRINT
        printf("dis: ");
        for (int i = 0; i < numConstraints*2; i++)
            printf("%u ", (uint)*(d_index1 + i));
        printf("\n");

        printf("ds: ");
        for (int i = 0; i < numConstraints; i++)
            printf("%u ", (uint)*(d_distances + i));
        printf("\n");
#endif

        // keeps track of the index
        thrust::counting_iterator<uint> first(0);
        thrust::counting_iterator<uint> last = first + numConstraints;

        thrust::device_ptr<float> d_B(B.data());
        float *pJT = thrust::raw_pointer_cast(JT.data());

        // fill J and b
        thrust::for_each(
                    thrust::make_zip_iterator(thrust::make_tuple(first, d_index1, d_index2, d_distances, d_B)),
                    thrust::make_zip_iterator(thrust::make_tuple(last, d_index1+numConstraints, d_index2+numConstraints, d_distances+numConstraints, d_B+numConstraints)),
                    gradient_functor((float4 *) particles, val, pJT, numParticles, numConstraints));

        float alpha = 1.f;
        float beta = 0.f;

        // raw pointers to data from vectors
        float *pJv = thrust::raw_pointer_cast(Jv.data());
        int *pJr = (int*)thrust::raw_pointer_cast(Jr.data());
        int *pJc = (int*)thrust::raw_pointer_cast(Jc.data());

        float *pDv = thrust::raw_pointer_cast(Dv.data());
        int *pDr = (int*)thrust::raw_pointer_cast(Dr.data());
        int *pDc = (int*)thrust::raw_pointer_cast(Dc.data());

        float *pA = thrust::raw_pointer_cast(A.data());
        float *pB = thrust::raw_pointer_cast(B.data());
        float *pB2 = thrust::raw_pointer_cast(B2.data());
        float *pC = thrust::raw_pointer_cast(C.data());

        // A = J * JT
        cusparseScsrmm(
                    cusparseHandle,
                    CUSPARSE_OPERATION_NON_TRANSPOSE,
                    numConstraints,      // m: rows of J
                    numConstraints,     // n: cols of JT and A
                    numParticles2,      // k: cols of J
                    numConstraints * 4, // nnz: num non-zero elements of J
                    &alpha,             // needs to be 1.f
                    matDescr,           // matrix descriptor for J
                    pJv,                // non zeros values of J
                    pJr,                // row offsets of J
                    pJc,                // col indices of J
                    pJT,                // dense matrix JT
                    numParticles2,      // leading dimension of JT
                    &beta,              // needs to be 0.f
                    pA,                 // matrix A: answer is stored here
                    numConstraints);    // leading dimension of A

#ifdef GAUSS
        float minimumError = 0.0001f;
        uint maximumIterations = gaussIters;

        uint numThreads, numBlocks;
        computeGridSize(numConstraints, 256, numBlocks, numThreads);

        thrust::fill(C.begin(), C.end(), 0.f);
//        thrust::fill(deltas.begin(), deltas.end(), 0.f);

//        float *pDeltas = thrust::raw_pointer_cast(deltas.data());
        uint *pOccurences = thrust::raw_pointer_cast(occurences.data());

        // execute the kernel
        constraintCentricSolveD<<< numBlocks, numThreads >>>(pA,
                                                             pB,
                                                             pJT,
                                                             pC,
                                                             particles,
                                                             indices,
                                                             pOccurences,
                                                             numParticles,
                                                             numConstraints,
                                                             minimumError,
                                                             maximumIterations,
                                                             omega);

#endif
#ifdef JACOBI

//        A[0] = 2;
//        A[1] = 5;
//        A[2] = 1;
//        A[3] = 7;

//        B[0] = 11;
//        B[1] = 13;

        // set D as diagonal from A and set A diag to zeros
        thrust::for_each(
                    thrust::make_zip_iterator(thrust::make_tuple(first, Dv.begin())),
                    thrust::make_zip_iterator(thrust::make_tuple(last, Dv.end())),
                    diag_extraction_functor(pA, numConstraints));

        // T = -D * A (A doubles as T after this computation)
        alpha = -1.f;
        cusparseScsrmm(
                    cusparseHandle,
                    CUSPARSE_OPERATION_NON_TRANSPOSE,
                    numConstraints,      // m: rows of D
                    numConstraints,     // n: cols of D and T (stored in A)
                    numConstraints,      // k: cols of D
                    numConstraints, // nnz: num non-zero elements of D
                    &alpha,             // -1.f
                    matDescr,           // matrix descriptor for D
                    pDv,                // non zeros values of D
                    pDr,                // row offsets of D
                    pDc,                // col indices of D
                    pA,                // dense matrix A
                    numConstraints,      // leading dimension of A
                    &beta,              // needs to be 0.f
                    pA,                 // matrix A: answer is stored here
                    numConstraints);    // leading dimension of A

        // c = D * b
        alpha = 1.f;
        cusparseScsrmv(cusparseHandle,
                       CUSPARSE_OPERATION_NON_TRANSPOSE,
                       numConstraints,
                       numConstraints,
                       numConstraints,
                       &alpha,
                       matDescr,
                       pDv,
                       pDr,
                       pDc,
                       pB,
                       &beta,
                       pC);

        thrust::fill(B.begin(), B.end(), 1.f);
        beta = 1.f;

        for (int i = 0; i < 100; i++)
        {
            thrust::copy(C.begin(), C.end(), B2.begin());
            cublasSgemv(cublasHandle,
                        CUBLAS_OP_N,
                        numConstraints,
                        numConstraints,
                        &alpha,
                        pA,
                        numConstraints,
                        pB,
                        1,
                        &beta,
                        pB2,
                        1);

            thrust::copy(C.begin(), C.end(), B.begin());
            cublasSgemv(cublasHandle,
                        CUBLAS_OP_N,
                        numConstraints,
                        numConstraints,
                        &alpha,
                        pA,
                        numConstraints,
                        pB2,
                        1,
                        &beta,
                        pB,
                        1);
        }
#endif
#ifdef SOLVER
        cusolverDnSpotrf(
                    cusolverHandle,
                    CUBLAS_FILL_MODE_LOWER,
                    numConstraints,
                    pA,
                    numConstraints,
                    d_work,
                    workSize,
                    d_devInfo);

        // overwrite b
        cusolverDnSpotrs(
                    cusolverHandle,
                    CUBLAS_FILL_MODE_LOWER,
                    numConstraints,
                    1,
                    pA,
                    numConstraints,
                    pB,
                    numConstraints,
                    d_devInfo);
#endif
#ifndef GAUSS
        // update particles
        thrust::for_each(
                    thrust::make_zip_iterator(thrust::make_tuple(first, W.begin(), occurences.begin())),
                    thrust::make_zip_iterator(thrust::make_tuple(last, W.end(), occurences.end())),
                    col_reduction((float4 *) particles, pB, pJT, numParticles2, numConstraints));
#endif




#ifdef PRINT

        thrust::device_ptr<float> d_deltas(deltas.data());
        printf("Deltas:\n");
        for (uint c = 0; c < numParticles2; c++)
        {
            printf("%.2f ", (float)*(d_deltas+c));
        }
        printf("\n");

        CSRView J(numConstraints, numParticles2, numConstraints * 4,
                cusp::make_array1d_view(Jr.begin(), Jr.end()),
                cusp::make_array1d_view(Jc.begin(), Jc.end()),
                cusp::make_array1d_view(Jv.begin(), Jv.end()));

        printf("J:\n");
        cusp::print(J);

        thrust::device_ptr<float> d_JT(JT.data());
        printf("MATRIX JT:\n");
        for (uint r = 0; r < numParticles2; r++)
        {
            for (uint c = 0; c < numConstraints; c++)
            {
                uint index = c * numParticles2 + r;
                printf("%.2f ", (float)*(d_JT+index));
            }
            printf("\n");
        }
        printf("\n");

        thrust::device_ptr<float> d_A(A.data());
        printf("MATRIX A:\n");
        for (uint r = 0; r < numConstraints; r++)
        {
            for (uint c = 0; c < numConstraints; c++)
            {
                uint index = c * numConstraints + r;
                printf("%.2f ", (float)*(d_A+index));
            }
            printf("\n");
        }
        printf("\n");

        CSRView D(numConstraints, numConstraints, numConstraints,
                cusp::make_array1d_view(Dr.begin(), Dr.end()),
                cusp::make_array1d_view(Dc.begin(), Dc.end()),
                cusp::make_array1d_view(Dv.begin(), Dv.end()));
        printf("MATRIX D:\n");
        cusp::print(D);

        printf("Vector B:\n");
        for (uint c = 0; c < numConstraints; c++)
        {
            printf("%.2f ", (float)*(d_B+c));
        }
        printf("\n");

        thrust::device_ptr<float> d_B2(B2.data());
        printf("Vector B2:\n");
        for (uint c = 0; c < numConstraints; c++)
        {
            printf("%.2f ", (float)*(d_B2+c));
        }
        printf("\n");

        thrust::device_ptr<float> d_C(C.data());
        printf("Vector C:\n");
        for (uint c = 0; c < numConstraints; c++)
        {
            printf("%.2f ", (float)*(d_C+c));
        }
        printf("\n");

        thrust::device_ptr<float4> d_pos4((float4 *)particles);
        for (uint i = 0; i < numParticles; i++)
        {
            printf("Pos : (%.2f, %.2f, %.2f, %.2f)\n",
                   ((float4)*(d_pos4+i)).x,((float4)*(d_pos4+i)).y,
                   ((float4)*(d_pos4+i)).z,((float4)*(d_pos4+i)).w);
        }

#endif
    }


    void addRigidBody(float *hpos, int width, int height, float2 center, float particleRadius)
    {
        float diameter = particleRadius * 2;
        float startX = center.x - (width * .5f) * diameter + particleRadius;
        float startY = center.y - (height * .5f) * diameter + particleRadius;

        float2 posSum = make_float2(0.f);
        int numParticles = 0;

        d_rbAngle.push_back(0);
        d_rbRStart.push_back(d_Rs.size());
        d_rbPartStart.push_back(0);

        float2 pos;
        int index = 0;
        for (int i = 0; i < width; i++)
        {
            for (int j = 0; j < height; j++)
            {
                pos = make_float2(startX + i * diameter, startY + j * diameter);
                hpos[index]   = pos.x;
                hpos[index+1] = pos.y;
                hpos[index+2] = 0.f;
                hpos[index+3] = 1.f;

                d_Rs.push_back(pos.x - center.x);
                d_Rs.push_back(pos.y - center.y);

                posSum += pos;
                numParticles++;
                index += 4;
            }
        }

        d_rbCenterMass.push_back(posSum.x / numParticles);
        d_rbCenterMass.push_back(posSum.y / numParticles);
        d_rbSize.push_back(numParticles);

    }

    void solveRigidBodies()
    {
        float *dRs = thrust::raw_pointer_cast(d_Rs.data());

        thrust::for_each(
                    thrust::make_zip_iterator(thrust::make_tuple(d_rbAngle.begin(), d_rbCenterMass.begin(), d_rbPartStart.begin(), d_rbRStart.begin(), d_rbSize.begin())),
                    thrust::make_zip_iterator(thrust::make_tuple(d_rbAngle.end(), d_rbCenterMass.end(), d_rbPartStart.end(), d_rbRStart.end(), d_rbSize.end())),
                    rigid_body_functor((float2*) dRs));
    }

//    thrust::device_vector<float> d_rbAngle;
//    thrust::device_vector<float> d_rbCenterMass;    // dim2
//    thrust::device_vector<uint> d_rbPartStart;
//    thrust::device_vector<uint> d_rbRStart;
//    thrust::device_vector<uint> d_rbSize;

//    thrust::device_vector<float> d_Rs;              // dim 2
}
