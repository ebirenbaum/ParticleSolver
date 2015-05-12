
//#define PRINT


#include <cuda_runtime.h>
#include <cuda_gl_interop.h>

#include <thrust/device_ptr.h>
#include <thrust/device_vector.h>
#include <thrust/for_each.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/sort.h>
#include <thrust/reduce.h>
#include <thrust/transform.h>

#include <stdio.h>


//#include <cublas_v2.h>
//#include <cusparse.h>

#include "helper_cuda.h"
#include "solver_kernel.cuh"
#include "util.cuh"
#include "shared_variables.cuh"


//cublasHandle_t cublasHandle;
//cusparseHandle_t cusparseHandle;
//cusparseMatDescr_t matDescr;

thrust::device_vector<uint> distsI;
thrust::device_vector<float> dists;

thrust::device_vector<uint> pointsI;
thrust::device_vector<float> points;

thrust::device_vector<uint> sortedI;
thrust::device_vector<float> deltas;

thrust::device_vector<uint> occurences;     // number of constraints affecting a particle

extern "C"
{

//    void initHandles()
//    {

////        checkCudaErrors(cublasCreate(&cublasHandle));
////        checkCudaErrors(cusparseCreate(&cusparseHandle));

////        checkCudaErrors(cusparseCreateMatDescr(&matDescr));
////        cusparseSetMatType(matDescr, CUSPARSE_MATRIX_TYPE_GENERAL);
////        cusparseSetMatIndexBase(matDescr,CUSPARSE_INDEX_BASE_ZERO);

//    }

//    void destroyHandles()
//    {
////        checkCudaErrors(cublasDestroy(cublasHandle));
////        checkCudaErrors(cusparseDestroy(cusparseHandle));
//    }

    void appendSolverParticle(uint numParticles)
    {
        uint sizeO = occurences.size();
        occurences.resize(sizeO + numParticles);
        uint *dOcc = thrust::raw_pointer_cast(occurences.data());
        checkCudaErrors(cudaMemset(dOcc + sizeO, 0, numParticles * sizeof(uint)));
    }

    void updateOccurences(uint *index, uint num)
    {
        thrust::device_vector<uint> dvIndex(index, index + num);
        thrust::device_vector<uint> dvSorted(num);
        thrust::device_vector<uint> dvOnes(num, 1);

        thrust::device_ptr<uint> d_Index(dvIndex.data());
        thrust::device_ptr<uint> d_Ones(dvOnes.data());
        thrust::device_ptr<uint> d_Sorted(dvSorted.data());

        thrust::sort(dvIndex.begin(), dvIndex.end());

//        printf("sorted: %u\n", dvIndex.size());
//        for (uint i = 0; i < dvIndex.size(); i++)
//        {
//            printf("%u\n", (uint)*(d_Index + i));
//        }

        thrust::pair<thrust::device_ptr<uint>,thrust::device_ptr<uint> > new_end;
        new_end = thrust::reduce_by_key(d_Index, d_Index + num, d_Ones, d_Sorted, d_Ones);

//        printf("reduced: %u\n", dvIndex.size());
//        for (uint i = 0; i < dvIndex.size(); i++)
//        {
//            printf("%u\n", (uint)*(d_Index + i));
//        }


        uint *dOcc = thrust::raw_pointer_cast(occurences.data());

        thrust::for_each(
            thrust::make_zip_iterator(thrust::make_tuple(d_Sorted, d_Ones)),
            thrust::make_zip_iterator(thrust::make_tuple(new_end.first, new_end.second)),
            occurence_functor(dOcc));
    }

    void addPointConstraint(uint *index, float *point, uint numConstraints)
    {
        uint sizeP = points.size();
        uint sizeI = pointsI.size();

        points.resize(sizeP + 3 * numConstraints);
        pointsI.resize(sizeI + numConstraints);

        float *dPoints = thrust::raw_pointer_cast(points.data());
        uint *dPointsI = thrust::raw_pointer_cast(pointsI.data());

        copyArrayToDevice(dPoints + sizeP, point, 0, 3 * numConstraints * sizeof(float));
        copyArrayToDevice(dPointsI + sizeI, index, 0, numConstraints * sizeof(uint));

        updateOccurences(index, numConstraints);
    }

    void addDistanceConstraint(uint *index, float *distance, uint numConstraints)
    {
        uint sizeD = dists.size();
        uint sizeI = distsI.size();

        dists.resize(sizeD + numConstraints);
        distsI.resize(sizeI + 2 * numConstraints);

        float *dDists = thrust::raw_pointer_cast(dists.data());
        uint *dDistsI = thrust::raw_pointer_cast(distsI.data());

        copyArrayToDevice(dDists + sizeD, distance, 0, numConstraints * sizeof(float));
        copyArrayToDevice(dDistsI + sizeI, index, 0, 2 * numConstraints * sizeof(uint));

//        thrust::device_ptr<uint> dOcc(occurences.data());
//        printf("before: \n");
//        for (uint i = 0; i < occurences.size(); i++)
//        {
//            printf("%u\n", (uint)*(dOcc + i));
//        }

        sortedI.resize(distsI.size());
        deltas.resize(8 * dists.size());

        updateOccurences(index, 2 * numConstraints);

//        printf("after: \n");
//        for (uint i = 0; i < occurences.size(); i++)
//        {
//            printf("%u\n", (uint)*(dOcc + i));
//        }
    }

    void freeSolverVectors()
    {
        // set size to zero
        distsI.clear();
        dists.clear();
        pointsI.clear();
        points.clear();
        sortedI.clear();
        deltas.clear();
        occurences.clear();

        // free memory
        distsI.shrink_to_fit();
        dists.shrink_to_fit();
        pointsI.shrink_to_fit();
        points.shrink_to_fit();
        sortedI.shrink_to_fit();
        deltas.shrink_to_fit();
        occurences.shrink_to_fit();

    }

    void solvePointConstraints(float *particles)
    {
        uint numConstraints = pointsI.size();

        if (numConstraints == 0)
            return;

        thrust::device_ptr<uint> d_indices(pointsI.data());
        thrust::device_ptr<float3> d_points((float3*) thrust::raw_pointer_cast(points.data()));

        thrust::for_each(
            thrust::make_zip_iterator(thrust::make_tuple(d_indices, d_points)),
            thrust::make_zip_iterator(thrust::make_tuple(d_indices+numConstraints, d_points+numConstraints)),
            point_constraint_functor((float4 *)particles));
    }

    void solveDistanceConstraints(float *particles)
    {
        uint numConstraints = dists.size();

        if (numConstraints == 0)
            return;

        thrust::device_ptr<float4> d_pos4((float4*)particles);

        thrust::device_ptr<uint2> d_indices((uint2*)thrust::raw_pointer_cast(distsI.data()));
        thrust::device_ptr<float> d_dists(dists.data());
        thrust::device_ptr<uint> d_sortedI1(sortedI.data());
        thrust::device_ptr<uint> d_sortedI2 = d_sortedI1 + numConstraints;
        thrust::device_ptr<float4> d_deltas1((float4*) thrust::raw_pointer_cast(deltas.data()));
        thrust::device_ptr<float4> d_deltas2 = d_deltas1 + numConstraints;

        thrust::for_each(
                    thrust::make_zip_iterator(thrust::make_tuple(d_indices, d_dists, d_sortedI1, d_sortedI2, d_deltas1, d_deltas2)),
                    thrust::make_zip_iterator(thrust::make_tuple(d_indices+numConstraints, d_dists+numConstraints,
                                                                 d_sortedI1+numConstraints, d_sortedI2+numConstraints,
                                                                 d_deltas1+numConstraints, d_deltas2+numConstraints)),
            delta_computing_functor((float4 *)particles));

        thrust::sort_by_key(sortedI.begin(), sortedI.end(), d_deltas1);

        thrust::pair<thrust::device_ptr<uint>,thrust::device_ptr<float4> > new_end;
        new_end = thrust::reduce_by_key(d_sortedI1, d_sortedI1+numConstraints*2, d_deltas1, d_sortedI1, d_deltas1);

        uint size = new_end.first - d_sortedI1;
        thrust::device_ptr<uint> d_occ(occurences.data());

        thrust::for_each(
            thrust::make_zip_iterator(thrust::make_tuple(d_pos4, d_deltas1, d_occ)),
            thrust::make_zip_iterator(thrust::make_tuple(d_pos4+size, d_deltas1 + size, d_occ+size)),
            distance_constraint_functor());
    }

}
