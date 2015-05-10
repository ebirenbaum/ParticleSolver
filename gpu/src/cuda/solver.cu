
//#define PRINT


#include <cuda_runtime.h>
#include <cuda_gl_interop.h>


#include <cublas_v2.h>
#include <cusparse.h>

#include "helper_cuda.h"
#include "solver_kernel.cuh"
#include "util.cuh"
#include "shared_variables.cuh"

#include "thrust/device_ptr.h"
#include "thrust/device_vector.h"
#include "thrust/for_each.h"
#include "thrust/iterator/zip_iterator.h"
#include "thrust/sort.h"
#include "thrust/reduce.h"
#include "thrust/transform.h"

#include <stdio.h>


cublasHandle_t cublasHandle;
cusparseHandle_t cusparseHandle;
cusparseMatDescr_t matDescr;

thrust::device_vector<uint> distsI;
thrust::device_vector<float> dists;

thrust::device_vector<uint> pointsI;
thrust::device_vector<float> points;

thrust::device_vector<uint> sortedI;
thrust::device_vector<float> deltas;

thrust::device_vector<uint> occurences;     // number of constraints affecting a particle

extern "C"
{

    void initHandles()
    {

        checkCudaErrors(cublasCreate(&cublasHandle));
        checkCudaErrors(cusparseCreate(&cusparseHandle));

        checkCudaErrors(cusparseCreateMatDescr(&matDescr));
        cusparseSetMatType(matDescr, CUSPARSE_MATRIX_TYPE_GENERAL);
        cusparseSetMatIndexBase(matDescr,CUSPARSE_INDEX_BASE_ZERO);

    }

    void destroyHandles()
    {
        checkCudaErrors(cublasDestroy(cublasHandle));
        checkCudaErrors(cusparseDestroy(cusparseHandle));
    }

    void appendSolverParticle()
    {
        occurences.push_back(0);
    }

    void addPointConstraint(uint index, float3 point)
    {
        pointsI.push_back(index);
        points.push_back(point.x);
        points.push_back(point.y);
        points.push_back(point.z);
        occurences[index] += 1;
    }

    void addDistanceConstraint(uint2 index, float distance)
    {
        distsI.push_back(index.x);
        distsI.push_back(index.y);
        dists.push_back(distance);
        occurences[index.x] += 1;
        occurences[index.y] += 1;

        sortedI.push_back(0);
        sortedI.push_back(0);

        deltas.push_back(0.f);
        deltas.push_back(0.f);
        deltas.push_back(0.f);
        deltas.push_back(0.f);
        deltas.push_back(0.f);
        deltas.push_back(0.f);
        deltas.push_back(0.f);
        deltas.push_back(0.f);
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
