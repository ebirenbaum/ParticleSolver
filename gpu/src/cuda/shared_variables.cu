

#include "thrust/device_vector.h"
#include "helper_cuda.h"
#include "cuda_runtime.h"
#include "util.cuh"

thrust::device_vector<float> Xstar;	// guess vectors
thrust::device_vector<float> W;     // vector of inverse masses
thrust::device_vector<int> phase;


// textures
texture<float4, 1, cudaReadModeElementType> oldPosTex;
texture<float4, 1, cudaReadModeElementType> oldVelTex;


extern "C"
{

	void freeSharedVectors()
	{
        Xstar.clear();
        W.clear();
        phase.clear();

        Xstar.shrink_to_fit();
        W.shrink_to_fit();
        phase.shrink_to_fit();
	}

    void appendPhaseAndMass(int *fase, float *w, uint numParticles)
    {
        int sizeW = W.size();

        // resize the vectors
        phase.resize(sizeW + numParticles);
        W.resize(sizeW + numParticles);

        // get raw pointers to the data
        int *dPhase = thrust::raw_pointer_cast(phase.data());
        float *dW = thrust::raw_pointer_cast(W.data());

        // copy the new data to the gpu
        copyArrayToDevice(dPhase + sizeW, fase, 0, numParticles * sizeof(int));
        copyArrayToDevice(dW + sizeW, w, 0, numParticles * sizeof(float));

        // resize but don't neet to fill
        Xstar.resize(4 * W.size());
    }

	void copyToXstar(float *pos, uint numParticles)
	{
        // copy X to X*
        float *dXstar = thrust::raw_pointer_cast(Xstar.data());
        checkCudaErrors(cudaMemcpy((void*)dXstar, (void*)pos, numParticles*4*sizeof(float), cudaMemcpyDeviceToDevice));
    }

    int *getPhaseRawPtr()
    {
        return thrust::raw_pointer_cast(phase.data());
    }

    float *getXstarRawPtr()
    {
        return thrust::raw_pointer_cast(Xstar.data());
    }

    float *getWRawPtr()
    {
        return thrust::raw_pointer_cast(W.data());
    }

    void printXstar()
    {
        printf("Xstar: size: %u\n", (uint)Xstar.size());
        thrust::device_ptr<float> d_Xstar(Xstar.data());
        uint index;
        for (uint i = 0; i < Xstar.size(); i++)
        {
            index = i * 4;
            printf("i: %u: %.2f, %.2f, %.2f\n", i, (float)*(d_Xstar + index + 0), (float)*(d_Xstar + index + 1), (float)*(d_Xstar + index + 2));
        }
        printf("\n");
    }

}
