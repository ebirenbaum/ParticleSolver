#ifndef UTIL_CUH
#define UTIL_CUH

typedef unsigned int uint;

extern "C"
{
    void cudaInit();

    void allocateArray(void **devPtr, int size);
    void freeArray(void *devPtr);

    void copyArrayToDevice(void *device, const void *host, int offset, int size);

    void registerGLBufferObject(unsigned int vbo, struct cudaGraphicsResource **cuda_vbo_resource);
    void unregisterGLBufferObject(struct cudaGraphicsResource *cuda_vbo_resource);

    void *mapGLBufferObject(struct cudaGraphicsResource **cuda_vbo_resource);
    void unmapGLBufferObject(struct cudaGraphicsResource *cuda_vbo_resource);

    void copyArrayFromDevice(void *host, const void *device, int size);

    uint iDivUp(uint a, uint b);
    void computeGridSize(uint n, uint blockSize, uint &numBlocks, uint &numThreads);
}

#endif // UTIL_CUH
