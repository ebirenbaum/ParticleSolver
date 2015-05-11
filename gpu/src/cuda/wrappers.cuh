#ifndef WRAPPERS_CUH
#define WRAPPERS_CUH

typedef unsigned int uint;

extern "C"
{

    /*
     *   INTEGRATION
     */
    void initIntegration();

    void appendIntegrationParticle(float *v, float *ro, uint iterations);

    void freeIntegrationVectors();

    void setParameters(SimParams *hostParams);

    void integrateSystem(float *pos,
                         float deltaTime,
                         uint numParticles);

    void calcHash(uint  *gridParticleHash,
                  uint  *gridParticleIndex,
                  float *pos,
                  int    numParticles);

    void sortParticles(uint *dGridParticleHash, uint *dGridParticleIndex, uint numParticles);

    void reorderDataAndFindCellStart(uint  *cellStart,
                                     uint  *cellEnd,
                                     float *sortedPos,
                                     float *sortedW,
                                     int   *sortedPhase,
                                     uint  *gridParticleHash,
                                     uint  *gridParticleIndex,
                                     float *oldPos,
                                     uint   numParticles,
                                     uint   numCells);

    void collideWorld(float *pos,
                      float *sortedPos,
                      uint numParticles,
                      int3 minBounds,
                      int3 maxBounds);

    void collide(float *particles,
                 float *sortedPos,
                 float *sortedW,
                 int   *sortedPhase,
                 uint  *gridParticleIndex,
                 uint  *cellStart,
                 uint  *cellEnd,
                 uint   numParticles,
                 uint   numCells);

    void sortByType(float *dPos, uint numParticles);


    void calcVelocity(float *dpos, float deltaTime, uint numParticles);


    /*
     * SOLVER
     */
    void initHandles();
    void destroyHandles();

    void appendSolverParticle(uint numParticles);

    void addPointConstraint(uint *index, float *point, uint numConstraints);
    void addDistanceConstraint(uint *index, float *distance, uint numConstraints);

    void freeSolverVectors();

    void solvePointConstraints(float *particles);

    void solveDistanceConstraints(float *particles);

    ////////////////////////////////// FLUIDS ////////////////////////
    void solveFluids(float *sortedPos,
                     int   *sortedPhase,
                     uint  *gridParticleIndex,
                     uint  *cellStart,
                     uint  *cellEnd,
                     float *particles,
                     uint   numParticles,
                     uint   numCells,
                     float4 mousePos);
}

#endif // WRAPPERS_CUH
