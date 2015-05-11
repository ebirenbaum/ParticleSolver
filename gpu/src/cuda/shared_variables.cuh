#ifndef DEVICE_VARIABLES_H
#define DEVICE_VARIABLES_H

#define NO_COLLIDE -1
#define FLUID 0
#define GAS 1
#define CLOTH 2
#define SOLID 3
#define RIGID 4

#include "thrust/device_vector.h"

extern "C"
{
	void freeSharedVectors();

    void appendPhaseAndMass(int *fase, float *w, uint numParticles);

	void copyToXstar(float *pos, uint numParticles);
	
	int *getPhaseRawPtr();

    float *getXstarRawPtr();

    float *getWRawPtr();

    void printXstar();
    
    // void bindOldPos();
       
    // void bindOldVel();

    // void unbindOldPos();
       
    // void unbindOldVel();
}


#endif // DEVICE_VARIABLES_H
