#ifndef INCLUDES_H
#define INCLUDES_H

#define GLM_FORCE_RADIANS

// Standard includes
#include <stdlib.h>
#include <iostream>

// GL includes
#define GL_GLEXT_PROTOTYPES
#include <qgl.h>
#include <GL/glu.h>

// GLM includes
#include <glm/vec2.hpp>
#include <glm.hpp>
#include <glm/gtx/rotate_vector.hpp>

// Qt data includes
#include <QList>
#include <QHash>

// Generally helpful functions
inline float frand() { return (double)rand() / (double)RAND_MAX; }
inline float urand(double a, double b) { return a + (b - a) * frand(); }

using namespace std;

inline void printVec(const glm::dvec2 v) {
    cout << "(" << v.x << ", " << v.y << ")" << endl;
}

#define EPSILON .0001

#define D2R(d) (d * M_PI / 180)
#define R2D(r) (r * 180 / M_PI)

#endif // INCLUDES_H
