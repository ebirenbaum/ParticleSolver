/**
 * @file   OrbitingCamera.cpp
 *
 * (See the header file.) You don't need to be poking around in this file unless you're interested
 * in how an orbiting camera works.
 *
 * The way we have implemented this class is NOT how you should be implementing your Camtrans. This
 * camera is a DIFFERENT TYPE of camera which we're providing so you can easily view your Shapes
 * and to make sure your scene graph is working if your camera isn't.
 *
 * In the Camtrans lab, you'll be implementing your own perspective camera from scratch! This one
 * uses OpenGL.
 */

#include "orbitingcamera.h"
#include "kernel.cuh"

#define GLM_FORCE_RADIANS
//#include <gtc/matrix_transform.hpp>
#include <gtx/transform.hpp>

#include <float.h>
#include <math.h>

#define NEAR 0.01f
#define FAR 1000.f
#define FOV 90.f

OrbitingCamera::OrbitingCamera()
{
    m_aspectRatio = 1;
    m_angleX = m_angleY = 0;
#ifdef TWOD
    m_zoomZ = -20;
    m_eyeHeight = 10.f;
#else
    m_eyeHeight = 5.f;
    m_zoomZ = -20;
    m_angleY = -45.f;
#endif
    m_eyePan = 0.f;

    updateMatrices();
}

OrbitingCamera::~OrbitingCamera()
{
}

void OrbitingCamera::setAspectRatio(float aspectRatio)
{
    m_aspectRatio = aspectRatio;

    updateProjectionMatrix();
}

glm::mat4 OrbitingCamera::getProjectionMatrix() const
{
    return m_projectionMatrix;
}

glm::mat4 OrbitingCamera::getViewMatrix() const
{
    return m_viewMatrix;
}


glm::mat4 OrbitingCamera::getScaleMatrix() const
{
    float far = glm::max(FAR, NEAR + 100.f * FLT_EPSILON);
    float h = far * glm::tan(glm::radians(FOV / 2.f));
    float w = m_aspectRatio * h;

    glm::mat4 scale = glm::mat4(1.0 / w,       0.0,          0.0,     0.0,
                                    0.0,       1.0 / h,        0.0,     0.0,
                                    0.0,         0.0,       1.0 / far,  0.0,
                                    0.0,         0.0,          0.0,     1.0);
    return glm::transpose(scale);
}

glm::mat4 OrbitingCamera::getProjectionViewMatrix() const
{
    return m_projView;
}

float OrbitingCamera::getFovY() const
{
    return /*glm::radians(*/FOV / m_aspectRatio;
}

void OrbitingCamera::mouseDown(int x, int y)
{
    m_oldX = x;
    m_oldY = y;
}

void OrbitingCamera::mouseDragged(int x, int y)
{
    m_angleY += (x - m_oldX) * .5f;
    m_angleX += (y - m_oldY) * .5f;

#ifdef TWOD
    m_eyeHeight += (y - m_oldY) * .05f;
    m_eyePan += (x - m_oldX) * .05f;
#else
    if (m_angleX < -90) m_angleX = -90;
    if (m_angleX > 90) m_angleX = 90;
#endif

    m_oldX = x;
    m_oldY = y;

    updateViewMatrix();
}

void OrbitingCamera::mouseScrolled(int delta)
{
    // Use an exponential factor so the zoom increments are small when we are
    // close to the object and large when we are far away from the object
    m_zoomZ *= powf(0.999f, delta);

    updateViewMatrix();
}

void OrbitingCamera::updateMatrices()
{
    updateProjectionMatrix();
    updateViewMatrix();
//    updateScaleMatrix();
}

void OrbitingCamera::updateProjectionMatrix()
{
    // Make sure glm gets a far value that is greater than the near value.
    // Thanks Windows for making far a keyword!
    float farPlane = glm::max(FAR, NEAR + 100.f * FLT_EPSILON);

    m_projectionMatrix = glm::perspective(
            glm::radians(FOV), m_aspectRatio, NEAR, farPlane);

    m_projView = m_projectionMatrix * m_viewMatrix;
}

void OrbitingCamera::updateViewMatrix()
{

#ifdef TWOD
    m_viewMatrix = glm::translate(glm::mat4(), glm::vec3(m_eyePan, -m_eyeHeight, m_zoomZ));
#else
    m_viewMatrix =
            glm::translate(glm::vec3(0.f, -m_eyeHeight, m_zoomZ)) *
            glm::rotate(glm::radians(m_angleY), glm::vec3(0.f, 1.f, 0.f)) *
            glm::rotate(glm::radians(m_angleX), glm::vec3(1.f, 0.f, 0.f));
#endif

    m_projView = m_projectionMatrix * m_viewMatrix;
}
