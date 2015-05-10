#ifndef ORBITINGCAMERA_H
#define ORBITINGCAMERA_H

#define GLM_FORCE_RADIANS
#include <glm.hpp>

/**
 * @class OrbitingCamera
 *
 * Camera that move in response to mouse interaction.
 *
 * You shouldn't need to work on this class. It's there for your convenience, really,
 * and the way this camera is implemented is NOT the way you should be implementing
 * your camera in the Camtrans lab. We hide the real implementation by using OpenGL to
 * perform the camera calculations.
 *
 */
class OrbitingCamera
{
public:
    OrbitingCamera();
    virtual ~OrbitingCamera();

    virtual void setAspectRatio(float aspectRatio);

    virtual glm::mat4 getProjectionMatrix() const;
    virtual glm::mat4 getViewMatrix() const;
    virtual glm::mat4 getScaleMatrix() const;
    virtual glm::mat4 getProjectionViewMatrix() const;

    virtual float getFovY() const;

    virtual void mouseDown(int x, int y);
    virtual void mouseDragged(int x, int y);
    virtual void mouseScrolled(int delta);

    void updateMatrices();

private:
    void updateProjectionMatrix();
    void updateViewMatrix();

    glm::mat4 m_viewMatrix;
    glm::mat4 m_projectionMatrix;
    glm::mat4 m_scaleMatrix;
    glm::mat4 m_projView;
    float m_aspectRatio, m_angleX, m_angleY, m_zoomZ;
    int m_oldX, m_oldY;
    float m_eyeHeight;
    float m_eyePan;
};

#endif // ORBITINGCAMERA_H
