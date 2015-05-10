#ifndef CAMERA_H
#define CAMERA_H

#include <glm.hpp>

class Camera
{
public:
    Camera();
    virtual ~Camera();

    glm::mat4 getProjectionViewMatrix();
    glm::mat4 getProjectionMatrix();
    glm::mat4 getViewMatrix();
    glm::mat4 getScaleMatrix();
    glm::mat4 getFrustumMatrix();
    glm::vec4 getLook();
    glm::vec4 getUp();
    glm::vec4 getEye();
    float getAspectRatio();

    void setAspectRatio(float a);
    void orientLook(glm::vec4 &eye, glm::vec4 &look, glm::vec4 &up);

    void moveHorizontal(glm::vec2 dir);
    void moveAlongU(float mag);
    void moveAlongUp(float mag);
    void moveAlongLook(float mag);

    void pitch(float degrees);
    void yaw(float degrees);
    void roll(float degrees);

protected:
    void setCameraSpace();
    void setViewMatrix();
    void setProjectionMatrix();
    void setFrustumMatrix();

    glm::vec4 m_u, m_v, m_w;
    glm::vec4 m_eye, m_look, m_up;

    glm::mat4 m_view, m_proj, m_frustum;
    glm::mat4 m_scale;

    // View variables
    float m_near, m_far, m_heightDegrees, m_aspectRatio;

    float m_thirdDist;

};

#endif // CAMERA_H
