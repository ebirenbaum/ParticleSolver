#include "camera.h"

#define GLM_FORCE_RADIANS
#include <gtx/rotate_vector.hpp>
#include <gtx/vector_angle.hpp>

Camera::Camera()
{
    glm::vec4 eye = glm::vec4(0, 0, 5, 0);
    glm::vec4 look = -eye;
    glm::vec4 up = glm::vec4(0, 1, 0, 0);
    orientLook(eye, look, up);

    // Projection Defaults
    m_heightDegrees = 90.f;
    m_aspectRatio = 1.0f;
    m_near = 0.1f;
    m_far = 1000.0f;
    setProjectionMatrix();

    m_thirdDist = 0.f;
    setFrustumMatrix();
}

Camera::~Camera()
{
}

glm::mat4 Camera::getProjectionViewMatrix()
{
    return m_proj * m_view;
}

glm::mat4 Camera::getProjectionMatrix()
{
    return m_proj;
}

glm::mat4 Camera::getViewMatrix()
{
    return m_view;
}

glm::mat4 Camera::getScaleMatrix()
{
    return m_scale;
}

glm::mat4 Camera::getFrustumMatrix()
{
    return m_frustum;
}

glm::vec4 Camera::getLook()
{
    return m_look;
}

glm::vec4 Camera::getUp()
{
    return m_up;
}

glm::vec4 Camera::getEye()
{
    return m_eye;
}

float Camera::getAspectRatio()
{
    return m_aspectRatio;
}

void Camera::setAspectRatio(float a)
{
    m_aspectRatio = a;
    setProjectionMatrix();
    setFrustumMatrix();
}

void Camera::orientLook(glm::vec4 &eye, glm::vec4 &look, glm::vec4 &up)
{
    // Camera vecs
    m_eye = eye;
    m_look = look;
    m_up = up;

    setCameraSpace();
    setViewMatrix();
    setFrustumMatrix();
}

void Camera::moveHorizontal(glm::vec2 dir)
{
    m_eye += glm::normalize(glm::vec4(m_look.x, 0.f, m_look.z, 0.f)) * dir.x;
    m_eye += glm::normalize(glm::vec4(-m_look.z, 0.f, m_look.x, 0.f)) * dir.y;
    setViewMatrix();
    setFrustumMatrix();
}

void Camera::moveAlongU(float mag)
{
    m_eye += m_u * mag;
    setViewMatrix();
    setFrustumMatrix();
}

void Camera::moveAlongUp(float mag)
{
    m_eye += m_up * mag;
    setViewMatrix();
    setFrustumMatrix();
}

void Camera::moveAlongLook(float mag)
{
    m_eye += m_look * mag;
    setViewMatrix();
    setFrustumMatrix();
}

void Camera::pitch(float degrees)
{
    glm::vec4 oldLook = m_look;
    m_look = glm::rotate(m_look, glm::radians(degrees), glm::cross(glm::vec3(m_up), glm::vec3(m_look)));
    setCameraSpace();
    if (m_up.y < 0)
    {
        m_look = oldLook;
        setCameraSpace();
    }
    setViewMatrix();
    setFrustumMatrix();
}

void Camera::yaw(float degrees)
{
    float radians = glm::radians(degrees);
    glm::vec3 vec = glm::vec3(0.f, -1.f, 0.f);
    m_look = glm::rotate(m_look, radians, vec);
    m_up = glm::rotate(m_up, radians, vec);
    setCameraSpace();
    setViewMatrix();
    setFrustumMatrix();
}

void Camera::roll(float degrees)
{
    m_up = glm::rotate(m_up, glm::radians(degrees), glm::vec3(m_look));
    setCameraSpace();
    setViewMatrix();
    setFrustumMatrix();
}

void Camera::setCameraSpace()
{
    // Camera coordinate space
    m_look = glm::normalize(m_look);
    m_w = -m_look;
    m_v = glm::normalize(m_up - (glm::dot(m_up, m_w) * m_w));
    m_u = glm::vec4(
                glm::normalize(
                        glm::cross(glm::vec3(m_v.x, m_v.y, m_v.z),
                                   glm::vec3(m_w.x, m_w.y, m_w.z))),
                    0);
    m_up = m_v;
}

void Camera::setViewMatrix()
{
    // View Matrices
    glm::mat4 trans = glm::mat4();
    trans[3][0] = -m_eye[0];
    trans[3][1] = -m_eye[1];
    trans[3][2] = -m_eye[2];

    glm::mat4 rot = glm::mat4(m_u.x, m_u.y, m_u.z, 0.0,
                              m_v.x, m_v.y, m_v.z, 0.0,
                              m_w.x, m_w.y, m_w.z, 0.0,
                               0.0,   0.0,   0.0,  1.0);
    rot = glm::transpose(rot);
    m_view = rot * trans;
}

void Camera::setProjectionMatrix()
{
    // Projection Matrices
    float h = m_far * glm::tan(glm::radians(m_heightDegrees / 2.0));
    float w = m_aspectRatio * h;

    m_scale = glm::mat4(1.0 / (w / 2.0),       0.0,          0.0,     0.0,
                                     0.0,        1.0 / (h / 2.0),    0.0,     0.0,
                                     0.0,              0.0,       1.0 / m_far,  0.0,
                                     0.0,              0.0,          0.0,     1.0);
    m_scale = glm::transpose(m_scale);

    float c = -m_near / m_far;
    glm::mat4 perspective = glm::mat4(1.0,   0.0,      0.0,         0.0,
                                      0.0,   1.0,      0.0,         0.0,
                                      0.0,   0.0, -1.0/(1.0+c),  c/(1.0+c),
                                      0.0,   0.0,     -1.0,         0.0);
    perspective = glm::transpose(perspective);
    m_proj = perspective * m_scale;
}

void Camera::setFrustumMatrix()
{
    m_frustum = glm::transpose(m_proj * m_view);
}



