#ifndef ACTIONCAMERA_H
#define ACTIONCAMERA_H

#include "camera.h"

class ActionCamera : public Camera
{
public:
    ActionCamera();
    ~ActionCamera();

    float getOffset();

    void setCenter(glm::vec3 pos);
    void setOffset(float offset);
    void setLook(glm::vec4 look);
    void setOffsetHeight(float height);

    void moveRelativeToLook(glm::vec3 dir);

private:
    glm::vec4 m_pos;//, m_offsetVec;
    float m_offset;
    float m_offsetHeight;
};

#endif // ACTIONCAMERA_H
