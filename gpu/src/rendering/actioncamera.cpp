#include "actioncamera.h"

ActionCamera::ActionCamera()
{
    m_pos = glm::vec4();
    m_offset = 0.f;
    m_offsetHeight = 0.f;
}

ActionCamera::~ActionCamera()
{
}


float ActionCamera::getOffset()
{
    return m_offset;
}

void ActionCamera::setLook(glm::vec4 look)
{
    m_look = look;
    glm::vec4 up = glm::vec4(0, 1, 0, 0);
    orientLook(m_eye, m_look, up);
}

void ActionCamera::setCenter(glm::vec3 pos)
{
    m_pos = glm::vec4(pos, 1.f);
    glm::vec4 eye = m_pos - m_look * m_offset + m_up * m_offsetHeight * m_offset;
    orientLook(eye, m_look, m_up);
}


void ActionCamera::setOffsetHeight(float height)
{
    m_offsetHeight = height;
    glm::vec4 eye = m_pos - m_look * m_offset + m_up * m_offsetHeight * m_offset;
    orientLook(eye, m_look, m_up);
}

void ActionCamera::setOffset(float offset)
{
    m_offset = offset;
    glm::vec4 eye = m_pos - m_look * m_offset + m_up * m_offsetHeight * m_offset;
    orientLook(eye, m_look, m_up);
}


void ActionCamera::moveRelativeToLook(glm::vec3 dir)
{
    moveAlongU(dir.x);
    moveAlongUp(dir.y);
    moveAlongLook(dir.z);
}
