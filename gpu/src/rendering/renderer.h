#ifndef RENDERER_H
#define RENDERER_H

#include <glm/glm.hpp>
#include <vector_types.h>
#include <vector>

typedef unsigned int GLuint;
typedef unsigned int uint;
class ActionCamera;
class QMouseEvent;
class QKeyEvent;

class Renderer
{
public:
    Renderer();
    ~Renderer();

    void createVAO(GLuint vbo, float radius);

    void setVBO(GLuint vbo, uint numParticles);
    void render(std::vector<int2> colorIndices, std::vector<float4> colors);

    float4 raycast2XYPlane(float x, float y);
    float3 getDir(float x, float y);
    float3 getEye();

//    void mousePressed(QMouseEvent *e);
    void mouseMoved(QMouseEvent *, float deltaX, float deltaY);
//    void mouseScrolled(QWheelEvent *e);

    void keyPressed(QKeyEvent* e);
    void keyReleased(QKeyEvent* e);

    void update(float secs);

    void resize(int w, int h);

private:
    void _initGL();
    GLuint _compileProgram(const char *vertex_file_path, const char *fragment_file_path);

    GLuint m_program;
    GLuint m_vbo;
    GLuint m_vao;

    ActionCamera *m_camera;
    glm::ivec2 m_screenSize;

    float m_particleRadius;
    int m_wsadeq;

};

#endif // RENDERER_H
