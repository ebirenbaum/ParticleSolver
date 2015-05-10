#include "renderer.h"
#include <GL/glew.h>
#include <QFile>
#include <QTextStream>
#include <QImage>
#include <QDir>
#include <QGLWidget>
#include <QMouseEvent>
#include <QWheelEvent>
#include <string>
#include <sstream>
#include <vector>
#include "actioncamera.h"
#include "helper_math.h"
#include "kernel.cuh"

#define GLM_FORCE_RADIANS
#include <gtc/type_ptr.hpp>
#include <gtx/norm.hpp>

#include"debugprinting.h"

Renderer::Renderer()
    : m_program(0),
      m_vbo(0),
      m_vao(0),
      m_camera(NULL),
      m_particleRadius(0),
      m_wsadeq(0)
{
    _initGL();
    m_camera = new ActionCamera();
    m_camera->setOffset(0.f);
    m_camera->setOffsetHeight(0.f);
    m_camera->setCenter(glm::vec3(0, 10, 25));
//#ifdef TWOD
//    // 2D
//    glm::vec4 eye = glm::vec4(0,10,200,0);
//    glm::vec4 look = glm::vec4(0,.2,-1,0);
//#else
//    // 3D
//    glm::vec4 eye = glm::vec4(0,10,50,0);
//    glm::vec4 look = glm::vec4(0,-.15,-1,0);
//#endif
//    // 3D cloth
////    glm::vec4 eye = glm::vec4(30,20,30,0);
////    glm::vec4 look = glm::vec4(-1,-.3,-1,0);
//    glm::vec4 up = glm::vec4(0,1,0,0);
//    m_camera->orientLook(eye, look , up);
}


Renderer::~Renderer()
{
    if (m_vbo)
        glDeleteBuffers(1, &m_vbo);
    if (m_vao)
        glDeleteVertexArrays(1, &m_vao);
    if (m_camera)
        delete m_camera;
}


void Renderer::createVAO(GLuint vbo, float radius)
{
    if (m_vbo && m_vbo != vbo)
        glDeleteBuffers(1, &m_vbo);
    if (m_vao)
        glDeleteVertexArrays(1, &m_vao);

    m_vbo = vbo;
    m_particleRadius = radius;

    // Initialize the vertex array object.
    glGenVertexArrays(1, &m_vao);
    glBindVertexArray(m_vao);

    // bind vertex buffer object.
    glBindBuffer(GL_ARRAY_BUFFER, m_vbo);

    GLuint position = glGetAttribLocation(m_program, "position");

    glEnableVertexAttribArray(position);
    glVertexAttribPointer(
        position,
        4,                   // Num coordinates per position
        GL_FLOAT,            // Type
        GL_FALSE,            // Normalized
        sizeof(GLfloat) * 4, // Stride
        (void*) 0            // Array buffer offset
    );

    // Unbind buffers.
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
}


void Renderer::render(std::vector<int2> colorIndices, std::vector<float4> colors)
{
    glEnable(GL_BLEND); //Enable blending.
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glm::mat4 projView = m_camera->getProjectionViewMatrix();
    glm::mat4 proj = m_camera->getProjectionMatrix();
    glm::mat4 view = m_camera->getViewMatrix();

    glUseProgram(m_program);

    glUniformMatrix4fv(glGetUniformLocation(m_program, "pv"), 1, GL_FALSE, glm::value_ptr(projView));
    glUniformMatrix4fv(glGetUniformLocation(m_program, "projection"), 1, GL_FALSE, glm::value_ptr(proj));
    glUniformMatrix4fv(glGetUniformLocation(m_program, "view"), 1, GL_FALSE, glm::value_ptr(view));
    glUniform1f(glGetUniformLocation(m_program, "particleRadius"), m_particleRadius);
    glUniform1i(glGetUniformLocation(m_program, "screenHeight"), m_screenSize.y);


    glEnable(GL_PROGRAM_POINT_SIZE);

    glBindVertexArray(m_vao);

    int size = colors.size();
    int2 index;
    float4 color;
    for (int i = 0; i < size; i++)
    {
        index = colorIndices.at(i);
        color = colors.at(i);
        glUniform4f(glGetUniformLocation(m_program, "color"), color.x, color.y, color.z, color.w);
        glDrawArrays(GL_POINTS, index.x, index.y - index.x);
    }
    glBindVertexArray(0);

    glDisable(GL_PROGRAM_POINT_SIZE);

    glDisable(GL_BLEND); //Enable blending.
    glUseProgram(0);
}


GLuint Renderer::_compileProgram(const char *vertex_file_path, const char *fragment_file_path)
{
    // NOTE: MUST INIT GLEW BEFORE USING THIS CODE

    // Create the shaders
    GLuint VertexShaderID = glCreateShader(GL_VERTEX_SHADER);
    GLuint FragmentShaderID = glCreateShader(GL_FRAGMENT_SHADER);

    // Read the Vertex Shader code from the file
    std::string VertexShaderCode;
    QString vertFilePath = QString(vertex_file_path);
    QFile vertFile(vertFilePath);
    if (vertFile.open(QIODevice::ReadOnly | QIODevice::Text)){
        QTextStream vertStream(&vertFile);
        VertexShaderCode = vertStream.readAll().toStdString();
    }


    // Read fragment shader code from file
    std::string FragmentShaderCode;
    QString fragFilePath = QString(fragment_file_path);
    QFile fragFile(fragFilePath);
    if (fragFile.open(QIODevice::ReadOnly | QIODevice::Text)){
        QTextStream fragStream(&fragFile);
        FragmentShaderCode = fragStream.readAll().toStdString();
    }

    GLint Result = GL_FALSE;
    int InfoLogLength;

    // Compile Vertex Shader
    char const * VertexSourcePointer = VertexShaderCode.c_str();
    glShaderSource(VertexShaderID, 1, &VertexSourcePointer , NULL);
    glCompileShader(VertexShaderID);

    // Check Vertex Shader
    glGetShaderiv(VertexShaderID, GL_COMPILE_STATUS, &Result);
    glGetShaderiv(VertexShaderID, GL_INFO_LOG_LENGTH, &InfoLogLength);
    std::vector<char> VertexShaderErrorMessage(InfoLogLength);
    glGetShaderInfoLog(VertexShaderID, InfoLogLength, NULL, &VertexShaderErrorMessage[0]);
    if (!Result)
        fprintf(stderr, "Error compiling shader: %s\n%s\n",
                vertex_file_path, &VertexShaderErrorMessage[0]);

    // Compile Fragment Shader
    char const * FragmentSourcePointer = FragmentShaderCode.c_str();
    glShaderSource(FragmentShaderID, 1, &FragmentSourcePointer , NULL);
    glCompileShader(FragmentShaderID);

    // Check Fragment Shader
    glGetShaderiv(FragmentShaderID, GL_COMPILE_STATUS, &Result);
    glGetShaderiv(FragmentShaderID, GL_INFO_LOG_LENGTH, &InfoLogLength);
    std::vector<char> FragmentShaderErrorMessage(InfoLogLength);
    glGetShaderInfoLog(FragmentShaderID, InfoLogLength, NULL, &FragmentShaderErrorMessage[0]);
    if (!Result)
        fprintf(stderr, "Error compiling shader: %s\n%s\n",
                fragment_file_path, &FragmentShaderErrorMessage[0]);

    // Link the program
    GLuint programId = glCreateProgram();
    glAttachShader(programId, VertexShaderID);
    glAttachShader(programId, FragmentShaderID);
    glLinkProgram(programId);

    // Check the program
    glGetProgramiv(programId, GL_LINK_STATUS, &Result);
    glGetProgramiv(programId, GL_INFO_LOG_LENGTH, &InfoLogLength);
    std::vector<char> ProgramErrorMessage( std::max(InfoLogLength, int(1)) );
    glGetProgramInfoLog(programId, InfoLogLength, NULL, &ProgramErrorMessage[0]);
    if (!Result)
        fprintf(stderr, "Error linking shader: %s\n", &ProgramErrorMessage[0]);

    glDeleteShader(VertexShaderID);
    glDeleteShader(FragmentShaderID);

    return programId;
}


float4 Renderer::raycast2XYPlane(float x, float y)
{
    glm::mat4 ftw = glm::inverse(m_camera->getScaleMatrix() * m_camera->getViewMatrix());
    glm::vec4 farWorld = ftw * glm::vec4(x, y, -1, 1);

    glm::vec4 p = glm::inverse(m_camera->getViewMatrix()) * glm::vec4(0,0,0,1);
    glm::vec4 d = glm::normalize(farWorld - p);

#ifdef TWOD
    float t = -ZPOS-p.z / d.z;
#else
    float t = -p.z / d.z;
#endif
    glm::vec4 point = p + d * t;
    return make_float4(point.x, point.y, point.z, point.w);
}

float3 Renderer::getDir(float x, float y)
{
    glm::mat4 ftw = glm::inverse(m_camera->getScaleMatrix() * m_camera->getViewMatrix());
    glm::vec4 farWorld = ftw * glm::vec4(x, y, -1, 1);

    glm::vec4 p = glm::inverse(m_camera->getViewMatrix()) * glm::vec4(0,0,0,1);
    glm::vec4 d = glm::normalize(farWorld - p);

    return make_float3(d.x, d.y, d.z);
}

float3 Renderer::getEye()
{
    glm::vec4 p = glm::inverse(m_camera->getViewMatrix()) * glm::vec4(0,0,0,1);

    return make_float3(p.x, p.y, p.z);
}

void Renderer::_initGL()
{
    m_program = _compileProgram(":/shaders/default.vert", ":/shaders/default.frag");
}

void Renderer::mouseMoved(QMouseEvent *, float deltaX, float deltaY)
{
    m_camera->yaw(deltaX * 30.f);
    m_camera->pitch(deltaY * 30.f);
}

void Renderer::keyPressed(QKeyEvent *e)
{
    switch (e->key())
    {
    case Qt::Key_W:
        m_wsadeq |= 0b100000;
        break;
    case Qt::Key_S:
        m_wsadeq |= 0b010000;
        break;
    case Qt::Key_A:
        m_wsadeq |= 0b001000;
        break;
    case Qt::Key_D:
        m_wsadeq |= 0b000100;
        break;
    case Qt::Key_E:
        m_wsadeq |= 0b000010;
        break;
    case Qt::Key_Q:
        m_wsadeq |= 0b000001;
        break;
    default:
        break;
    }
}


void Renderer::keyReleased(QKeyEvent *e)
{
    switch (e->key())
    {
    case Qt::Key_W:
        m_wsadeq &= 0b011111;
        break;
    case Qt::Key_S:
        m_wsadeq &= 0b101111;
        break;
    case Qt::Key_A:
        m_wsadeq &= 0b110111;
        break;
    case Qt::Key_D:
        m_wsadeq &= 0b111011;
        break;
    case Qt::Key_E:
        m_wsadeq &= 0b111101;
        break;
    case Qt::Key_Q:
        m_wsadeq &= 0b111110;
        break;
    default:
        break;
    }
}

void Renderer::update(float secs)
{
    float forceAmt = 8.f;
    glm::vec3 force = glm::vec3();
    if (m_wsadeq & 0b100000)
        force.z += 1;
    if (m_wsadeq & 0b010000)
        force.z -= 1;
    if (m_wsadeq & 0b001000)
        force.x -= 1;
    if (m_wsadeq & 0b000100)
        force.x += 1;
    if (m_wsadeq & 0b000010)
        force.y += 1;
    if (m_wsadeq & 0b000001)
        force.y -= 1;

    glm::vec4 look = m_camera->getLook();

    glm::vec3 thrust = glm::normalize(glm::vec3(look.x, 0.f, look.z)) * force.z;
    thrust += glm::normalize(glm::vec3(-look.z, 0.f, look.x)) * force.x;
    if (glm::length2(thrust) > 0.00001)
        thrust = glm::normalize(thrust) * forceAmt;
    thrust.y = force.y * forceAmt;

    glm::vec3 eye = glm::vec3(m_camera->getEye());
    eye += thrust * secs;
    m_camera->setCenter(eye);
}


void Renderer::resize(int w, int h)
{
    m_camera->setAspectRatio(w * 1.f / h);
    m_screenSize = glm::ivec2(w, h);
}

