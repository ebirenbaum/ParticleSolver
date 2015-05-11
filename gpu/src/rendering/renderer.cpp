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

//#include"debugprinting.h"

Renderer::Renderer(int3 minBounds, int3 maxBounds)
    : m_program(0),
      m_vbo(0),
      m_vao(0),
      m_vboGrid(0),
      m_vaoGrid(0),
      m_numGridVerts(0),
      m_camera(NULL),
      m_particleRadius(0),
      m_wsadeq(0)
{
    _initGL();
    m_camera = new ActionCamera();
    m_camera->setOffset(0.f);
    m_camera->setOffsetHeight(0.f);
    m_camera->setCenter(glm::vec3(0, 10, 25));

    _buildGrid(minBounds, maxBounds);
}


Renderer::~Renderer()
{
    if (m_vbo)
        glDeleteBuffers(1, &m_vbo);
    if (m_vao)
        glDeleteVertexArrays(1, &m_vao);

    if (m_vboGrid)
        glDeleteBuffers(1, &m_vboGrid);
    if (m_vaoGrid)
        glDeleteVertexArrays(1, &m_vaoGrid);

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
    glUniform1f(glGetUniformLocation(m_program, "particleRadius"), -1.f);
    glUniform1i(glGetUniformLocation(m_program, "screenHeight"), m_screenSize.y);

    GLuint colorLoc = glGetUniformLocation(m_program, "color");

    // Draw floor
    glBindVertexArray(m_vaoGrid);
    glUniform4f(colorLoc, .1f, .1f, .1f, 1.f);
    glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);

    // Draw grid and walls
    glUniform4f(colorLoc, .5f, .5f, .5f, 1.f);
    glDrawArrays(GL_LINES, 4, m_numGridVerts);


    // Draw points
    glUniform1f(glGetUniformLocation(m_program, "particleRadius"), m_particleRadius);
    glEnable(GL_PROGRAM_POINT_SIZE);

    glBindVertexArray(m_vao);

    int size = colors.size();
    int2 index;
    float4 color;
    for (int i = 0; i < size; i++)
    {
        index = colorIndices.at(i);
        color = colors.at(i);
        glUniform4f(colorLoc, color.x, color.y, color.z, color.w);
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

void Renderer::_buildGrid(int3 minBounds, int3 maxBounds)
{
    int spacing = 4;
    int maxDistX = maxBounds.x - minBounds.x;
    int maxDistZ = maxBounds.z - minBounds.z;
    int numX = maxDistX / spacing;
    int numZ = maxDistZ / spacing;

    if (maxDistX % spacing != 0)
        numX++;
    if (maxDistZ % spacing != 0)
        numZ++;

    numX = 2 * (numX + 1);
    numZ = 2 * (numZ + 1);

    m_numGridVerts =  numX + numZ + 16; // 16 for wall lines
    int size = (m_numGridVerts + 4) * 3;
    float *data = new float[size];

    data[0] = maxBounds.x; data[1] = minBounds.y; data[2] = minBounds.z;
    data[3] = minBounds.x; data[4] = minBounds.y; data[5] = minBounds.z;
    data[6] = maxBounds.x; data[7] = minBounds.y; data[8] = maxBounds.z;
    data[9] = minBounds.x; data[10] = minBounds.y; data[11] = maxBounds.z;

    int index = 12;
    int i;
    bool flop = true;

    // X
    for (i = minBounds.x; i <= maxBounds.x; i+=spacing)
    {
        if (flop)
        {
            data[index++] = i; data[index++] = minBounds.y, data[index++] = minBounds.z;
            data[index++] = i; data[index++] = minBounds.y, data[index++] = maxBounds.z;
        }
        else
        {
            data[index++] = i; data[index++] = minBounds.y, data[index++] = maxBounds.z;
            data[index++] = i; data[index++] = minBounds.y, data[index++] = minBounds.z;
        }
        flop = !flop;
    }
    if (i-spacing != maxBounds.x)
    {
        if (flop)
        {
            data[index++] = maxBounds.x; data[index++] = minBounds.y, data[index++] = minBounds.z;
            data[index++] = maxBounds.x; data[index++] = minBounds.y, data[index++] = maxBounds.z;
        }
        else
        {
            data[index++] = maxBounds.x; data[index++] = minBounds.y, data[index++] = maxBounds.z;
            data[index++] = maxBounds.x; data[index++] = minBounds.y, data[index++] = minBounds.z;
        }
        flop = !flop;
    }


    // Z
    flop = true;
    for (i = minBounds.z; i <= maxBounds.z; i+=spacing)
    {
        if (flop)
        {
            data[index++] = minBounds.x; data[index++] = minBounds.y, data[index++] = i;
            data[index++] = maxBounds.x; data[index++] = minBounds.y, data[index++] = i;
        }
        else
        {
            data[index++] = maxBounds.x; data[index++] = minBounds.y, data[index++] = i;
            data[index++] = minBounds.x; data[index++] = minBounds.y, data[index++] = i;
        }
        flop = !flop;
    }
    if (i-spacing != maxBounds.z)
    {
        if (flop)
        {
            data[index++] = minBounds.x; data[index++] = minBounds.y, data[index++] = maxBounds.z;
            data[index++] = maxBounds.x; data[index++] = minBounds.y, data[index++] = maxBounds.z;
        }
        else
        {
            data[index++] = maxBounds.x; data[index++] = minBounds.y, data[index++] = maxBounds.z;
            data[index++] = minBounds.x; data[index++] = minBounds.y, data[index++] = maxBounds.z;
        }
        flop = !flop;
    }


    // Walls
    data[index++] = minBounds.x; data[index++] = minBounds.y; data[index++] = minBounds.z;
    data[index++] = minBounds.x; data[index++] = maxBounds.y; data[index++] = minBounds.z;
    data[index++] = minBounds.x; data[index++] = maxBounds.y; data[index++] = minBounds.z;

    data[index++] = maxBounds.x; data[index++] = maxBounds.y; data[index++] = minBounds.z;
    data[index++] = maxBounds.x; data[index++] = minBounds.y; data[index++] = minBounds.z;
    data[index++] = maxBounds.x; data[index++] = maxBounds.y; data[index++] = minBounds.z;
    data[index++] = maxBounds.x; data[index++] = maxBounds.y; data[index++] = minBounds.z;

    data[index++] = maxBounds.x; data[index++] = maxBounds.y; data[index++] = maxBounds.z;
    data[index++] = maxBounds.x; data[index++] = minBounds.y; data[index++] = maxBounds.z;
    data[index++] = maxBounds.x; data[index++] = maxBounds.y; data[index++] = maxBounds.z;
    data[index++] = maxBounds.x; data[index++] = maxBounds.y; data[index++] = maxBounds.z;

    data[index++] = minBounds.x; data[index++] = maxBounds.y; data[index++] = maxBounds.z;
    data[index++] = minBounds.x; data[index++] = minBounds.y; data[index++] = maxBounds.z;
    data[index++] = minBounds.x; data[index++] = maxBounds.y; data[index++] = maxBounds.z;
    data[index++] = minBounds.x; data[index++] = maxBounds.y; data[index++] = maxBounds.z;

    data[index++] = minBounds.x; data[index++] = maxBounds.y; data[index++] = minBounds.z;

//    cout << index << ", " << size << endl;

//    m_numGridVerts = 0;
//    size = 12;

    _setGridBuffer(data, size * sizeof(float));

    delete [] data;
}

void Renderer::_setGridBuffer(float *data, int memsize)
{
    if (m_vboGrid)
        glDeleteBuffers(1, &m_vboGrid);
    if (m_vaoGrid)
        glDeleteVertexArrays(1, &m_vaoGrid);

    // set and bind vertex buffer object.
    glGenBuffers(1, &m_vboGrid);
    glBindBuffer(GL_ARRAY_BUFFER, m_vboGrid);

    // Initialize the vertex array object.
    glGenVertexArrays(1, &m_vaoGrid);
    glBindVertexArray(m_vaoGrid);

    glBufferData(GL_ARRAY_BUFFER, memsize, data, GL_STATIC_DRAW);

    GLuint position = glGetAttribLocation(m_program, "position");

    glEnableVertexAttribArray(position);
    glVertexAttribPointer(
        position,
        3,                   // Num coordinates per position
        GL_FLOAT,            // Type
        GL_FALSE,            // Normalized
        sizeof(GLfloat) * 3, // Stride
        (void*) 0            // Array buffer offset
    );

    // Unbind buffers.
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
}

