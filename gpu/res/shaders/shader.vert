#version 410 core

in vec4 position;           // position of point in world space

uniform mat4 pv;            // projection * view matrix
uniform mat4 projection;
uniform mat4 view;
uniform float particleRadius = -1.0; // world space particle size

uniform int screenHeight;

const float PI = 3.1415926535;

void main()
{
    gl_Position = pv * vec4(position.xyz, 1); // storing inverse mass in w position

    if (particleRadius > -.5)
        gl_PointSize = screenHeight * projection[1][1] * particleRadius / gl_Position.w;
}
