#version 410 core

uniform vec4 color;
uniform mat4 view;

out vec4 fragColor;

void main()
{
    vec3 camPos = -vec3(view[3]);
    vec3 lightDir = camPos + vec3(-camPos.z, camPos.y, camPos.x);
    lightDir = normalize(lightDir);
    // calculate normal from texture coordinates
    vec3 N;
    N.xy = gl_PointCoord * vec2(2.0, -2.0) + vec2(-1.0, 1.0);
    float mag = dot(N.xy, N.xy);

    if (mag > 1.0) discard;   // kill pixels outside circle

    N.z = sqrt(1.0-mag);

    // calculate lighting
    vec3 diffuse = vec3(max(0.0, dot(lightDir, N)));
    vec3 shadingColor = diffuse + vec3(.2); // plus ambient

//    fragColor = vec4(clamp(shadingColor, 0, 1), 1);
//    fragColor = vec4(color);
    fragColor = vec4(color.xyz * shadingColor, color.w);
//    fragColor = vec4(0,0,color.z,.5);
}
