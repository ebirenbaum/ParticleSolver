#ifndef PARTICLEAPP_H
#define PARTICLEAPP_H

class QMouseEvent;
class QWheelEvent;
class QKeyEvent;
class ParticleSystem;
class Renderer;

class ParticleApp
{
public:
    ParticleApp();
    ~ParticleApp();

    void render();
    void tick(float secs);

    void mousePressed(QMouseEvent *e, float x, float y);
    void mouseReleased(QMouseEvent *e, float x, float y);
    void mouseMoved(QMouseEvent *e, float deltaX, float deltaY);
    void mouseScrolled(QWheelEvent *e);

    void keyPressed(QKeyEvent *e);
    void keyReleased(QKeyEvent *e);

    void resize(int w, int h);

private:
    void makeInitScene();

    ParticleSystem *m_particleSystem;
    Renderer *m_renderer;

    bool m_mouseDownL;
    bool m_mouseDownR;

    bool m_fluidEmmiterOn;
    float m_timer;
};

#endif // PARTICLEAPP_H
