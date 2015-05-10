#ifndef VIEW_H
#define VIEW_H

#include <GL/glew.h>
#include <qgl.h>
#include <QTime>
#include <QTimer>

class ParticleApp;

class View : public QGLWidget
{
    Q_OBJECT

public:
    View(QGLFormat format, QWidget *parent);
    ~View();

private:
    ParticleApp *m_app;

    QTime time;
    QTimer timer;
    float fps;
    float fpsTimer;
    float avgFps;
    int fpsTicks;

    void initializeGL();
    void paintGL();
    void resizeGL(int w, int h);

    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    void mouseReleaseEvent(QMouseEvent *event);

    void wheelEvent(QWheelEvent *event);

    void keyPressEvent(QKeyEvent *event);
    void keyReleaseEvent(QKeyEvent *event);

private slots:
    void tick();

signals:
    void changeTitle(const QString);

private:
    int m_width;
    int m_height;

};

#endif // VIEW_H

