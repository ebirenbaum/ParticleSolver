#include "view.h"
#include <QMainWindow>
#include <QApplication>
#include <QKeyEvent>
#include <iostream>
#include "particleapp.h"

View::View(QGLFormat format, QWidget *parent)
    : QGLWidget(format, parent),
      m_width(parent->width()),
      m_height(parent->height())
{
    // View needs all mouse move events, not just mouse drag events
    setMouseTracking(true);

    // Hide the cursor since this is a fullscreen app
    setCursor(Qt::BlankCursor);

    // View needs keyboard focus
    setFocusPolicy(Qt::StrongFocus);

    // The game loop is implemented using a timer
    connect(&timer, SIGNAL(timeout()), this, SLOT(tick()));

    m_app = NULL;
}

View::~View()
{
    if (m_app)
        delete m_app;
}

void View::initializeGL()
{
    // All OpenGL initialization *MUST* be done during or after this
    // method. Before this method is called, there is no active OpenGL
    // context and all OpenGL calls have no effect.

    glewExperimental = GL_TRUE;
    GLenum err = glewInit();
    glGetError(); // Clear errors after call to glewInit
    if (GLEW_OK != err)
    {
      // Problem: glewInit failed, something is seriously wrong.
      fprintf(stderr, "Error initializing glew: %s\n", glewGetErrorString(err));
    }

    // init the App object.
    m_app = new ParticleApp();

    // Enable depth testing, so that objects are occluded based on depth instead of drawing order.
    glEnable(GL_DEPTH_TEST);

    // Move the polygons back a bit so lines are still drawn even though they are coplanar with the
    // polygons they came from, which will be drawn before them.
    glEnable(GL_POLYGON_OFFSET_LINE);
    glPolygonOffset(-1, -1);

    // Enable back-face culling, meaning only the front side of every face is rendered.
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);

    // Specify that the front face is represented by vertices in counterclockwise order (this is
    // the default).
    glFrontFace(GL_CCW);

    // Start a timer that will try to get 60 frames per second (the actual
    // frame rate depends on the operating system and other running programs)
    time.start();
//    timer.start();
    timer.start(1000 / 60);

    // Center the mouse, which is explained more in mouseMoveEvent() below.
    // This needs to be done here because the mouse may be initially outside
    // the fullscreen window and will not automatically receive mouse move
    // events. This occurs if there are two monitors and the mouse is on the
    // secondary monitor.
//    QCursor::setPos(mapToGlobal(QPoint(width() / 2, height() / 2)));
}

void View::paintGL()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glClearColor(0.f, 0.f, 0.f, 1.f);

    QString title = "PARTICLES TEST      FPS: " + QString::number((int) fps);
    emit changeTitle(title);

    // TODO: call your game rendering code here
    m_app->render();
}

void View::resizeGL(int w, int h)
{
    glViewport(0, 0, w, h);

    m_width = w;
    m_height = h;

    m_app->resize(w, h);
}

void View::mousePressEvent(QMouseEvent *e)
{
    float x = (e->x() * 2.f / m_width) - 1.f;
    float y = 1.f - (e->y() * 2.f / m_height);
    m_app->mousePressed(e, x, y);
}

void View::mouseMoveEvent(QMouseEvent *event)
{
    // This starter code implements mouse capture, which gives the change in
    // mouse position since the last mouse movement. The mouse needs to be
    // recentered after every movement because it might otherwise run into
    // the edge of the screen, which would stop the user from moving further
    // in that direction. Note that it is important to check that deltaX and
    // deltaY are not zero before recentering the mouse, otherwise there will
    // be an infinite loop of mouse move events.
    int halfWidthI = width() / 2;
    int halfHightI = height() / 2;

    int deltaXI = event->x() - halfWidthI;
    int deltaYI = event->y() - halfHightI;
    if (!deltaXI && !deltaYI) return;

    QCursor::setPos(mapToGlobal(QPoint(halfWidthI, halfHightI)));

    // sets mouse deltas between -1 and 1
    float halfWidth = width() * .5f;
    float halfHeight = height() * .5f;

    float deltaX = (event->x() - halfWidth) / halfWidth;
    float deltaY = (event->y() - halfHeight) / halfHeight;

    m_app->mouseMoved(event, deltaX, deltaY);

}

void View::mouseReleaseEvent(QMouseEvent *e)
{
    float x = (e->x() * 2.f / m_width) - 1.f;
    float y = 1.f - (e->y() * 2.f / m_height);
    m_app->mouseReleased(e, x, y);
}

void View::wheelEvent(QWheelEvent *e)
{
    m_app->mouseScrolled(e);
}

void View::keyPressEvent(QKeyEvent *event)
{
    if (event->key() == Qt::Key_Escape) QApplication::quit();

    // TODO: Handle keyboard presses here
    m_app->keyPressed(event);
}

void View::keyReleaseEvent(QKeyEvent *e)
{
    m_app->keyReleased(e);
}

void View::tick()
{
    // Get the number of seconds since the last tick (variable update rate)
    float seconds = time.restart() * 0.001f;
    fps = .02f / seconds + .98f * fps;

    // TODO: Implement the game update here
    m_app->tick(seconds);

    // don't show cursor
    QApplication::setOverrideCursor(Qt::CrossCursor);

    // Flag this view for repainting (Qt will call paintGL() soon after)
    update();
}
