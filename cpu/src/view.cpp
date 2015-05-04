#include "view.h"
#include <QApplication>
#include <QKeyEvent>

View::View(QWidget *parent) : QGLWidget(parent)
{
    // View needs all mouse move events, not just mouse drag events
    setMouseTracking(true);

    // Hide the cursor since this is a fullscreen app
//    setCursor(Qt::BlankCursor);

    // View needs keyboard focus
    setFocusPolicy(Qt::StrongFocus);

    // The game loop is implemented using a timer
    connect(&timer, SIGNAL(timeout()), this, SLOT(tick()));

    fps = 60;
    scale = 10;
    tickTime = 0.0;
    timestepMode = true;
    current = FLUID_TEST;
}

View::~View()
{
}

void View::initializeGL()
{
    // All OpenGL initialization *MUST* be done during or after this
    // method. Before this method is called, there is no active OpenGL
    // context and all OpenGL calls have no effect.

    // Start a timer that will try to get 60 frames per second (the actual
    // frame rate depends on the operating system and other running programs)
    time.start();
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
    glClearColor(0,0,0,1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-scale, scale, -scale, scale, -1, 1);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    sim.draw();

    glColor3f(1,1,1);
    renderText(10, 20, "FPS: " + QString::number((int) (fps)), this->font());
    renderText(10, 40, "# Particles: " + QString::number(sim.getNumParticles()), this->font());
    renderText(10, 60, "Kinetic Energy: " + QString::number(sim.getKineticEnergy()), this->font());
}

void View::resizeGL(int w, int h)
{
    glViewport(0, 0, w, h);
    dimensions.x = w;
    dimensions.y = h;

    sim.resize(glm::ivec2(scale, scale));
}



void View::mousePressEvent(QMouseEvent *event)
{
    glm::dvec2 screen(event->x(), height() - event->y());
    glm::dvec2 world = (scale * 2.) * (screen / glm::dvec2(width(), height())) - glm::dvec2((double)scale, (double)scale);
    sim.mousePressed(world);
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
    int deltaX = event->x() - width() / 2;
    int deltaY = event->y() - height() / 2;
    if (!deltaX && !deltaY) return;
//    QCursor::setPos(mapToGlobal(QPoint(width() / 2, height() / 2)));

    // TODO: Handle mouse movements here
}

void View::mouseReleaseEvent(QMouseEvent *event)
{

}

void View::wheelEvent(QWheelEvent *event)
{
    if (event->delta() > 0) {
        scale /= 1.2;
    } else {
        scale *= 1.2;
    }
    sim.resize(glm::ivec2(scale, scale));
}

void View::keyPressEvent(QKeyEvent *event)
{
    if (event->key() == Qt::Key_Escape) QApplication::quit();
    if (event->key() == Qt::Key_R) sim.init(current);
    if (event->key() == Qt::Key_T) timestepMode = !timestepMode;
    if (event->key() == Qt::Key_Space) tickTime = .01;
    if (event->key() == Qt::Key_Backspace) tickTime = -.01;
    if (event->key() == Qt::Key_C) sim.debug = !sim.debug;

    if (event->key() == Qt::Key_1) {
        current = GRANULAR_TEST;
        sim.init(current);
    } else if (event->key() == Qt::Key_2) {
        current = STACKS_TEST;
        sim.init(current);
    } else if (event->key() == Qt::Key_3) {
        current = WALL_TEST;
        sim.init(current);
    } else if (event->key() == Qt::Key_4) {
        current = PENDULUM_TEST;
        sim.init(current);
    } else if (event->key() == Qt::Key_5) {
        current = ROPE_TEST;
        sim.init(current);
    } else if (event->key() == Qt::Key_6) {
        current = FLUID_TEST;
        sim.init(current);
    } else if (event->key() == Qt::Key_7) {
        current = FLUID_SOLID_TEST;
        sim.init(current);
    } else if (event->key() == Qt::Key_8) {
        current = GAS_TEST;
        sim.init(current);
    } else if (event->key() == Qt::Key_9) {
        current = FRICTION_TEST;
        sim.init(current);
    } else if (event->key() == Qt::Key_0) {
        current = WATER_BALLOON_TEST;
        sim.init(current);
    } else if (event->key() == Qt::Key_N) {
        current = CRADLE_TEST;
        sim.init(current);
    } else if (event->key() == Qt::Key_Period) {
        current = SDF_TEST;
        sim.init(current);
    }
}

void View::keyReleaseEvent(QKeyEvent *event)
{
}

void View::tick()
{
    // Get the number of seconds since the last tick (variable update rate)
    double seconds = time.restart() * 0.001;
    fps = .02 / seconds + .98 * fps;

    if (timestepMode) {
        if (tickTime != 0.0) {
            sim.tick(tickTime);
            tickTime = 0.0;
        }
    } else {
        sim.tick(.01f);
    }

    // Flag this view for repainting (Qt will call paintGL() soon after)
    update();
}
