#include <QtGui/QApplication>
#include "mainwindow.h"

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow w;

    // We cannot use w.showFullscreen() here because on Linux that creates the
    // window behind all other windows, so we have to set it to fullscreen after
    // it has been shown.

    w.setFixedSize(800, 800);
    w.show();
//    w.setWindowState(w.windowState() | Qt::WindowFullScreen); // Comment out this line to have a windowed 800x600 game on startup.

    return a.exec();
}

