#include <QApplication>
#include <iostream>
#include "mainwindow.h"

int main(int argc, char *argv[])
{
#ifdef CUDA
    QApplication a(argc, argv);
    MainWindow w;

    // We cannot use w.showFullscreen() here because on Linux that creates the
    // window behind all other windows, so we have to set it to fullscreen after
    // it has been shown. 
    w.show();
//    w.setWindowState(w.windowState() | Qt::WindowFullScreen); // Comment out this line to have a windowed 800x600 game on startup.

    return a.exec();
#else
    std::cout << "The CUDA directory required to run this program cannot be located.\nExiting..." << std::endl;
    return 0;
#endif

}

