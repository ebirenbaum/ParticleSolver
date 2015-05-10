#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

class View;

namespace Ui {
    class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private:
    Ui::MainWindow *ui;

private slots:
    void changeTitle(const QString &title);

};

#endif // MAINWINDOW_H

