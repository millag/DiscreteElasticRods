#include <QApplication>
#include "MainWindow.h"
#include <cmath>

/* this code runs the basic main window and is created by the Qt Creator app */
int main(int argc, char *argv[])
{
    // make an instance of the QApplication
    QApplication a(argc, argv);
    // Create a new MainWindow
    MainWindow w;
    // show it
    w.show();
    // hand control over to Qt framework
    return a.exec();

    // this is bullshit
//    ngl::Mat4 m1;
//    m1.rotateZ(90);
//    ngl::Mat4 m2;
//    m2.rotateX(90);
//    ngl::Mat4 m = m1 * m2;
//    m1 *= m2;
//    std::cout << "m = \n" << m <<std::endl;
//    std::cout << "m1 = \n" << m1 <<std::endl;
}
