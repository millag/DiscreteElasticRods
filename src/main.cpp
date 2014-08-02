#include <QApplication>
#include "MainWindow.h"
#include <cmath>

#include <vector>
#include <iostream>
#include "Types.h"

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

//    std::vector<mg::Vec4D> verts(3);

//    for (unsigned i = 0; i < verts.size(); ++i)
//    {
//        std::cout << verts[i] << std::endl;
//    }
}
