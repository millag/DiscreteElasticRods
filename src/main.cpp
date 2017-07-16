#include "MainWindow.h"

#include <QApplication>


void InitDefaultGLSurfaceFormat()
{
    auto format = QSurfaceFormat::defaultFormat();
    format.setMajorVersion( 4 );
    format.setMinorVersion( 1 );
    format.setProfile( QSurfaceFormat::CoreProfile );
    format.setSwapBehavior( QSurfaceFormat::DoubleBuffer );
    format.setDepthBufferSize( 24 );
    format.setStencilBufferSize( 8 );
    format.setSamples( 4 );
    QSurfaceFormat::setDefaultFormat( format );
}

int main(int argc, char *argv[])
{
    InitDefaultGLSurfaceFormat();

    QApplication app(argc, argv);

    try
    {
        MainWindow w;
        w.show();
        return app.exec();
    }
    catch( std::exception& e )
    {
        qInfo() << "Fatal error: "
                << e.what();
    }

    return 1;
}
