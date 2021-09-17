#include "MainWindow.h"
#include "Scene.h"
#include "SceneLoader.h"
#include "BasicParser.h"

#include <QApplication>

void InitGLDefaultSurfaceFormat()
{
    auto format = QSurfaceFormat::defaultFormat();
    format.setProfile( QSurfaceFormat::CoreProfile );
    format.setSwapBehavior( QSurfaceFormat::DoubleBuffer );
    format.setDepthBufferSize( 24 );
    format.setStencilBufferSize( 8 );
    format.setSamples( 4 );
    QSurfaceFormat::setDefaultFormat( format );
}

void test() {
    char szOrbits[] = "365.24 29.53";
//    char* pEnd;
    qDebug() << BasicParser::parseDouble(szOrbits);
//    double d1, d2;
//    d1 = std::strtod (szOrbits, &pEnd);
//    d2 = std::strtod (pEnd, NULL);
//    qDebug("The moon completes %.2f orbits per Earth year.\n", d1/d2);

}

int main(int argc, char *argv[])
{
//    test();
//    return 0;

    qInfo() << "Multi-threading" << ( ( MULTI_THREADING_ON )? "enabled" : "disabled" );

    try
    {
        InitGLDefaultSurfaceFormat();
        QApplication app(argc, argv);

        // create our scene
        Scene scene;
        SceneLoader loader;
        // loader.loadTestScene( scene );
        loader.loadScene( "assets/scene1_long_curly.mg", scene );
        scene.initialize();

        MainWindow window;
        window.setScene( scene );
        window.show();
        return app.exec();
    }
    catch( std::exception& e )
    {
        qFatal("Fatal error: %s", e.what());
    }

    return 0;
}
