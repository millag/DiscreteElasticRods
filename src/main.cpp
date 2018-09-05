#include "MainWindow.h"
#include "Scene.h"
#include "SceneLoader.h"

#include <QApplication>

void InitializeGLDefaultSurfaceFormat()
{
	auto format = QSurfaceFormat::defaultFormat();
	format.setProfile( QSurfaceFormat::CoreProfile );
	format.setSwapBehavior( QSurfaceFormat::DoubleBuffer );
	format.setDepthBufferSize( 24 );
	format.setStencilBufferSize( 8 );
	format.setSamples( 4 );
	QSurfaceFormat::setDefaultFormat( format );
}

int main(int argc, char *argv[])
{
	InitializeGLDefaultSurfaceFormat();

	QApplication app(argc, argv);

	if ( MULTI_THREAD )
	{
		qDebug() << "Multi-threading enabled";
	}
	else
	{
		qDebug() << "Multi-threading disabled";
	}

	try
	{
//		create our scene
		Scene scene;
		SceneLoader loader;
//		loader.loadTestScene( scene );
		loader.loadScene( "assets/scene1_long_curly.mg", scene );
		scene.initialize();

		MainWindow window;
		window.setScene( scene );
		window.show();
		return app.exec();
	}
	catch( std::exception& e )
	{
		qInfo() << "Fatal error: "
		        << e.what();
	}

	return 0;
}
