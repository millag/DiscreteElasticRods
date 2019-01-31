#include "MainWindow.h"
#include "ui_MainWindow.h"
#include "Exporter.h"
#include <QFileDialog>

MainWindow::MainWindow( QWidget *parent ):
    QMainWindow( parent )
{
//	 setup the user interface
	m_ui = std::make_unique<Ui::MainWindow>();
	m_ui->setupUi( this );
//	 create our GL window for drawing the scene
	m_gl = new GLViewport( this );
//	 add glWindow to the UI
	m_ui->s_mainWindowGridLayout->addWidget( m_gl, 0, 0, 9, 1 );

	m_exportDir = "animation";
	m_gl->setSelection( false );
	m_animationBuffer.setCapacity( 500 );

//	wire up UI component slots
	connect( m_ui->m_selected, SIGNAL(currentIndexChanged(int)), this, SLOT(selectRenderObject(int)) );

	connect( m_ui->m_bendStiffness, SIGNAL(valueChanged(double)), this, SLOT(setBendingStiffness(double)) );
	connect( m_ui->m_twistStiffness, SIGNAL(valueChanged(double)), this, SLOT(setTwistingStiffness(double)) );
	connect( m_ui->m_maxForce, SIGNAL(valueChanged(double)), this, SLOT(setMaxElasticForce(double)) );
	connect( m_ui->m_drag, SIGNAL(valueChanged(double)), this, SLOT(setDrag(double)) );
	connect( m_ui->m_pbdIter, SIGNAL(valueChanged(int)), this, SLOT(setPBDIter(int)) );

	connect( m_ui->m_minimizationMethod, SIGNAL(currentIndexChanged(int)), this, SLOT(selectMinimizationMethod(int)) );
	connect( m_ui->m_minTolerance, SIGNAL(valueChanged(double)), this, SLOT(setMinimizationTolerance(double)) );
	connect( m_ui->m_minMaxIter, SIGNAL(valueChanged(int)), this, SLOT(setMinimizationMaxIter(int)) );

	connect( m_ui->m_collisions, SIGNAL(clicked(bool)), this, SLOT(toggleCollisions(bool)) );
	connect( m_ui->m_selfInteractions, SIGNAL(clicked(bool)), this, SLOT(toggleSelfInterations(bool)) );
	connect( m_ui->m_stiction, SIGNAL(valueChanged(double)), this, SLOT(setSelfStiction(double)) );
	connect( m_ui->m_repulsion, SIGNAL(valueChanged(double)), this, SLOT(setSelfRepusion(double)) );

	connect( m_ui->m_timerUpdate, SIGNAL(valueChanged(int)), this, SLOT(setTimerUpdateDuration(int)) );
	connect( m_ui->m_simBtn, SIGNAL(clicked(bool)), this, SLOT(toggleSim(bool)) );
	connect( m_ui->m_stepForward, SIGNAL(clicked()), this, SLOT(updateEvent()) );
	connect( m_ui->m_recordBtn, SIGNAL(clicked(bool)), this, SLOT(toggleRecord(bool)) );

	connect( m_ui->m_selectDirBtn, SIGNAL(clicked()), this, SLOT(selectExportDirectory()) );
	connect( m_ui->m_expBtn, SIGNAL(clicked()), this, SLOT(exportSim()) );

	connect( &m_updateTimer, SIGNAL(timeout()), this, SLOT(updateEvent()) );
	connect( &m_recordTimer, SIGNAL(timeout()), this, SLOT(recordEvent()) );

	toggleSim( m_ui->m_simBtn->isChecked() );
}

MainWindow::~MainWindow() = default;

void MainWindow::setScene( Scene& scene )
{
	m_scene = &scene;
	m_gl->setScene( m_scene );
	updateUI();
}

void MainWindow::toggleSim(bool s)
{
	if(!s)
	{
		m_updateTimer.stop();
		return;
	}

	m_updateTimer.start(m_ui->m_timerUpdate->value());
}

void MainWindow::updateEvent()
{
	if ( m_selectedObject )
	{
		m_selectedObject->setTransform( m_gl->getSelectionTransform() );
	}

	for ( auto i = 0; i < m_ui->m_simIter->value(); ++i )
	{
		m_scene->update( m_ui->m_timeStep->value() );
	}

	m_gl->update();
}

void MainWindow::toggleRecord(bool s)
{
	if(!s)
	{
		m_recordTimer.stop();
		return;
	}

	auto hair = m_scene->getHairById(0);
	assert( hair != nullptr );
	assert( hair->getSkin() != nullptr );

	if ( m_animationBuffer.size() )
	{
		auto object = m_scene->getRenderObjects()[hair->getSkin()->getId()];
		object->setTransform( m_animationBuffer.getFrame( 0 ) );
		m_gl->setSelectionTransform( m_animationBuffer.getFrame( 0 ) );
	}

	m_animationBuffer.clear();
	m_animationBuffer.saveHairState( *hair );
	m_recordTimer.start( m_ui->m_timerUpdate->value() );
}

void MainWindow::recordEvent()
{
	auto hair = m_scene->getHairById(0);
	assert( hair != nullptr );
	assert( hair->getSkin() != nullptr );

	std::cout << "Recording frame: " << m_animationBuffer.size() << std::endl;

	const auto transform = m_gl->getSelectionTransform();
	auto object = m_scene->getRenderObjects()[hair->getSkin()->getId()];

	m_animationBuffer.saveFrame( transform );
	object->setTransform( transform );
}

void MainWindow::setTimerUpdateDuration( int ms )
{
	UNUSED_VALUE( ms );
	m_updateTimer.setInterval( m_ui->m_timerUpdate->value() );
}

void MainWindow::selectRenderObject( int idx )
{
	idx = idx - 1;
	if ( idx < 0 )
	{
		m_selectedObject.reset();
		m_gl->setSelection( false );
		return;
	}

	m_selectedObject = m_scene->getRenderObjects()[idx];
	m_gl->setSelection( true );
	m_gl->setSelectionTransform( m_selectedObject->getTransform() );
}

void MainWindow::selectExportDirectory()
{
	QString dir = QFileDialog::getExistingDirectory( this,
	                                                 "Select Directory",
	                                                 m_exportDir,
	                                                 QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks );

	if ( !dir.isEmpty() )
	{
		m_exportDir = dir;
	}
}

void MainWindow::exportSim()
{
	if ( m_ui->m_filePrefix->text().isEmpty() )
	{
		std::cerr << "Specify export file name" << std::endl;
		return;
	}

	if ( m_animationBuffer.size() == 0 )
	{
		std::cerr << "No data to export" << std::endl;
		return;
	}

	auto hair = m_scene->getHairById(0);
	assert( hair != nullptr );
	assert( hair->getSkin() != nullptr );

	QDir parentDir( m_exportDir );
	m_exportGeoDir = m_ui->m_filePrefix->text().append( "geo" );
	if ( !parentDir.exists( m_exportGeoDir ) )
	{
		parentDir.mkdir( m_exportGeoDir);
		m_exportGeoDir = parentDir.filePath( m_exportGeoDir );
	}

	m_exportCurvDir = m_ui->m_filePrefix->text().append( "curv" );
	if ( !parentDir.exists (m_exportCurvDir ) )
	{
		parentDir.mkdir( m_exportCurvDir );
		m_exportCurvDir = parentDir.filePath( m_exportCurvDir );
	}

	m_animationBuffer.restoreHairState( *m_scene->getHairById(0) );
	auto object = m_scene->getRenderObjects()[ hair->getSkin()->getId() ];
	for ( auto frame = 0u; frame < m_animationBuffer.size(); ++frame )
	{
		std::cout << "Exporting frame: " << frame << std::endl;

		object->setTransform( m_animationBuffer.getFrame( frame ) );
		for ( auto i = 0; i < m_ui->m_simIter->value(); ++i )
		{
			m_scene->update( m_ui->m_timeStep->value() );
		}

		exportFrame( frame );
	}
}

void MainWindow::exportFrame(unsigned frame)
{
	auto hair = m_scene->getHairById(0);
	assert( hair != nullptr );
	assert( hair->getSkin() != nullptr );

	QDir dir( m_exportGeoDir );
	QString filename = dir.filePath( m_ui->m_filePrefix->text().append( "geo_%1.obj" ).arg( frame ) );
	Exporter exporter;

	exporter.exportGeometry( filename.toLocal8Bit(), *hair->getSkin() );

	dir.setPath( m_exportCurvDir );
	filename = dir.filePath( m_ui->m_filePrefix->text().append( "curv_%1.obj" ).arg( frame ) );
	exporter.exportCurves( filename.toLocal8Bit(), hair->getStrands() );
}

void MainWindow::updateUI()
{
	m_ui->m_minimizationMethod->clear();
	m_ui->m_minimizationMethod->addItem( QString( "NONE" ) );
	m_ui->m_minimizationMethod->addItem( QString( "NEWTON" ) );
	m_ui->m_minimizationMethod->addItem( QString( "BFGS" ) );
	m_ui->m_minimizationMethod->addItem( QString( "BFGS NUMERIC" ) );

	m_ui->m_selected->clear();
	m_ui->m_selected->addItem( QString( "None" ) );

	if ( m_scene )
	{
		for ( auto i = 0u; i < m_scene->getRenderObjects().size(); ++i )
		{
			m_ui->m_selected->addItem( QString( "RenderObject_%1" ).arg( i ) );
		}

		auto hair = m_scene->getHairById(0);
		if ( hair )
		{
			m_ui->m_minimizationMethod->setCurrentIndex( static_cast<int>( hair->m_params.m_rodParams.m_strategy ) );
			m_ui->m_minTolerance->setValue( hair->m_params.m_rodParams.m_tolerance );
			m_ui->m_minMaxIter->setValue( hair->m_params.m_rodParams.m_maxIter );

			m_ui->m_bendStiffness->setValue( hair->m_params.m_rodParams.m_B( 0, 0 ) );
			m_ui->m_twistStiffness->setValue( hair->m_params.m_rodParams.m_beta );
			m_ui->m_maxForce->setValue( hair->m_params.m_rodParams.m_maxElasticForce );
			m_ui->m_drag->setValue( hair->m_params.m_drag);
			m_ui->m_pbdIter->setValue( hair->m_params.m_pbdIter );

			m_ui->m_collisions->setChecked( hair->m_params.m_resolveCollisions );
			m_ui->m_selfInteractions->setChecked( hair->m_params.m_resolveSelfInterations );
			m_ui->m_stiction->setValue( hair->m_params.m_selfStiction );
			m_ui->m_repulsion->setValue( hair->m_params.m_selfRepulsion );
		}
	}
}

void MainWindow::selectMinimizationMethod(int idx)
{
	assert( m_scene != nullptr );
	auto hair = m_scene->getHairById(0);
	assert( hair != nullptr );

	hair->m_params.m_rodParams.m_strategy = ( idx < static_cast<int>( ElasticRodParams::MinimizationStrategy::SENTINEL ) )?
	                                            static_cast<ElasticRodParams::MinimizationStrategy>( idx ) :
	                                            ElasticRodParams::MinimizationStrategy::SENTINEL;
}

void MainWindow::setMinimizationTolerance(double val)
{
	assert( m_scene != nullptr );
	auto hair = m_scene->getHairById(0);
	assert( hair != nullptr );

	hair->m_params.m_rodParams.m_tolerance = val;
}

void MainWindow::setMinimizationMaxIter(int val)
{
	assert( m_scene != nullptr );
	auto hair = m_scene->getHairById(0);
	assert( hair != nullptr );

	hair->m_params.m_rodParams.m_maxIter = val;
}

void MainWindow::setBendingStiffness(double val)
{
	assert( m_scene != nullptr );
	auto hair = m_scene->getHairById(0);
	assert( hair != nullptr );

	hair->m_params.m_rodParams.setBendStiffness(val);
}

void MainWindow::setTwistingStiffness(double val)
{
	assert( m_scene != nullptr );
	auto hair = m_scene->getHairById(0);
	assert( hair != nullptr );

	hair->m_params.m_rodParams.setTwistStiffness(val);
}

void MainWindow::setMaxElasticForce(double val)
{
	assert( m_scene != nullptr );
	auto hair = m_scene->getHairById(0);
	assert( hair != nullptr );

	hair->m_params.m_rodParams.m_maxElasticForce = val;
}

void MainWindow::setDrag(double val)
{
	assert( m_scene != nullptr );
	auto hair = m_scene->getHairById(0);
	assert( hair != nullptr );

	hair->m_params.m_drag = val;
}

void MainWindow::setPBDIter(int val)
{
	assert( m_scene != nullptr );
	auto hair = m_scene->getHairById(0);
	assert( hair != nullptr );

	hair->m_params.m_pbdIter = val;
}

void MainWindow::toggleCollisions(bool val)
{
	assert( m_scene != nullptr );
	auto hair = m_scene->getHairById(0);
	assert( hair != nullptr );

	hair->m_params.m_resolveCollisions = val;
}

void MainWindow::toggleSelfInterations(bool val)
{
	assert( m_scene != nullptr );
	auto hair = m_scene->getHairById(0);
	assert( hair != nullptr );

	hair->m_params.m_resolveSelfInterations = val;
}

void MainWindow::setSelfStiction(double val)
{
	assert( m_scene != nullptr );
	auto hair = m_scene->getHairById(0);
	assert( hair != nullptr );

	hair->m_params.m_selfStiction = val;
}

void MainWindow::setSelfRepusion(double val)
{
	assert( m_scene != nullptr );
	auto hair = m_scene->getHairById(0);
	assert( hair != nullptr );

	hair->m_params.m_selfRepulsion = val;
}
