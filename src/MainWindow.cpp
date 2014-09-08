#include "MainWindow.h"
#include "ui_MainWindow.h"
#include "SceneLoader.h"
#include "Exporter.h"


SceneLoader loader;
Exporter exporter;

//----------------------------------------------------------------------------------------------------------------------
MainWindow::MainWindow(QWidget *parent) :QMainWindow(parent), m_ui(new Ui::MainWindow)
{
    // create our scene
    m_scene = loader.loadScene("assets/scene2_long_wavy.mg");
//    m_scene = loader.loadTestScene();
    m_scene->initialize();

    // create an openGL format and pass to the new GLWidget
    QGLFormat format;
    format.setVersion(3, 2);
    format.setProfile(QGLFormat::CoreProfile);

	// setup the user interface
	m_ui->setupUi(this);
    // create our GL window for drawing the scene
    m_gl = new GLWindow(format, this);
    m_gl->setScene(m_scene);
    // add glWindow to the UI
    m_ui->s_mainWindowGridLayout->addWidget(m_gl, 0, 0, 6, 1);

    m_exportDir = "animation";
    m_selectedObject = NULL;
    m_gl->setSelection(false);
    m_animationBuffer.setCapacity(500);

    populateUI();
    // now we wire up the UI components to the slots
    connect(m_ui->m_selected, SIGNAL(currentIndexChanged(int)), this, SLOT(selectRenderObject(int)));
    connect(m_ui->m_minimizationMethod, SIGNAL(currentIndexChanged(int)), this, SLOT(selectMinimizationMethod(int)));
    connect(m_ui->m_timerUpdate, SIGNAL(valueChanged(int)), m_gl, SLOT(setTimerUpdateDuration(int)));
    connect(m_ui->m_simBtn, SIGNAL(clicked(bool)), this, SLOT(toggleSim(bool)));
    connect(m_ui->m_recordBtn, SIGNAL(clicked(bool)), this, SLOT(toggleRecord(bool)));
    connect(m_ui->m_selectDirBtn, SIGNAL(clicked()), this, SLOT(selectExportDirectory()));
    connect(m_ui->m_expBtn, SIGNAL(clicked()), this, SLOT(exportSim()));

    connect(&m_updateTimer, SIGNAL(timeout()), this, SLOT(updateEvent()));
    connect(&m_recordTimer, SIGNAL(timeout()), this, SLOT(recordEvent()));

    m_ui->m_minimizationMethod->setCurrentIndex(m_scene->getHairById(0)->m_params->m_rodParams->m_strategy);
    m_chronometer.start();
    toggleSim(m_ui->m_simBtn->isChecked());
}

MainWindow::~MainWindow()
{
    delete m_scene;
    delete m_ui;
}

void MainWindow::toggleSim(bool s)
{
    if(!s)
    {
        m_updateTimer.stop();
        m_chronometer.restart();
        return;
    }

    m_updateTimer.start(m_ui->m_timerUpdate->value());
    m_chronometer.restart();
}

void MainWindow::updateEvent()
{
    std::cout << "TIME: " << m_chronometer.elapsed() << std::endl;
    std::cout << "FPS: " << (float)mg::SEC / m_chronometer.restart() << std::endl;

    if (m_selectedObject != NULL)
    {
        m_selectedObject->setTransform(m_gl->getSelectionTransform());
    }

    for (int i = 0; i < m_ui->m_simIter->value(); ++i)
    {
        m_scene->update(m_ui->m_timeStep->value());
    }

    m_gl->updateGL();
}

void MainWindow::toggleRecord(bool s)
{
    if(!s)
    {
        m_recordTimer.stop();
        return;
    }

    if (m_animationBuffer.size())
    {
        RenderObject* object = m_scene->getRenderObjects()[ m_scene->getHairById(0)->m_object->getId() ];
        object->setTransform(m_animationBuffer.getFrame(0));
        m_gl->setSelectionTransform(m_animationBuffer.getFrame(0));
    }

    m_animationBuffer.clear();
    m_animationBuffer.saveHairState(m_scene->getHairById(0));
    m_recordTimer.start(m_ui->m_timerUpdate->value());
}

void MainWindow::recordEvent()
{
    if (m_scene->getHairById(0) == NULL)
    {
        std::cerr << "No hair object found. Skipping..." << std::endl;
        return;
    }
    std::cout << "Recording frame: " << m_animationBuffer.size() << std::endl;

    mg::Matrix4D transform = m_gl->getSelectionTransform();
    RenderObject* object = m_scene->getRenderObjects()[ m_scene->getHairById(0)->m_object->getId() ];

    m_animationBuffer.saveFrame(transform);
    object->setTransform(transform);
}

void MainWindow::setTimerUpdateDuration(int ms)
{
    m_updateTimer.setInterval(m_ui->m_timerUpdate->value());
}

void MainWindow::selectRenderObject(int index)
{
    index = index - 1;
    if (index < 0)
    {
        m_selectedObject = NULL;
        m_gl->setSelection(false);
        return;
    }
    m_selectedObject = m_scene->getRenderObjects()[index];
    m_gl->setSelection(true);
    m_gl->setSelectionTransform(m_selectedObject->getTransform());
}

void MainWindow::selectMinimizationMethod(int index)
{
    if (m_scene->getHairById(0) == NULL)
    {
        std::cerr << "No hair object found. Skipping..." << std::endl;
        return;
    }

    switch (index)
    {
        case 0:
            m_scene->getHairById(0)->m_params->m_rodParams->m_strategy = ElasticRodParams::NONE;
            break;
        case 1:
            m_scene->getHairById(0)->m_params->m_rodParams->m_strategy = ElasticRodParams::NEWTON;
            break;
        case 2:
            m_scene->getHairById(0)->m_params->m_rodParams->m_strategy = ElasticRodParams::BFGS;
            break;
        case 3:
            m_scene->getHairById(0)->m_params->m_rodParams->m_strategy = ElasticRodParams::BFGS_NUMERIC;
            break;
        default:
            break;
    }
}

void MainWindow::selectExportDirectory()
{
    QString dir = QFileDialog::getExistingDirectory(this,
                                                    "Select Directory",
                                                    m_exportDir,
                                                    QFileDialog::ShowDirsOnly
                                                    | QFileDialog::DontResolveSymlinks);

    if (dir.isEmpty())
    {
        return;
    }
    m_exportDir = dir;
}

void MainWindow::exportSim()
{
    if (m_ui->m_filePrefix->text().isEmpty())
    {
        std::cerr << "Specify export file name" << std::endl;
        return;
    }

    if(m_animationBuffer.size() == 0)
    {
        std::cerr << "No data to export" << std::endl;
        return;
    }

    QDir parentDir(m_exportDir);
    m_exportGeoDir = m_ui->m_filePrefix->text().append("geo");
    if (!parentDir.exists(m_exportGeoDir))
    {
        parentDir.mkdir(m_exportGeoDir);
        m_exportGeoDir = parentDir.filePath(m_exportGeoDir);
    }

    m_exportCurvDir = m_ui->m_filePrefix->text().append("curv");
    if (!parentDir.exists(m_exportCurvDir))
    {
        parentDir.mkdir(m_exportCurvDir);
        m_exportCurvDir = parentDir.filePath(m_exportCurvDir);
    }

    m_animationBuffer.restoreHairState(m_scene->getHairById(0));
    RenderObject* object = m_scene->getRenderObjects()[ m_scene->getHairById(0)->m_object->getId() ];
    for (unsigned frame = 0; frame < m_animationBuffer.size(); ++frame)
    {
        std::cout << "Exporting frame: " << frame << std::endl;

        object->setTransform(m_animationBuffer.getFrame(frame));

        for (int i = 0; i < m_ui->m_simIter->value(); ++i)
        {
            m_scene->update(m_ui->m_timeStep->value());
        }

        exportFrame(frame);
    }
}

void MainWindow::exportFrame(unsigned frame)
{
    QDir dir(m_exportGeoDir);
    QString filename = dir.filePath(m_ui->m_filePrefix->text().append("geo_%1.obj").arg(frame));
    exporter.exportGeometry(filename.toLocal8Bit(), *m_scene->getHairById(0)->m_object);

    dir.setPath(m_exportCurvDir);
    filename = dir.filePath(m_ui->m_filePrefix->text().append("curv_%1.obj").arg(frame));
    exporter.exportCurves(filename.toLocal8Bit(), m_scene->getHairById(0)->m_strands);
}

void MainWindow::populateUI()
{
    m_ui->m_selected->addItem(QString("None"));
    for (unsigned i = 0; i < m_scene->getRenderObjects().size(); ++i)
    {
        QString label = QString("RenderObject_%1").arg(i);
        m_ui->m_selected->addItem(label);
    }

    m_ui->m_minimizationMethod->addItem(QString("NONE"));
    m_ui->m_minimizationMethod->addItem(QString("NEWTON"));
    m_ui->m_minimizationMethod->addItem(QString("BFGS"));
    m_ui->m_minimizationMethod->addItem(QString("BFGS NUMERIC"));
}
