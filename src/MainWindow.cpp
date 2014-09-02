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
    m_scene = loader.loadScene("assets/scene1.mg");
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
    m_ui->s_mainWindowGridLayout->addWidget(m_gl, 0, 0, 5, 1);

    m_exportDir = "animation";

    populateUI();

    // now we wire up the UI components to the slots
    connect(m_ui->m_selected, SIGNAL(currentIndexChanged(int)), this, SLOT(selectRenderObject(int)));
    connect(m_ui->m_timerUpdate, SIGNAL(valueChanged(int)), m_gl, SLOT(setTimerUpdateDuration(int)));
    connect(m_ui->m_simBtn, SIGNAL(clicked(bool)), this, SLOT(toggleSim(bool)));
    connect(m_ui->m_expBtn, SIGNAL(clicked(bool)), this, SLOT(exportSim(bool)));
    connect(m_ui->m_selectDirBtn, SIGNAL(clicked()), this, SLOT(selectExportDirectory()));

    connect(&m_updateTimer, SIGNAL(timeout()), this, SLOT(timerEvent()));

    m_chronometer.start();
    toggleSim(m_ui->m_simBtn->isChecked());
    exportSim(m_ui->m_expBtn->isChecked());
}

MainWindow::~MainWindow()
{
    delete m_scene;
    delete m_ui;
}

void MainWindow::toggleSim(bool s)
{
    if(s == true)
    {
        startSim();
    }
    else
    {
        stopSim();
    }
}

void MainWindow::startSim()
{
    m_updateTimer.start(m_ui->m_timerUpdate->value());
    m_chronometer.restart();
}

void MainWindow::stopSim()
{
    m_updateTimer.stop();
    m_chronometer.restart();
}

void MainWindow::timerEvent()
{
    std::cout << "TIME: " << m_chronometer.elapsed() << std::endl;
    std::cout << "FPS: " << (float)mg::SEC / m_chronometer.restart() << std::endl;

    for (int i = 0; i < m_ui->m_simIter->value(); ++i)
    {
        m_scene->update(m_ui->m_timeStep->value());
//        m_scene->update(0.01);
    }
    if (m_frame >= 0)
    {
        exportFrame();
    }
    m_gl->updateGL();
}

void MainWindow::setTimerUpdateDuration(int ms)
{
    m_updateTimer.setInterval(m_ui->m_timerUpdate->value());
}

void MainWindow::selectRenderObject(int index)
{
    index = index - 1;
    RenderObject* object = (index < 0)? NULL : m_scene->getRenderObjects()[index];
    m_gl->setSelectedObject(object);
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

void MainWindow::exportSim(bool v)
{
    if (!v)
    {
        m_frame = -1;
        return;
    }
    if (m_ui->m_filePrefix->text().isEmpty())
    {
        std::cerr << "Pls specify export file name" << std::endl;
        return;
    }

    m_frame = 0;

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
}

void MainWindow::exportFrame()
{
    QDir dir(m_exportGeoDir);
    QString filename = dir.filePath(m_ui->m_filePrefix->text().append("geo_%1.obj").arg(m_frame));
    exporter.exportGeometry(filename.toLocal8Bit(), *m_scene->getHair()->m_object);

    dir.setPath(m_exportCurvDir);
    filename = dir.filePath(m_ui->m_filePrefix->text().append("curv_%1.obj").arg(m_frame));
    exporter.exportCurves(filename.toLocal8Bit(), m_scene->getHair()->m_strands);
    ++m_frame;
}

void MainWindow::populateUI()
{
    m_ui->m_selected->addItem(QString("None"));
    for (unsigned i = 0; i < m_scene->getRenderObjects().size(); ++i)
    {
        QString label = QString("RenderObject_%1").arg(i);
        m_ui->m_selected->addItem(label);
    }
}
