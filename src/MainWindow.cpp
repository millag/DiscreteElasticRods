#include "MainWindow.h"
#include "ui_MainWindow.h"
#include "ElasticRod.h"

//----------------------------------------------------------------------------------------------------------------------
MainWindow::MainWindow(QWidget *parent) :QMainWindow(parent), m_ui(new Ui::MainWindow)
{
    // create our scene
    m_scene = new Scene();
    m_scene->initialize();

    // create an openGL format and pass to the new GLWidget
    QGLFormat format;
    format.setVersion(3,2);
    format.setProfile(QGLFormat::CoreProfile);


	// setup the user interface
	m_ui->setupUi(this);
    // create our GL window for drawing the scene
    m_gl = new  GLWindow(format, m_ui->m_timerValue->value(), this);
    m_gl->setScene(m_scene);
    // add glWindow to the UI
    m_ui->s_mainWindowGridLayout->addWidget(m_gl,0,0,4,1);

    populateUI();

    // now we wire up the UI components to the slots
    connect(m_ui->m_selected, SIGNAL(currentIndexChanged(int)), this, SLOT(selectRenderObject(int)));
    connect(m_ui->m_move, SIGNAL(valueChanged(double)), this, SLOT(setPosition(double)));
    connect(m_ui->m_sim, SIGNAL(clicked(bool)), this, SLOT(toggleSim(bool)));
    connect(m_ui->m_timerValue, SIGNAL(valueChanged(int)), m_gl, SLOT(setTimerDuration(int)));

//    toggleSim(false);
}

//----------------------------------------------------------------------------------------------------------------------
MainWindow::~MainWindow()
{
    delete m_scene;
    delete m_ui;
}

void MainWindow::setPosition(double _v)
{
    std::cout << "Move position" << std::endl;
    ElasticRod* strand = m_scene->getStrands()[0];
    strand->m_ppos[0][0] += _v;
}

//----------------------------------------------------------------------------------------------------------------------
void MainWindow::selectRenderObject(int index)
{
    index = index - 1;
    RenderObject* object = (index < 0)? NULL : m_scene->getRenderObjects()[index];
    m_gl->setSelectedObject(object);
}

//----------------------------------------------------------------------------------------------------------------------
void MainWindow::toggleSim(bool s)
{
    if(s == true)
	{
		m_gl->startSimTimer();
	}
	else
	{
		m_gl->stopSimTimer();
	}
}

void MainWindow::populateUI()
{
    m_ui->m_selected->addItem(QString("None"));
    int i = 1;
    typedef std::vector<RenderObject*>::const_iterator Iter;
    for (Iter it = m_scene->getRenderObjects().begin(); it != m_scene->getRenderObjects().end(); ++it)
    {
        QString label = QString("RenderObject_%1").arg(i);
        m_ui->m_selected->addItem(label);
        ++i;
    }
}
