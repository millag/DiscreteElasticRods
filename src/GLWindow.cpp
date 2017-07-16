#include "GLWindow.h"
#include "Utils.h"

#include <QMouseEvent>
#include <QOpenGLContext>
#include <QOpenGLFunctions>
#include <QOpenGLExtraFunctions>
#include <QOpenGLShaderProgram>


GLViewport::GLViewport(QWidget* parent , Qt::WindowFlags f):
    QOpenGLWidget(parent, f),
    m_transformHdl(&m_cam)
{
    // re-size the widget to that of the parent (in this case the GLFrame passed in on construction)
    if (parent)
    {
        resize(parent->size());
    }

    // set this widget to have the initial keyboard focus
    setFocus();
}

GLViewport::~GLViewport()
{
    qInfo() << "Shutting down GL viewport";
    m_drawList.clear();
}

void GLViewport::initializeGL()
{
//	This virtual function is called once before the first call to paintGL() or resizeGL(),
//	and then once whenever the widget has been assigned a new QOpenGLContext.
//	This function should set up any required OpenGL context rendering flags, defining display lists, etc.

    auto gl = QOpenGLContext::currentContext()->functions();

    gl->glClearColor(0.2f, 0.2f, 0.2f, 1.0f);
    gl->glEnable(GL_DEPTH_TEST);
    gl->glEnable(GL_LINE_SMOOTH);

    auto shaderMan = GLShaderManager::getInstance();
    auto shader = shaderMan->loadShader("Color",
                                        "shaders/ColorVert.glsl",
                                        "shaders/ColorFrag.glsl");

    shaderMan->loadShader("Phong",
                          "shaders/PhongVertex.glsl",
                          "shaders/PhongFragment.glsl");

    shaderMan->loadShader("Tube",
                          "shaders/TubeVert.glsl",
                          "shaders/TubeFrag.glsl",
                          "",
                          "shaders/TubeTCS.glsl",
                          "shaders/TubeTES.glsl"
                          );

    shaderMan->loadShader("DebugRod",
                          "shaders/DebugVert.glsl",
                          "shaders/DebugFrag.glsl",
                          "shaders/DebugGeom.glsl",
                          "",
                          ""
                          );

    m_cam.lookAt( mg::Vec3D(1.f, 2.f, 3.f), mg::Vec3D(0.f, 0.f, 0.f), mg::Vec3D(0.f, 1.f, 0.f));
    GLDrawable::createGrid(m_refGrid, 10, 10, *shader);

    if (m_scene)
    {
        buildVAOs(m_scene->getMeshes(), m_drawList);
    }

//	m_text = new ngl::Text(QFont("Arial",13));
//	m_text->setColour(0.7,0.7,0.7);
//	m_text->setScreenSize(width(), height());

//	// now bind the shader attributes for most NGL primitives we use the following
//	// layout attribute 0 is the vertex data (x,y,z)
//	shader->bindAttribute("Phong",0,"inVert");
//	shader->bindAttribute("Phong",1,"inUV");
//	shader->bindAttribute("Phong",2,"inNormal");
//	// now we have associated this data we can link the shader
//	shader->linkProgramObject("Phong");
//	// and make it active ready to load values
//	(*shader)["Phong"]->use();
//	shader->setShaderParam1i("Normalize",1);

//	// now pass the modelView and projection values to the shader
//	// the shader will use the currently active material and light0 so set them
//	ngl::Material m(ngl::SILVER);
//	m.loadToShader("material");
//	ngl::Light light(ngl::Vec3(0,0,0),ngl::Colour(1,1,1,1),ngl::Colour(1,1,1,1),ngl::POINTLIGHT);

//	// now create our light this is done after the camera so we can pass the
//	// transpose of the projection matrix to the light to do correct eye space
//	// transformations
//	ngl::Mat4 iv=m_cam->getViewMatrix();
//	iv.transpose();
//	light.setTransform(iv);
//	light.setAttenuation(1,0,0);
//	light.enable();
//	// load these values to the shader as well
//	light.loadToShader("light");

//	shader->bindAttribute("Tube",0,"inVert");
//	shader->bindAttribute("Tube",1,"inNormal");

//	shader->linkProgramObject("Tube");
//	(*shader)["Tube"]->use();
//	light.loadToShader("light");
//	m.loadToShader("material");

//	shader->bindAttribute("DebugRod",0,"inVert");
//	shader->bindAttribute("DebugRod",1,"inKB");
//	shader->bindAttribute("DebugRod",2,"inM1");
//	shader->bindAttribute("DebugRod",3,"inM2");
//	shader->bindAttribute("DebugRod",4,"inForce");

//	shader->linkProgramObject("DebugRod");
//	(*shader)["DebugRod"]->use();
//	light.loadToShader("light");
//	m.loadToShader("material");
}

void GLViewport::resizeGL(int w, int h)
{
//	This virtual function is called whenever the widget has been resized.
//	The new size is passed in width and height.
    const auto aspect = (mg::Real)(w) / h;
    m_cam.perspective( mg::Constants::pi() / 3, aspect, 0.001f, 1000.f );
}

//void GLViewport::loadMatricesToShader()
//{
//	ngl::ShaderLib *shader=ngl::ShaderLib::instance();

//	ngl::Mat4 MV;
//	ngl::Mat4 MVP;
//	ngl::Mat3 normalMatrix;
//	ngl::Mat4 M;
//	M = m_transform.getMatrix();
//	MV = M * m_cameraTransform.getMatrix() * m_cam->getViewMatrix();
//	MVP = MV * m_cam->getProjectionMatrix();
//	normalMatrix = MV;
//	normalMatrix.inverse();
//	shader->setShaderParamFromMat4("MV",MV);
//	shader->setShaderParamFromMat4("MVP",MVP);
//	shader->setShaderParamFromMat3("normalMatrix",normalMatrix);
//	shader->setShaderParamFromMat4("M",M);
//}

//void GLViewport::loadMatricesToHairShader()
//{
//	ngl::ShaderLib *shader=ngl::ShaderLib::instance();

//	ngl::Mat4 M;
//	ngl::Mat4 MV;
//	ngl::Mat4 MVP;
//	M = m_transform.getMatrix();
//	MV = M  * m_cameraTransform.getMatrix() * m_cam->getViewMatrix();
//	MVP = MV * m_cam->getProjectionMatrix();
//	shader->setShaderParamFromMat4("MVP",MVP);
//	shader->setShaderParamFromMat4("M",M);
//}

//void GLViewport::loadMatricesToHairShader2()
//{
//	ngl::ShaderLib *shader=ngl::ShaderLib::instance();

//	ngl::Mat4 MV;
//	ngl::Mat4 MVP;
//	ngl::Mat3 normalMatrix;
//	ngl::Mat4 M;
//	M = m_transform.getMatrix();
//	MV = M * m_cameraTransform.getMatrix() * m_cam->getViewMatrix();
//	MVP = MV * m_cam->getProjectionMatrix();
//	normalMatrix = MV;
//	normalMatrix.inverse();
//	shader->setShaderParamFromMat4("MV",MV);
//	shader->setShaderParamFromMat4("MVP",MVP);
//	shader->setShaderParamFromMat3("normalMatrix",normalMatrix);
//}

void GLViewport::paintGL()
{
//	This virtual function is called whenever the widget needs to be painted.
//	this is our main drawing routine

    auto gl = QOpenGLContext::currentContext()->functions();
//	clear the screen and depth buffer
    gl->glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    auto vp = m_cam.getVPMatrix();
//	draw reference grid
    if (m_refGrid.isValid())
    {
        auto shader = m_refGrid.getShader();
        shader->bind();

        shader->setUniformValue("mvp", QMatrix4x4((vp * m_refGrid.getTransform()).data()));
        m_refGrid.draw();
    }

    if (m_scene)
    {
        auto roList = m_scene->getRenderObjects();
        for (auto it = roList.begin(); it != roList.end(); ++it)
        {
            RenderObject* ro = (*it);
            if (!ro)
            {
                continue;
            }

            auto& drawable =  m_drawList.at(ro->getMeshId());
            if (!drawable || !drawable->isValid())
            {
                continue;
            }

            auto shader = drawable->getShader();
            shader->bind();

            shader->setUniformValue( "mvp", QMatrix4x4((vp * drawable->getTransform()).data()));
            drawable->draw();
        }
    }

//	if (m_scene->getHairById(0) == NULL)
//	{
//		return;
//	}

//#ifndef DBUGG
//	shader->use("Tube");
//	shader->setShaderParam4f("Col1", 0.47, 0.4, 0.0, 1.0);
//	shader->setShaderParam4f("Col2", 0.0, 0.0, 0.0, 1.0);
//	shader->setShaderParam1f("Radius", 0.1);
//	glPatchParameteri(GL_PATCH_VERTICES, 4);
//	// load transform stack
//	m_transform.reset();
//	loadMatricesToHairShader2();
//#endif

//#ifdef DBUGG
//	glLineWidth(0.05);
//	shader->use("DebugRod");
//	shader->setShaderParam4f("Colour", 0.8, 0.8, 0.0, 1.0);
//	// load transform stack
//	m_transform.reset();
//	loadMatricesToHairShader();
//#endif

//	const std::vector<ElasticRod*>& strands = m_scene->getHairById(0)->m_strands;
////    const std::vector<ElasticRod*>& strands = m_scene->getStrands();
//	typedef std::vector<ElasticRod*>::const_iterator SIter;
//	for (SIter it = strands.begin(); it != strands.end(); ++it)
//	{
//		drawHairStrand(**it);
//	}

//	m_text->renderText(10, 10, QString("TIME: %1ms").arg(chronometer.elapsed()));
//	m_text->renderText(10, 30, QString("FPS: %1").arg((float)mg::SEC / chronometer.restart()));
//	m_text->renderText(10, 50, QString("strands: %1").arg(m_scene->getHairById(0)->m_strands.size()));
//	m_text->renderText(10, 70, QString("points per strand: %1").arg(m_scene->getHairById(0)->m_params->m_nParticles));
}


// mouse controls ----------------------------------
void GLViewport::mousePressEvent(QMouseEvent* event)
{
    m_mouseX = event->x();
    m_mouseY = event->y();

    TransformHandle::TransformMode mode = TransformHandle::TM_None;
    switch (event->button())
    {
        case Qt::LeftButton:
        {
            mode = TransformHandle::TM_Rotation;
            break;
        }
        case Qt::MidButton:
        {
            mode = TransformHandle::TM_Translation;
            break;
        }
        case Qt::RightButton:
        {
            mode = TransformHandle::TM_Scale;
            break;
        }
    }

    m_transformHdl.setMode(mode);
}

void GLViewport::mouseMoveEvent(QMouseEvent* event)
{
    if (m_transformHdl.isActive())
    {
        const mg::Real dx = (mg::Real) (event->x() - m_mouseX) / width();
        const mg::Real dy = (mg::Real) (event->y() - m_mouseY) / height();

        m_transformHdl.update(dx, dy);

        m_mouseX = event->x();
        m_mouseY = event->y();
        update();
    }
}

void GLViewport::mouseReleaseEvent(QMouseEvent* event)
{
    m_transformHdl.setMode(TransformHandle::TM_None);
}


// utility functions -------------------------------
void GLViewport::buildVAOs(const std::vector< Mesh* >& meshList, DrawList& o_drawList) const
{
    auto shaderMan = GLShaderManager::getInstance();
    auto shader = shaderMan->getShader("Color");
    if (!shader)
    {
        return;
    }

    o_drawList.resize(o_drawList.size() + meshList.size());

    std::size_t i = 0;
    for (auto it = meshList.begin(); it != meshList.end(); ++it, ++i)
    {
        if ((*it) && (*it)->getNPrimitives() <= 0)
        {
            continue;
        }

        auto drawablePtr = new GLDrawable();
        std::unique_ptr< GLDrawable > drawable(drawablePtr);
        if (GLDrawable::createFrom(*drawable, **it, *shader))
        {
            o_drawList[i] = std::move(drawable);
        }
    }
}


#ifdef DBUGG
void GLWindow::drawHairStrand(const ElasticRod& strand)
{
    if (m_strandVAO != NULL)
    {
        m_strandVAO->bind();

        m_strandVAO->updateIndexedData(0, strand.m_ppos.size() * sizeof(mg::Vec3D),
                                       strand.m_ppos[0][0]);
        m_strandVAO->setVertexAttributePointer(0, 3, GL_FLOAT, 0, 0);

        m_strandVAO->updateIndexedData(1, strand.m_kb.size() * sizeof(mg::Vec3D),
                                       strand.m_kb[0][0]);
        m_strandVAO->setVertexAttributePointer(1, 3, GL_FLOAT, 0, 0);

        m_strandVAO->updateIndexedData(2, strand.m_m1.size() * sizeof(mg::Vec3D),
                                       strand.m_m1[0][0]);
        m_strandVAO->setVertexAttributePointer(2, 3, GL_FLOAT, 0, 0);

        m_strandVAO->updateIndexedData(3, strand.m_m2.size() * sizeof(mg::Vec3D),
                                       strand.m_m2[0][0]);
        m_strandVAO->setVertexAttributePointer(3, 3, GL_FLOAT, 0, 0);

        m_strandVAO->updateIndexedData(4, strand.m_elasticForce.size() * sizeof(mg::Vec3D),
                                       strand.m_elasticForce[0][0]);
        m_strandVAO->setVertexAttributePointer(4, 3, GL_FLOAT, 0, 0);

        m_strandVAO->draw();
        m_strandVAO->unbind();
        return;
    }

    std::vector<unsigned> indices(strand.m_ppos.size() + 2);
    unsigned i = 0;
    indices[i] = 0;
    for (i = 0; i < strand.m_ppos.size(); ++i)
    {
        indices[i + 1] = i;
    }
    indices[i + 1] = i - 1;

    m_strandVAO = ngl::VertexArrayObject::createVOA(GL_LINE_STRIP_ADJACENCY);
    m_strandVAO->bind();

    m_strandVAO->setIndexedData(strand.m_ppos.size() * sizeof(mg::Vec3D),
                                strand.m_ppos[0][0],
                                indices.size(),
                                &indices[0],
                                GL_UNSIGNED_INT);
    m_strandVAO->setVertexAttributePointer(0, 3, GL_FLOAT, 0, 0);

    m_strandVAO->setIndexedData(strand.m_kb.size() * sizeof(mg::Vec3D),
                                strand.m_kb[0][0],
                                indices.size(),
                                &indices[0],
                                GL_UNSIGNED_INT);
    m_strandVAO->setVertexAttributePointer(1, 3, GL_FLOAT, 0, 0);

    m_strandVAO->setIndexedData(strand.m_m1.size() * sizeof(mg::Vec3D),
                                strand.m_m1[0][0],
                                indices.size(),
                                &indices[0],
                                GL_UNSIGNED_INT);
    m_strandVAO->setVertexAttributePointer(2, 3, GL_FLOAT, 0, 0);

    m_strandVAO->setIndexedData(strand.m_m2.size() * sizeof(mg::Vec3D),
                                strand.m_m2[0][0],
                                indices.size(),
                                &indices[0],
                                GL_UNSIGNED_INT);
    m_strandVAO->setVertexAttributePointer(3, 3, GL_FLOAT, 0, 0);

    m_strandVAO->setIndexedData(strand.m_elasticForce.size() * sizeof(mg::Vec3D),
                                strand.m_elasticForce[0][0],
                                indices.size(),
                                &indices[0],
                                GL_UNSIGNED_INT);
    m_strandVAO->setVertexAttributePointer(4, 3, GL_FLOAT, 0, 0);

    m_strandVAO->setNumIndices(indices.size());

    m_strandVAO->draw();
    m_strandVAO->unbind();
}
#endif

#ifndef DBUGG
void GLViewport::drawHairStrand(const ElasticRod& strand)
{
//	if (m_strandVAO != NULL)
//	{
//		m_strandVAO->bind();

//		m_strandVAO->updateIndexedData(0, strand.m_ppos.size() * sizeof(mg::Vec3D),
//									   strand.m_ppos[0][0]);
//		m_strandVAO->setVertexAttributePointer(0, 3, GL_FLOAT, 0, 0);

//		m_strandVAO->updateIndexedData(1, strand.m_m1.size() * sizeof(mg::Vec3D),
//									   strand.m_m1[0][0]);
//		m_strandVAO->setVertexAttributePointer(1, 3, GL_FLOAT, 0, 0);

//		m_strandVAO->draw();
//		m_strandVAO->unbind();
//		return;
//	}

//	std::vector<unsigned> indices( (strand.m_ppos.size() - 1) * 4 );
//	int nPoints = (int)strand.m_ppos.size();
//	for (int i = -1; i < nPoints - 2; ++i)
//	{
//		for (int j = 0; j < 4; ++j)
//		{
//			assert( (i + 1) >= 0 && (4*(i + 1) + j) < (int)indices.size());
//			indices[4*(i + 1) + j] = std::min(std::max(0, i + j), nPoints - 1);
//		}
//	}

//	m_strandVAO = ngl::VertexArrayObject::createVOA(GL_PATCHES);
//	m_strandVAO->bind();

//	m_strandVAO->setIndexedData(strand.m_ppos.size() * sizeof(mg::Vec3D),
//								strand.m_ppos[0][0],
//								indices.size(),
//								&indices[0],
//								GL_UNSIGNED_INT);
//	m_strandVAO->setVertexAttributePointer(0, 3, GL_FLOAT, 0, 0);

//	m_strandVAO->setIndexedData(strand.m_m1.size() * sizeof(mg::Vec3D),
//								strand.m_m1[0][0],
//								indices.size(),
//								&indices[0],
//								GL_UNSIGNED_INT);
//	m_strandVAO->setVertexAttributePointer(1, 3, GL_FLOAT, 0, 0);

//	m_strandVAO->setNumIndices(indices.size());

//	m_strandVAO->draw();
//	m_strandVAO->unbind();
}
#endif
