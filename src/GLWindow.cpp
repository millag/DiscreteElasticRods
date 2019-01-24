#include "GLWindow.h"
#include "Utils.h"

#include <QMouseEvent>
#include <QOpenGLContext>
#include <QOpenGLFunctions>
#include <QOpenGLExtraFunctions>
#include <QOpenGLShaderProgram>

namespace
{

static void buildVAOs( const std::vector<Mesh*>& meshList, QOpenGLShaderProgram& shader, DrawList& o_drawList )
{
	auto cnt = o_drawList.size();
	o_drawList.resize( cnt + meshList.size() );

	for ( auto it = meshList.begin(); it != meshList.end(); ++it, ++cnt )
	{
		if ( ( *it ) && ( *it )->getNPrimitives() <= 0 )
		{
			continue;
		}

		auto drawable = std::make_unique<GLDrawable>();
		if ( GLDrawable::createFrom( **it, shader, *drawable ) )
		{
			o_drawList[cnt] = std::move( drawable );
		}
	}
}

}

GLViewport::GLViewport( QWidget* parent , Qt::WindowFlags format ):
    QOpenGLWidget( parent, format ),
    m_transformHdl( &m_cam )
{
	// re-size the widget to that of the parent
	if ( parent )
	{
		resize( parent->size() );
	}

	// set this widget to have the initial keyboard focus
	setFocus();
}

GLViewport::~GLViewport()
{
	qInfo() << "Shutting down viewport";
	m_drawList.clear();
}

void GLViewport::initializeGL()
{
//	this virtual function is called once before the first call to paintGL() or resizeGL(),
//	and then once whenever the widget has been assigned a new QOpenGLContext.
//	this function should set up any required OpenGL context rendering flags, defining display lists, etc.

	if ( !m_renderer.initialize() )
	{
		qWarning() << "Faled to initialize OpenGL renderer";
		return;
	}

	m_cam.lookAt( mg::Vec3D( 1.f, 2.f, 3.f ),
	              mg::Vec3D( 0.f, 0.f, 0.f ),
	              mg::Vec3D( 0.f, 1.f, 0.f ) );

	auto shaderMan = GLShaderManager::getInstance();
	auto shader = shaderMan->getShader( "Constant" );

//	reference grid
	GLDrawable::createGrid( 10, 10, *shader, m_refGrid );

//	scene
	if ( m_scene )
	{
		auto shader = shaderMan->getShader( "Phong" );
		buildVAOs( m_scene->getMeshes(), *shader, m_drawList );
	}

//	m_text = new ngl::Text(QFont("Arial",13));
//	m_text->setColour(0.7,0.7,0.7);
//	m_text->setScreenSize(width(), height());

//	// now create our light this is done after the camera so we can pass the
//	// transpose of the projection matrix to the light to do correct eye space
//	// transformations
//	ngl::Light light(ngl::Vec3(0,0,0),ngl::Colour(1,1,1,1),ngl::Colour(1,1,1,1),ngl::POINTLIGHT);
//	ngl::Mat4 iv=m_cam->getViewMatrix();
//	iv.transpose();
//	light.setTransform(iv);
//	light.setAttenuation(1,0,0);
//	light.enable();
//	light.loadToShader("light");

//	// the shader will use the currently active material and light0 so set them
//	ngl::Material m(ngl::SILVER);
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

void GLViewport::resizeGL( int w, int h )
{
	const auto aspect = static_cast<mg::Real>( w ) / h;
	m_cam.perspective( mg::Constants::pi() / 3, aspect, 0.001f, 1000.f );
}

void GLViewport::paintGL()
{
	auto gl = QOpenGLContext::currentContext()->functions();
//	clear viewport
	gl->glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	m_renderer.setCamera( m_cam );

//	draw reference frame
	m_renderer.beginDrawable();
	m_renderer.setColor( mg::Vec3D( 1.f, 0.f, 0.f ) );
	m_renderer.line( mg::Vec3D( 0.f, 0.f, 0.f ), mg::Vec3D( 1.f, 0.f, 0.f ) );
	m_renderer.setColor( mg::Vec3D( 0.f, 1.f, 0.f ) );
	m_renderer.line( mg::Vec3D( 0.f, 0.f, 0.f ), mg::Vec3D( 0.f, 1.f, 0.f ) );
	m_renderer.setColor( mg::Vec3D( 0.f, 0.f, 1.f ) );
	m_renderer.line( mg::Vec3D( 0.f, 0.f, 0.f ), mg::Vec3D( 0.f, 0.f, 1.f ) );
	m_renderer.endDrawable();

//	draw reference grid
	if ( m_refGrid.isValid() )
	{
		auto shader = m_refGrid.getShader();
		shader->bind();

		const mg::Matrix4D mvp = m_renderer.getViewProjectionMatrix() * m_refGrid.getTransform();
		shader->setUniformValue( "mvp", *reinterpret_cast<const GLMatrix4x4*>( mvp.data() ) );
		m_refGrid.draw();
	}

//	m_renderer.beginDrawable();
//	m_renderer.setColor(mg::Vec3D(0.3f, 0.5f, 0.1f));
//	m_renderer.box(mg::Vec3D(-2.f, 0.f, 0.f), mg::Vec3D(0.f, 1.f, 1.f),  mg::Vec3D(1.f, 1.f, 1.f));
//	m_renderer.setColor(mg::Vec3D(0.1f, 0.5f, 0.5f));
//	m_renderer.sphere(mg::Vec3D(2.f, 0.f, 0.f), 0.5f);
//	m_renderer.setColor(mg::Vec3D(0.3f, 0.1f, 0.7f));
//	m_renderer.cone(mg::Vec3D(0.f, 0.f, 0.f), mg::Vec3D(1.f, 0.f, 0.f), 0.9f, 0.2f);
//	m_renderer.endDrawable();

	if( !m_scene )
	{
		return;
	}

	auto roList = m_scene->getRenderObjects();
	for ( auto it = roList.begin(); it != roList.end(); ++it )
	{
		auto ro = ( *it );
		if ( !ro )
		{
			continue;
		}

		auto& drawable =  m_drawList.at( ro->getMeshId() );
		if ( !drawable || !drawable->isValid() )
		{
			continue;
		}

		auto shader = drawable->getShader();
		shader->bind();

		const mg::Matrix4D mv = m_renderer.getViewMatrix() * ro->getTransform();
		const mg::Matrix4D mvp = m_renderer.getViewProjectionMatrix() * ro->getTransform();
		shader->setUniformValue( "mv", *reinterpret_cast<const GLMatrix4x4*>( mv.data() ) );
		shader->setUniformValue( "mvp", *reinterpret_cast<const GLMatrix4x4*>( mvp.data() ) );
		drawable->draw();
	}

	auto hair = m_scene->getHairById( 0 );
	if ( !hair )
	{
		return;
	}

	auto shaderMan = GLShaderManager::getInstance();

#ifdef DBUGG
	auto shader = shaderMan->getShader( "DebugRod" );
	GLDrawable::createFrom( *hair, *shader, m_hair );

	shader->bind();

	gl->glLineWidth( 0.05f );
	const mg::Matrix4D mv = m_renderer.getViewMatrix() * m_hair.getTransform();
	const mg::Matrix4D mvp = m_renderer.getViewProjectionMatrix() * m_hair.getTransform();
	shader->setUniformValue( "mv", *reinterpret_cast<const GLMatrix4x4*>( mv.data() ) );
	shader->setUniformValue( "mvp", *reinterpret_cast<const GLMatrix4x4*>( mvp.data() ) );

	m_hair.draw();
#else
	auto shader = shaderMan->getShader( "Tube" );
	GLDrawable::createFrom( *hair, *shader, m_hair );

	shader->bind();

	const mg::Matrix4D mv = m_renderer.getViewMatrix() * m_hair.getTransform();
	const mg::Matrix4D mvp = m_renderer.getViewProjectionMatrix() * m_hair.getTransform();
	shader->setUniformValue( "mv", *reinterpret_cast<const GLMatrix4x4*>( mv.data() ) );
	shader->setUniformValue( "mvp", *reinterpret_cast<const GLMatrix4x4*>( mvp.data() ) );
	shader->setUniformValue( "radius", 0.1f );

	m_hair.draw();
#endif

//	m_text->renderText(10, 10, QString("TIME: %1ms").arg(chronometer.elapsed()));
//	m_text->renderText(10, 30, QString("FPS: %1").arg((float)mg::SEC / chronometer.restart()));
//	m_text->renderText(10, 50, QString("strands: %1").arg(m_scene->getHairById(0)->m_strands.size()));
//	m_text->renderText(10, 70, QString("points per strand: %1").arg(m_scene->getHairById(0)->m_params->m_nParticles));
}

void GLViewport::mousePressEvent( QMouseEvent* event )
{
	TransformHandle::TransformMode mode = TransformHandle::TM_None;
	switch ( event->button() )
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
	default:
		break;
	}

	m_mouseX = event->x();
	m_mouseY = event->y();
	m_transformHdl.setMode( mode );
}

void GLViewport::mouseMoveEvent( QMouseEvent* event )
{
	if ( m_transformHdl.isActive() )
	{
		const auto dx = static_cast<mg::Real>( event->x() - m_mouseX ) / width();
		const auto dy = static_cast<mg::Real>( event->y() - m_mouseY ) / height();
		m_mouseX = event->x();
		m_mouseY = event->y();

		m_transformHdl.update( dx, dy );
		update();
	}
}

void GLViewport::mouseReleaseEvent( QMouseEvent* event )
{
	UNUSED_VALUE( event );

	m_transformHdl.setMode( TransformHandle::TM_None );
}

#ifdef DBUGG
void GLWindow::drawHairStrand(const ElasticRod& strand)
{
	if (m_strandVAO != nullptr)
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
