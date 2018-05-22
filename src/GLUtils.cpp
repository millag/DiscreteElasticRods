#include "GLUtils.h"
#include "Mesh.h"

#include <QOpenGLContext>
#include <QOpenGLFunctions>

namespace
{
static inline GLenum ToGLPrimitive( Mesh::Primitive type )
{
	switch ( type )
	{
	case Mesh::TRIANGLES:
		return GL_TRIANGLES;
	case Mesh::LINES:
		return GL_LINE;
	case Mesh::POINTS:
		return GL_POINTS;
	}

	return GL_POINTS;
}
}

GLShaderManager* GLShaderManager::getInstance()
{
	static GLShaderManager shaderMan;
	return &shaderMan;
}

QOpenGLShaderProgram* GLShaderManager::loadShader(const std::string& shaderName,
                                                  const std::string& vertShaderPath,
                                                  const std::string& fragShaderPath,
                                                  bool reload)
{
	if (!reload)
	{
		auto result = getShader(shaderName);
		if (result)
		{
			return result;
		}
	}

	auto result = new QOpenGLShaderProgram();
	std::unique_ptr< QOpenGLShaderProgram > shaderProgram(result);
	if (!result->create())
	{
		qWarning() << "Unable to create shader program";
		return nullptr;
	}

	if (!result->addShaderFromSourceFile(QOpenGLShader::Vertex, vertShaderPath.c_str()))
	{
		qWarning() << "Vertex shader compilation failed";
		return nullptr;
	}

	if (!result->addShaderFromSourceFile(QOpenGLShader::Fragment, fragShaderPath.c_str()))
	{
		qWarning() << "Fragment shader compilation failed";
		return nullptr;
	}

	if (!result->link())
	{
		qWarning() << "Shader program not linked";
		return nullptr;
	}

	m_glShaderCache[ shaderName ] = std::move(shaderProgram);
	return result;
}

QOpenGLShaderProgram* GLShaderManager::loadShader(const std::string& shaderName,
                                                  const std::string& vertShaderPath,
                                                  const std::string& fragShaderPath,
                                                  const std::string& geomShaderPath,
                                                  const std::string& tessCtrlShaderPath,
                                                  const std::string& tessEvalShaderPath,
                                                  bool reload)
{
	if (!reload)
	{
		auto result = getShader(shaderName);
		if (result)
		{
			return result;
		}
	}

	auto result = new QOpenGLShaderProgram();
	std::unique_ptr< QOpenGLShaderProgram > shaderProgram(result);
	if (!result->create())
	{
		qWarning() << "Unable to create shader program";
		return nullptr;
	}

	if (!result->addShaderFromSourceFile(QOpenGLShader::Vertex, vertShaderPath.c_str()))
	{
		qWarning() << "Vertex shader compilation failed";
		return nullptr;
	}

	if (!result->addShaderFromSourceFile(QOpenGLShader::Fragment, fragShaderPath.c_str()))
	{
		qWarning() << "Fragment shader compilation failed";
		return nullptr;
	}

	if (   !geomShaderPath.empty()
	    && !result->addShaderFromSourceFile(QOpenGLShader::Geometry, geomShaderPath.c_str()))
	{
		qWarning() << "Geometry shader compilation failed";
		return nullptr;
	}

	if (   !tessCtrlShaderPath.empty()
	    && !result->addShaderFromSourceFile(QOpenGLShader::TessellationControl, tessCtrlShaderPath.c_str()))
	{
		qWarning() << "Tesselation control shader compilation failed";
		return nullptr;
	}

	if (   !tessEvalShaderPath.empty()
	    && !result->addShaderFromSourceFile(QOpenGLShader::TessellationEvaluation, tessEvalShaderPath.c_str()))
	{
		qWarning() << "Tesselation eval shader compilation failed";
		return nullptr;
	}

	if (!result->link())
	{
		qWarning() << "Shader program not linked";
		return nullptr;
	}

	m_glShaderCache[ shaderName ] = std::move(shaderProgram);
	return result;
}


GLDrawable::GLDrawable():
    m_vvbo( QOpenGLBuffer::VertexBuffer ),
    m_nvbo( QOpenGLBuffer::VertexBuffer ),
    m_ibo( QOpenGLBuffer::IndexBuffer )
{ }

GLDrawable::~GLDrawable()
{
	invalidate();
}

void GLDrawable::invalidate()
{
	if ( m_vao.isCreated() )
	{
		m_vao.destroy();
	}

	if ( m_vvbo.isCreated() )
	{
		m_vvbo.destroy();
	}

	if ( m_nvbo.isCreated() )
	{
		m_nvbo.destroy();
	}

	if ( m_ibo.isCreated() )
	{
		m_ibo.destroy();
	}

	m_shaderProgram = nullptr;
	m_nElements = 0;
}

void GLDrawable::draw()
{
	if ( !isValid() )
	{
		return;
	}

	auto gl = QOpenGLContext::currentContext()->functions();

	m_shaderProgram->bind();
	m_shaderProgram->setUniformValue( "color", m_color[0], m_color[1], m_color[2], 1.f );

	m_vao.bind();
	if ( m_ibo.isCreated() )
	{
		gl->glDrawElements( m_primitive, m_nElements, GL_UNSIGNED_INT, nullptr );
	}
	else
	{
		gl->glDrawArrays( m_primitive, 0, m_nElements );
	}

	m_vao.release();
}


bool GLDrawable::createGrid( unsigned usize, unsigned vsize, QOpenGLShaderProgram& shader, GLDrawable &o_out )
{
	o_out.invalidate();

	if ( !o_out.m_vao.create() )
	{
		qWarning() << "Unable to create grid VAO";
		return false;
	}

	if ( !o_out.m_vvbo.create() )
	{
		qWarning() << "Unable to create grid VBO";
		o_out.invalidate();
		return false;
	}

	// create grid geometry
	mg::Vec3D bleft(-1, 0, -1);
	mg::Vec3D bright(1, 0, -1);
	mg::Vec3D tleft(-1, 0, 1);
	mg::Vec3D tright(1, 0, 1);

	std::vector<mg::Vec3D> v( 2 * ( usize + vsize ) );
	for ( auto i = 0u; i < usize; ++i )
	{
		auto t = static_cast<float>( i ) / ( usize - 1 );
		v[2 * i] = mg::lerp( bleft, bright, t );
		v[2 * i + 1] = mg::lerp( tleft, tright, t );
	}

	const auto offset = 2 * usize;
	for( auto i = 0u; i < vsize; ++i )
	{
		auto t = static_cast<float>( i ) / ( vsize - 1 );
		v[offset + 2 * i] = mg::lerp( bleft, tleft, t );
		v[offset + 2 * i + 1] = mg::lerp( bright, tright, t );
	}

	o_out.m_vao.bind();

	o_out.m_shaderProgram = &shader;
	o_out.m_shaderProgram->bind();

	o_out.m_vvbo.bind();
	o_out.m_vvbo.setUsagePattern( QOpenGLBuffer::StaticDraw );
	o_out.m_vvbo.allocate( v.data(),
	                        static_cast<int>( v.size() ) * sizeof( v[0] ) );

	o_out.m_shaderProgram->enableAttributeArray( "position" );
	o_out.m_shaderProgram->setAttributeBuffer( "position",
	                                           GL_FLOAT,
	                                           0,
	                                           sizeof( v[0] ) / sizeof( v[0][0] ) );

	o_out.m_vao.release();

	o_out.m_nElements = static_cast<int>( v.size() );
	o_out.m_primitive = GL_LINES;
	o_out.m_color[0] = o_out.m_color[1] = o_out.m_color[2] = 0.8f;
	o_out.m_transform.identity();

	return true;
}

bool GLDrawable::createFrom( const Mesh& mesh, QOpenGLShaderProgram& shader, GLDrawable& o_out )
{
	o_out.invalidate();

	if ( !o_out.m_vao.create() )
	{
		qWarning() << "Unable to create VAO";
		return false;
	}

	if ( !o_out.m_vvbo.create() )
	{
		qWarning() << "Unable to create vertices VBO";
		o_out.invalidate();
		return false;
	}

	o_out.m_primitive = ToGLPrimitive( mesh.getPrimitiveType() );
	o_out.m_nElements = static_cast<int>( mesh.getNPrimitives() * mesh.getNVerticesPerPrimitive() );
	o_out.m_shaderProgram = &shader;
	o_out.m_color[0] = 0.9f, o_out.m_color[1] = 0.f, o_out.m_color[2] = 0.f;
	o_out.m_transform.identity();

	o_out.m_vao.bind();
	o_out.m_shaderProgram->bind();

	o_out.m_vvbo.bind();
	o_out.m_vvbo.setUsagePattern( QOpenGLBuffer::DynamicDraw );
	o_out.m_vvbo.allocate( mesh.m_vertices.data(),
	                       static_cast<int>( mesh.m_vertices.size() ) * sizeof( mesh.m_vertices[0] ) );

	o_out.m_shaderProgram->enableAttributeArray( "position" );
	o_out.m_shaderProgram->setAttributeBuffer( "position",
	                                           GL_FLOAT,
	                                           0,
	                                           sizeof( mesh.m_vertices[0] ) / sizeof( mesh.m_vertices[0][0] ) );

	if ( mesh.hasNormals() && o_out.m_nvbo.create() )
	{
		o_out.m_nvbo.bind();
		o_out.m_nvbo.setUsagePattern( QOpenGLBuffer::DynamicDraw );
		o_out.m_nvbo.allocate( mesh.m_normals.data(),
		                       static_cast<int>( mesh.m_normals.size() ) * sizeof( mesh.m_normals[0] ) );

		o_out.m_shaderProgram->enableAttributeArray( "normal" );
		o_out.m_shaderProgram->setAttributeBuffer( "normal",
		                                           GL_FLOAT,
		                                           0,
		                                           sizeof( mesh.m_normals[0] ) / sizeof( mesh.m_normals[0][0] ) );
	}

	if ( mesh.m_vindices.size() && o_out.m_ibo.create() )
	{
		o_out.m_ibo.bind();
		o_out.m_ibo.setUsagePattern( QOpenGLBuffer::DynamicDraw );
		o_out.m_ibo.allocate( mesh.m_vindices.data(),
		                       static_cast<int>( mesh.m_vindices.size() ) * sizeof( mesh.m_vindices[0] ) );
	}

	o_out.m_vao.release();
	return true;
}

GLLight::GLLight():
    m_ambient( 0.f, 0.f, 0.f ),
    m_diffuse( 1.f, 1.f, 1.f ),
    m_specular( 1.f, 1.f, 1.f ),
    m_position( 0.f, 0.f, 0.f )
{}

void GLLight::loadToShader( QOpenGLShaderProgram& shader ) const
{
	shader.setUniformValue( "light.ambient", m_ambient[0], m_ambient[1], m_ambient[2] );
	shader.setUniformValue( "light.diffuse", m_diffuse[0], m_diffuse[1], m_diffuse[2] );
	shader.setUniformValue( "light.specular", m_specular[0], m_specular[1], m_specular[2] );
	shader.setUniformValue( "light.position", m_position[0], m_position[1], m_position[2] );
}

GLMaterial::GLMaterial():
    m_ambient( 0.f, 0.f, 0.f ),
    m_diffuse( 0.7f, 0.f, 0.f ),
    m_specular( 1.f, 1.f, 1.f ),
    m_shininess( 16.f )
{}

void GLMaterial::loadToShader( QOpenGLShaderProgram& shader ) const
{
	shader.setUniformValue( "mtl.ambient", m_ambient[0], m_ambient[1], m_ambient[2] );
	shader.setUniformValue( "mtl.diffuse", m_diffuse[0], m_diffuse[1], m_diffuse[2] );
	shader.setUniformValue( "mtl.specular", m_specular[0], m_specular[1], m_specular[2] );
	shader.setUniformValue( "mtl.shininess", m_shininess );
}
