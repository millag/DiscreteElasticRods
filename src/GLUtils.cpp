#include "GLUtils.h"
#include "Mesh.h"
#include "Hair.h"
#include <QOpenGLContext>
#include <QOpenGLFunctions>

namespace
{
static inline GLenum ToGLPrimitive( Mesh::Primitive type )
{
	switch ( type )
	{
	case Mesh::Primitive::TRIANGLES:
		return GL_TRIANGLES;
	case Mesh::Primitive::LINES:
		return GL_LINE;
	case Mesh::Primitive::POINTS:
		return GL_POINTS;
	default:
		break;
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

	auto shaderProg = std::make_unique<QOpenGLShaderProgram>();
	if (!shaderProg->create())
	{
		qWarning() << "Unable to create shader program";
		return nullptr;
	}

	if (!shaderProg->addShaderFromSourceFile(QOpenGLShader::Vertex, vertShaderPath.c_str()))
	{
		qWarning() << "Vertex shader compilation failed";
		return nullptr;
	}

	if (!shaderProg->addShaderFromSourceFile(QOpenGLShader::Fragment, fragShaderPath.c_str()))
	{
		qWarning() << "Fragment shader compilation failed";
		return nullptr;
	}

	if (!shaderProg->link())
	{
		qWarning() << "Shader program not linked";
		return nullptr;
	}

	auto result = shaderProg.get();
	m_glShaderCache[ shaderName ] = std::move( shaderProg );

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

	auto shaderProg = std::make_unique<QOpenGLShaderProgram>();
	if (!shaderProg->create())
	{
		qWarning() << "Unable to create shader program";
		return nullptr;
	}

	if (!shaderProg->addShaderFromSourceFile(QOpenGLShader::Vertex, vertShaderPath.c_str()))
	{
		qWarning() << "Vertex shader compilation failed";
		return nullptr;
	}

	if (!shaderProg->addShaderFromSourceFile(QOpenGLShader::Fragment, fragShaderPath.c_str()))
	{
		qWarning() << "Fragment shader compilation failed";
		return nullptr;
	}

	if (   !geomShaderPath.empty()
	    && !shaderProg->addShaderFromSourceFile(QOpenGLShader::Geometry, geomShaderPath.c_str()))
	{
		qWarning() << "Geometry shader compilation failed";
		return nullptr;
	}

	if (   !tessCtrlShaderPath.empty()
	    && !shaderProg->addShaderFromSourceFile(QOpenGLShader::TessellationControl, tessCtrlShaderPath.c_str()))
	{
		qWarning() << "Tesselation control shader compilation failed";
		return nullptr;
	}

	if (   !tessEvalShaderPath.empty()
	    && !shaderProg->addShaderFromSourceFile(QOpenGLShader::TessellationEvaluation, tessEvalShaderPath.c_str()))
	{
		qWarning() << "Tesselation eval shader compilation failed";
		return nullptr;
	}

	if (!shaderProg->link())
	{
		qWarning() << "Shader program not linked";
		return nullptr;
	}

	auto result = shaderProg.get();
	m_glShaderCache[ shaderName ] = std::move(shaderProg);

	return result;
}

GLLight::GLLight():
    m_ambient( 0.f, 0.f, 0.f ),
    m_diffuse( 1.f, 1.f, 1.f ),
    m_specular( 1.f, 1.f, 1.f ),
    m_position( 0.f, 0.f, 0.f )
{ }

void GLLight::loadToShader( QOpenGLShaderProgram& shader ) const
{
	shader.setUniformValue( "light.ambient", m_ambient[0], m_ambient[1], m_ambient[2] );
	shader.setUniformValue( "light.diffuse", m_diffuse[0], m_diffuse[1], m_diffuse[2] );
	shader.setUniformValue( "light.specular", m_specular[0], m_specular[1], m_specular[2] );
	shader.setUniformValue( "light.position", m_position[0], m_position[1], m_position[2] );
}

GLMaterial::GLMaterial():
    m_ambient( 0.f, 0.f, 0.f ),
    m_diffuse( 0.7f, 0.7f, 0.7f ),
    m_specular( 1.f, 1.f, 1.f ),
    m_shininess( 16.f )
{ }

void GLMaterial::loadToShader( QOpenGLShaderProgram& shader ) const
{
	shader.setUniformValue( "mtl.ambient", m_ambient[0], m_ambient[1], m_ambient[2] );
	shader.setUniformValue( "mtl.diffuse", m_diffuse[0], m_diffuse[1], m_diffuse[2] );
	shader.setUniformValue( "mtl.specular", m_specular[0], m_specular[1], m_specular[2] );
	shader.setUniformValue( "mtl.shininess", m_shininess );
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
	m_shaderProgram->setUniformValue( "color", m_color[0], m_color[1], m_color[2] );

	QOpenGLVertexArrayObject::Binder raiiVAO( &m_vao );

	if ( m_ibo.isCreated() )
	{
		gl->glDrawElements( m_primitive, m_nElements, GL_UNSIGNED_INT, nullptr );
	}
	else
	{
		gl->glDrawArrays( m_primitive, 0, m_nElements );
	}
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
	mg::Vec3D bleft( -1.f, 0.f, -1.f );
	mg::Vec3D bright( 1.f, 0.f, -1.f );
	mg::Vec3D tleft( -1.f, 0.f, 1.f );
	mg::Vec3D tright( 1.f, 0.f, 1.f );

	std::vector<mg::Vec3D> v( 2 * ( usize + vsize ) );
	for ( auto i = 0u; i < usize; ++i )
	{
		auto t = static_cast<mg::Real>( i ) / ( usize - 1 );
		v[2 * i] = mg::lerp( bleft, bright, t );
		v[2 * i + 1] = mg::lerp( tleft, tright, t );
	}

	const auto offset = 2 * usize;
	for( auto i = 0u; i < vsize; ++i )
	{
		auto t = static_cast<mg::Real>( i ) / ( vsize - 1 );
		v[offset + 2 * i] = mg::lerp( bleft, tleft, t );
		v[offset + 2 * i + 1] = mg::lerp( bright, tright, t );
	}

	o_out.m_primitive = GL_LINES;
	o_out.m_nElements = static_cast<int>( v.size() );
	o_out.m_shaderProgram = &shader;
	o_out.m_color[0] = o_out.m_color[1] = o_out.m_color[2] = 0.8f;
	o_out.m_transform.identity();

	QOpenGLVertexArrayObject::Binder raiiVAO( &o_out.m_vao );

	o_out.m_vvbo.bind();
	o_out.m_vvbo.setUsagePattern( QOpenGLBuffer::StaticDraw );
	o_out.m_vvbo.allocate( v.data(),
	                        static_cast<int>( v.size() * sizeof( v[0] ) ) );

	o_out.m_shaderProgram->bind();
	o_out.m_shaderProgram->enableAttributeArray( "position" );
	o_out.m_shaderProgram->setAttributeBuffer( "position",
	                                           GL_FLOAT,
	                                           0,
	                                           sizeof( v[0] ) / sizeof( v[0][0] ) );

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

	QOpenGLVertexArrayObject::Binder raiiVAO( &o_out.m_vao );

	o_out.m_vvbo.bind();
	o_out.m_vvbo.setUsagePattern( QOpenGLBuffer::StaticDraw );
	o_out.m_vvbo.allocate( mesh.m_vertices.data(),
	                       static_cast<int>( mesh.m_vertices.size() * sizeof( mesh.m_vertices[0] ) ) );

	o_out.m_shaderProgram->bind();
	o_out.m_shaderProgram->enableAttributeArray( "position" );
	o_out.m_shaderProgram->setAttributeBuffer( "position",
	                                           GL_FLOAT,
	                                           0,
	                                           sizeof( mesh.m_vertices[0] ) / sizeof( mesh.m_vertices[0][0] ) );

	if ( mesh.hasNormals() && o_out.m_nvbo.create() )
	{
		o_out.m_nvbo.bind();
		o_out.m_nvbo.setUsagePattern( QOpenGLBuffer::StaticDraw );
		o_out.m_nvbo.allocate( mesh.m_normals.data(),
		                       static_cast<int>( mesh.m_normals.size() * sizeof( mesh.m_normals[0] ) ) );

		o_out.m_shaderProgram->enableAttributeArray( "normal" );
		o_out.m_shaderProgram->setAttributeBuffer( "normal",
		                                           GL_FLOAT,
		                                           0,
		                                           sizeof( mesh.m_normals[0] ) / sizeof( mesh.m_normals[0][0] ) );
	}

	if ( mesh.m_vindices.size() && o_out.m_ibo.create() )
	{
		o_out.m_ibo.bind();
		o_out.m_ibo.setUsagePattern( QOpenGLBuffer::StaticDraw );
		o_out.m_ibo.allocate( mesh.m_vindices.data(),
		                      static_cast<int>( mesh.m_vindices.size() * sizeof( mesh.m_vindices[0] ) ) );
	}

	return true;
}

bool GLDrawable::createFrom( const Hair& hair, QOpenGLShaderProgram& shader, GLDrawable& o_out )
{
	o_out.invalidate();

	if ( !o_out.m_vao.create() )
	{
		qWarning() << "Unable to create VAO";
		return false;
	}

	if ( !o_out.m_vvbo.create() )
	{
		qWarning() << "Unable to create vertices BO";
		o_out.invalidate();
		return false;
	}

	if ( !o_out.m_nvbo.create() )
	{
		qWarning() << "Unable to create normals BO";
		o_out.invalidate();
		return false;
	}

	if ( !o_out.m_ibo.create() )
	{
		qWarning() << "Unable to create indices BO";
		o_out.invalidate();
		return false;
	}

	const auto& strands =  hair.getStrands();
	auto nVertices = 0u;
	auto nPatches = 0u;
	for ( const auto& strand : strands )
	{
		if ( strand.m_ppos.size() > 0 )
		{
			nVertices += strand.m_ppos.size();
			nPatches += ( strand.m_ppos.size() - 1 ) * 4;
		}
	}

	std::vector<mg::Vec3D> vertices( nVertices );
	std::vector<mg::Vec3D> normals( nVertices );
	std::vector<unsigned> indices( nPatches );

	auto voffset = 0u;
	auto noffset = 0u;
	auto soffset = 0u;
	for ( const auto& strand : strands )
	{
		if ( strand.m_ppos.size() > 0 )
		{
			const auto maxPosIdx = static_cast<unsigned>( strand.m_ppos.size() - 1 );
			for ( auto i = 0u; i < maxPosIdx; ++i )
			{
				for ( auto j = 0u; j < 4; ++j )
				{
					const auto idx = 4 * i + j;
					assert( idx >= 0 && ( soffset + idx ) < indices.size() );

					indices[soffset + idx] = voffset + std::min( static_cast<unsigned>( std::max( static_cast<int>( i + j - 1 ), 0 ) ), maxPosIdx );
				}
			}
			soffset += maxPosIdx * 4;

			std::copy_n( strand.m_ppos.begin(), strand.m_ppos.size(), vertices.begin() + voffset );
			voffset += strand.m_ppos.size();

			const auto nCnt = std::min( strand.m_m1.size(), strand.m_ppos.size() );
			std::copy_n( strand.m_m1.begin(), nCnt, normals.begin() + noffset );
			noffset += strand.m_ppos.size();
		}
	}

	assert( voffset == vertices.size() );
	assert( noffset == normals.size() );
	assert( soffset == indices.size() );

	o_out.m_primitive = GL_PATCHES;
	o_out.m_nElements = static_cast<int>( indices.size() );
	o_out.m_shaderProgram = &shader;
	o_out.m_color[0] = 0.8f, o_out.m_color[1] = 0.8f, o_out.m_color[2] = 0.f;
	o_out.m_transform.identity();

	QOpenGLVertexArrayObject::Binder raiiVAO( &o_out.m_vao );
	o_out.m_shaderProgram->bind();

	o_out.m_vvbo.bind();
	o_out.m_vvbo.setUsagePattern( QOpenGLBuffer::StaticDraw );
	o_out.m_vvbo.allocate( vertices.data(),
	                       static_cast<int>( vertices.size() * sizeof( vertices[0] ) ) );

	o_out.m_shaderProgram->enableAttributeArray( "position" );
	o_out.m_shaderProgram->setAttributeBuffer( "position",
	                                           GL_FLOAT,
	                                           0,
	                                           sizeof( vertices[0] ) / sizeof( vertices[0][0] ) );

	o_out.m_nvbo.bind();
	o_out.m_nvbo.setUsagePattern( QOpenGLBuffer::StaticDraw );
	o_out.m_nvbo.allocate( normals.data(),
	                       static_cast<int>( normals.size() * sizeof( normals[0] ) ) );

	o_out.m_shaderProgram->enableAttributeArray( "normal" );
	o_out.m_shaderProgram->setAttributeBuffer( "normal",
	                                           GL_FLOAT,
	                                           0,
	                                           sizeof( normals[0] ) / sizeof( normals[0][0] ) );

	o_out.m_ibo.bind();
	o_out.m_ibo.setUsagePattern( QOpenGLBuffer::StaticDraw );
	o_out.m_ibo.allocate( indices.data(),
	                      static_cast<int>( indices.size() * sizeof( indices[0] ) ) );

	return true;
}
