#pragma once

#include "config.h"
#include <QOpenGLShaderProgram>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLBuffer>
#include <unordered_map>

class Mesh;
class Hair;

class GLShaderManager
{
public:
	static GLShaderManager* getInstance();

	inline QOpenGLShaderProgram* getShader(const std::string& shaderName) const
	{
		auto it = m_glShaderCache.find(shaderName);
		return (it != m_glShaderCache.end())? it->second.get() : nullptr;
	}

	QOpenGLShaderProgram* loadShader(const std::string& shaderName,
	                                 const std::string& vertShaderPath,
	                                 const std::string& fragShaderPath,
	                                 bool reload = false);

	QOpenGLShaderProgram* loadShader(const std::string& shaderName,
	                                 const std::string& vertShaderPath,
	                                 const std::string& fragShaderPath,
	                                 const std::string& geomShaderPath,
	                                 const std::string& tessCtrlShaderPath,
	                                 const std::string& tessEvalShaderPath,
	                                 bool reload = false);
private:
	GLShaderManager() = default;
	~GLShaderManager() = default;
	GLShaderManager(const GLShaderManager& other);
	GLShaderManager& operator =(const GLShaderManager& other);

	typedef std::unordered_map< std::string, std::unique_ptr< QOpenGLShaderProgram > > GLShaderCache;
	GLShaderCache m_glShaderCache;
};

class GLLight
{
public:
	GLLight();

	void loadToShader( QOpenGLShaderProgram& shader ) const;

	mg::Vec3D m_ambient;
	mg::Vec3D m_diffuse;
	mg::Vec3D m_specular;
	mg::Vec3D m_position;
};

class GLMaterial
{
public:
	GLMaterial();

	void loadToShader( QOpenGLShaderProgram& shader ) const;

	mg::Vec3D m_ambient;
	mg::Vec3D m_diffuse;
	mg::Vec3D m_specular;
	mg::Real m_shininess;
};

class GLDrawable
{
public:
	static bool createGrid( unsigned usize, unsigned vsize, QOpenGLShaderProgram& shader, GLDrawable& o_out );
	static bool createFrom( const Mesh& mesh, QOpenGLShaderProgram& shader, GLDrawable& o_out );
	static bool createFrom( const Hair& hair, QOpenGLShaderProgram& shader, GLDrawable& o_out );

	GLDrawable();
	~GLDrawable();

	inline bool isValid() const
	{
		return m_vao.isCreated() && m_nElements > 0 && m_shaderProgram && m_shaderProgram->isLinked();
	}

	inline const mg::Matrix4D& getTransform() const { return m_transform; }
	inline mg::Matrix4D& getTransform() { return m_transform; }
	inline QOpenGLShaderProgram* getShader() { return m_shaderProgram; }

	void invalidate();
	void draw();

private:
	QOpenGLVertexArrayObject m_vao;
	QOpenGLBuffer m_vvbo;
	QOpenGLBuffer m_nvbo;
	QOpenGLBuffer m_ibo;

	GLenum m_primitive = GL_POINTS;
	int m_nElements = 0;

	QOpenGLShaderProgram* m_shaderProgram = nullptr;
	mg::Vec3D m_color;
	mg::Matrix4D m_transform;
};

typedef std::vector< std::unique_ptr< GLDrawable > > DrawList;
