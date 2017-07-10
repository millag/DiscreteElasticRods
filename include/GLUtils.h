#ifndef GLUTILS_H
#define GLUTILS_H

#include "config.h"

#include <QOpenGLShaderProgram>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLBuffer>

#include <unordered_map>


class Mesh;


class GLShaderManager
{
public:
	static GLShaderManager* getInstance();

	GLShaderManager()
	{ }
	~GLShaderManager()
	{ }

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
	typedef std::unordered_map< std::string, std::unique_ptr< QOpenGLShaderProgram > > GLShaderCache;
	GLShaderCache m_glShaderCache;
};

class GLDrawable
{
public:
	static bool createGrid(GLDrawable& o_out, unsigned w, unsigned h, QOpenGLShaderProgram& shader);
	static bool createFrom(GLDrawable& o_out, const Mesh& mesh, QOpenGLShaderProgram& shader);

	GLDrawable()
	{ }
	~GLDrawable()
	{
		invalidate();
	}

	inline const mg::Matrix4D& getTransform() const { return m_transform; }
	inline mg::Matrix4D& getTransform() { return m_transform; }
	inline QOpenGLShaderProgram* getShader() { return m_shaderProgram; }

	inline bool isValid() const
	{
		return m_vao.isCreated() && m_nElements > 0 && m_shaderProgram && m_shaderProgram->isLinked();
	}

	void invalidate();
	void draw();

private:
	QOpenGLBuffer m_vvbo;
	QOpenGLBuffer m_nvbo;
	QOpenGLBuffer m_ibo;
	QOpenGLVertexArrayObject m_vao;
	GLenum m_glDrawMode = GL_TRIANGLES;
	int m_nElements = 0;

	QOpenGLShaderProgram* m_shaderProgram = nullptr;
	mg::Vec3D m_color;

	mg::Matrix4D m_transform;
};

#endif // GLUTILS_H
