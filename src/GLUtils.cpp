#include "GLUtils.h"
#include "Mesh.h"

#include <QOpenGLContext>
#include <QOpenGLFunctions>

static inline GLenum getGLDrawMode(Mesh::PrimitiveMode::Enum type)
{
    switch (type)
    {
        case Mesh::PrimitiveMode::TRIANGLES:
            return GL_TRIANGLES;
        case Mesh::PrimitiveMode::LINES:
            return GL_LINE;
        case Mesh::PrimitiveMode::POINTS:
            return GL_POINTS;
    }

    return GL_POINTS;
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


void GLDrawable::invalidate()
{
    if (m_vao.isCreated())
    {
        m_vao.destroy();
    }

    if (m_vvbo.isCreated())
    {
        m_vvbo.destroy();
    }

    if (m_nvbo.isCreated())
    {
        m_nvbo.destroy();
    }

    if (m_ibo.isCreated())
    {
        m_ibo.destroy();
    }

    m_shaderProgram = nullptr;
    m_nElements = 0;
}

void GLDrawable::draw()
{
    if (!isValid())
    {
        return;
    }

    auto gl = QOpenGLContext::currentContext()->functions();

    m_shaderProgram->bind();
    m_shaderProgram->setUniformValue("color", m_color[0], m_color[1], m_color[2], 1.f);

    m_vao.bind();
    if (m_ibo.isCreated() )
    {
        gl->glDrawElements(m_glDrawMode, m_nElements, GL_UNSIGNED_INT, nullptr);
    }
    else
    {
        gl->glDrawArrays(m_glDrawMode, 0, m_nElements);
    }

    m_vao.release();
}


bool GLDrawable::createGrid(GLDrawable &o_out, unsigned w, unsigned h, QOpenGLShaderProgram& shader)
{
    o_out.invalidate();

    if (!o_out.m_vao.create())
    {
        qWarning() << "Unable to create grid VAO";
        return false;
    }

    o_out.m_vao.bind();

    if (!o_out.m_vvbo.create())
    {
        qWarning() << "Unable to create grid VBO";
        o_out.invalidate();
        return false;
    }

    // create grid geometry
    const unsigned usize = 10;
    const unsigned vsize = 10;
    mg::Vec3D bleft(-1, 0, -1);
    mg::Vec3D bright(1, 0, -1);
    mg::Vec3D tleft(-1, 0, 1);
    mg::Vec3D tright(1, 0, 1);

    std::vector<mg::Vec3D> v(2 * ( usize + vsize ));
    for (auto i = 0u; i < usize; ++i)
    {
        auto t = (float)(i) / (usize - 1);
        v[2 * i] = mg::lerp(bleft, bright, t);
        v[2 * i + 1] = mg::lerp(tleft, tright, t);
    }

    const auto offset = 2 * usize;
    for( auto i = 0u; i < vsize; ++i )
    {
        auto t = (float)(i) / (vsize - 1);
        v[offset + 2 * i] = mg::lerp(bleft, tleft, t);
        v[offset + 2 * i + 1] = mg::lerp(bright, tright, t);
    }

    o_out.m_shaderProgram = &shader;
    o_out.m_shaderProgram->bind();

    o_out.m_vvbo.bind();
    o_out.m_vvbo.setUsagePattern(QOpenGLBuffer::StaticDraw);
    o_out.m_vvbo.allocate(v.data(), sizeof(v[0]) * (int)(v.size()));

    o_out.m_shaderProgram->enableAttributeArray("position");
    o_out.m_shaderProgram->setAttributeBuffer("position", GL_FLOAT, 0, 3);

    o_out.m_vao.release();

    o_out.m_nElements = (int)(v.size());
    o_out.m_glDrawMode = GL_LINES;
    o_out.m_color[0] = o_out.m_color[1] = o_out.m_color[2] = 0.8f;
    o_out.m_transform.identity();

    return true;
}

bool GLDrawable::createFrom(GLDrawable& o_out, const Mesh& mesh, QOpenGLShaderProgram& shader)
{
    o_out.invalidate();

    if (!o_out.m_vao.create())
    {
        qWarning() << "Unable to create VAO";
        return false;
    }

    o_out.m_vao.bind();

    if (!o_out.m_vvbo.create())
    {
        qWarning() << "Unable to create vertices VBO";
        o_out.invalidate();
        return false;
    }

    o_out.m_shaderProgram = &shader;
    o_out.m_shaderProgram->bind();

    o_out.m_vvbo.bind();
    o_out.m_vvbo.setUsagePattern(QOpenGLBuffer::DynamicDraw);
    o_out.m_vvbo.allocate(mesh.m_vertices.data(), sizeof(mesh.m_vertices[0]) * (int)(mesh.m_vertices.size()));

    o_out.m_shaderProgram->enableAttributeArray("position");
    o_out.m_shaderProgram->setAttributeBuffer("position", GL_FLOAT, 0, sizeof(mesh.m_vertices[0]) / sizeof(mesh.m_vertices[0][0]));

    if (mesh.hasNormals() && o_out.m_nvbo.create())
    {
        o_out.m_nvbo.bind();
        o_out.m_nvbo.setUsagePattern(QOpenGLBuffer::DynamicDraw);
        o_out.m_nvbo.allocate(mesh.m_normals.data(), sizeof(mesh.m_normals[0]) * (int)(mesh.m_normals.size()));

        o_out.m_shaderProgram->enableAttributeArray("normal");
        o_out.m_shaderProgram->setAttributeBuffer("normaml", GL_FLOAT, 0, sizeof(mesh.m_normals[0]) / sizeof(mesh.m_normals[0][0]));
    }

    o_out.m_vao.release();

    o_out.m_glDrawMode = getGLDrawMode(mesh.getPrimitiveType());
    o_out.m_nElements = (int)(mesh.getNPrimitives());
    o_out.m_color[0] = 0.9f, o_out.m_color[1] = o_out.m_color[2] = 0.f;
    o_out.m_transform.identity();

    return true;
}
