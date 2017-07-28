#include "GLRenderer.h"
#include "Utils.h"
#include "GLUtils.h"

#include <QOpenGLContext>
#include <QOpenGLFunctions>
#include <QOpenGLFunctions_4_1_Core>
#include <QOpenGLBuffer>

#include <cassert>

namespace {

const unsigned IdxBreak = 0xffffffff;

static PrimitiveType createCircle(unsigned udiv, std::vector<mg::Vec3D>& o_pos)
{
    o_pos.reserve( udiv );
    for (auto i = 0u; i < udiv; ++i)
    {
        mg::Real angle = i * mg::Constants::two_pi() / udiv;
        o_pos.emplace_back(std::cos(angle), 0.f, std::sin(angle));
    }

    return PrimitiveType::LineLoop;
}

static PrimitiveType createBox(std::vector<mg::Vec3D>& o_pos, std::vector<unsigned>& o_idx)
{
    o_pos.reserve( 4 + 4 );
    o_idx.reserve( 3 * 2 * 2 * 3);

    // front size
    o_pos.emplace_back(-0.5f, 0.5f, 0.5f);
    o_pos.emplace_back(-0.5f, -0.5f, 0.5f);
    o_pos.emplace_back(0.5f, -0.5f, 0.5f);
    o_pos.emplace_back(0.5f, 0.5f, 0.5f);

    // back side
    o_pos.emplace_back(-0.5f, 0.5f, -0.5f);
    o_pos.emplace_back(-0.5f, -0.5f, -0.5f);
    o_pos.emplace_back(0.5f, -0.5f, -0.5f);
    o_pos.emplace_back(0.5f, 0.5f, -0.5f);

    const unsigned idx[] = {
        0,1,2, 0,2,3,
        3,2,6, 6,7,3,
        7,6,5, 7,5,4,
        4,5,1, 4,1,0,
        0,3,4, 3,7,4,
        1,5,2, 2,5,6,
    };
    o_idx.assign(idx, idx + mg::CountOf(idx));

    return PrimitiveType::Triangle;
}

static PrimitiveType createSphere(unsigned udiv, std::vector<mg::Vec3D>& o_pos, std::vector<unsigned>& o_idx)
{
    // for reference see:
    // https://stackoverflow.com/questions/7687148/drawing-sphere-in-opengl-without-using-glusphere
    const unsigned n = udiv + 2;
    o_pos.reserve((udiv + 2) * (1 + (udiv + 1) * 1) / 2); // aritmetic series 1+2+3+...
    o_idx.reserve((udiv + 1) * (1 + udiv * 2) / 2); // aritmetic series: 1+3+5+ ...

    mg::Vec3D octahedron[] = {
        {0.5f, 0.f, 0.f},
        {0.f, 0.5f, 0.f},
        {0.f, 0.f, 0.5f},

        {0.5f, 0.f, 0.f},
        {0.f, -0.5f, 0.f},
        {0.f, 0.f, 0.5f},

        {0.5f, 0.f, 0.f},
        {0.f, -0.5f, 0.f},
        {0.f, 0.f, -0.5f},

        {0.5f, 0.f, 0.f},
        {0.f, 0.5f, 0.f},
        {0.f, 0.f, -0.5f},

        {-0.5f, 0.f, 0.f},
        {0.f, 0.5f, 0.f},
        {0.f, 0.f, 0.5f},

        {-0.5f, 0.f, 0.f},
        {0.f, -0.5f, 0.f},
        {0.f, 0.f, 0.5f},

        {-0.5f, 0.f, 0.f},
        {0.f, -0.5f, 0.f},
        {0.f, 0.f, -0.5f},

        {-0.5f, 0.f, 0.f},
        {0.f, 0.5f, 0.f},
        {0.f, 0.f, -0.5f},
    };

    for (int k = 0; k < mg::CountOf(octahedron); k += 3)
    {
        o_pos.push_back(octahedron[k]);
        for (auto i = 1u; i < n; ++i)
        {
            const auto idx1 = static_cast<unsigned>(o_pos.size());
            // compute positions
            const auto t = static_cast<mg::Real>(i)/(n-1);
            const auto v1 = mg::lerp(octahedron[k], octahedron[k+1], t);
            const auto v2 = mg::lerp(octahedron[k], octahedron[k+2], t);
            for (auto j = 0u; j <= i; ++j)
            {
                const auto k = static_cast<mg::Real>(j)/i;
                o_pos.push_back(mg::lerp(v1, v2, k));
            }
            // compute faces for the above positions
            const auto idx2 = static_cast<unsigned>(o_pos.size());
            for (auto j = 1u; j <= i; ++j)
            {
                o_idx.push_back(idx1-j);
                o_idx.push_back(idx2-j);
                o_idx.push_back(idx2-j-1);

                if (idx2-j-1 <= idx1)
                {
                    break;
                }

                o_idx.push_back(idx1-j);
                o_idx.push_back(idx2-j-1);
                o_idx.push_back(idx1-j-1);
            }
        }
    }

    for (auto& p : o_pos)
    {
        p.normalize();
    }

    return PrimitiveType::Triangle;
}

static PrimitiveType createCone(unsigned udiv, std::vector<mg::Vec3D>& o_pos, std::vector<unsigned>& o_idx)
{
    o_pos.reserve( udiv + 2 );
    o_idx.reserve( 2 * udiv + 2);

    o_pos.emplace_back(0.f, 0.f, 0.f);
    o_pos.emplace_back(0.f, 1.f, 0.f);
    for (auto i = 0u; i < udiv; ++i)
    {
        mg::Real angle = i * mg::Constants::two_pi() / udiv;
        o_pos.emplace_back(std::cos(angle), 0.f, std::sin(angle));
    }

    o_idx.push_back(0);
    for (auto i = 0u; i < udiv ; ++i)
    {
        o_idx.push_back(i + 2);
    }
    o_idx.push_back(2);

    o_idx.push_back(IdxBreak);
    o_idx.push_back(1);
    for (auto i = 0u; i < udiv ; ++i)
    {
        o_idx.push_back(i + 2);
    }
    o_idx.push_back(2);

    return PrimitiveType::TriangleFan;
}

}


bool GLRenderer::initialize()
{
    if (isValid())
    {
        return false;
    }

    m_context = QOpenGLContext::currentContext();
    if (m_context)
    {
//        auto gl = m_context->functions();
        auto gl = static_cast<QOpenGLFunctions_4_1_Core*>(m_context->versionFunctions());
        qDebug() << "Initializing OpenGL render";
        qDebug() << "Vendor:" << reinterpret_cast<const char*>( gl->glGetString( GL_VENDOR ) );
        qDebug() << "Renderer:" << reinterpret_cast<const char*>( gl->glGetString( GL_RENDERER ) );
        qDebug() << "Version:" << reinterpret_cast<const char*>( gl->glGetString( GL_VERSION ) );
        qDebug() << "Profile:" << m_context->format().profile();

        gl->glClearColor(0.2f, 0.2f, 0.2f, 1.0f);
        gl->glEnable(GL_DEPTH_TEST);
        gl->glEnable(GL_LINE_SMOOTH);
        gl->glFrontFace(GL_CCW);
        gl->glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

        auto shaderMan = GLShaderManager::getInstance();
        shaderMan->loadShader("Color",
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

        m_transformStack.push(mg::Matrix4D().identity() );
    }

    return isValid();
}

void GLRenderer::beginDrawable(PickMode mode , PickName pickName)
{
    if (!m_defaultVAO.isCreated())
    {
        m_defaultVAO.create();
    }
}

void GLRenderer::endDrawable()
{
}

void GLRenderer::line(const mg::Vec3D &start, const mg::Vec3D &end)
{
    assert( isValid() );

    const mg::Vec3D pos[2] = {start, end};
    polyline(pos, mg::CountOf(pos));
}

void GLRenderer::polyline(const mg::Vec3D pos[], int cnt, bool closed)
{
    assert( isValid() );

    auto shaderMan = GLShaderManager::getInstance();
    auto shader = shaderMan->getShader("Color");
    if(!shader)
    {
        return;
    }

    QOpenGLBuffer vbo(QOpenGLBuffer::VertexBuffer);
    vbo.setUsagePattern(QOpenGLBuffer::StaticDraw);
    if (!vbo.create())
    {
        return;
    }

    QOpenGLVertexArrayObject::Binder raiivao(&m_defaultVAO);
    vbo.bind();
    vbo.allocate(pos, cnt * sizeof(pos[0]));

    shader->bind();
    shader->setUniformValue("color", 0.f, 1.f, 0.f);
    shader->setUniformValue( "mvp", *reinterpret_cast<const GLMatrix4x4*>(getTransform().data()) );

    shader->enableAttributeArray("position");
    shader->setAttributeBuffer("position", GL_FLOAT, 0, 3);

    GLenum primType = (closed)? GL_LINE_LOOP : GL_LINE_STRIP;

    auto gl = m_context->functions();
    gl->glDrawArrays(primType, 0, cnt);
}

void GLRenderer::circle(const mg::Vec3D& base, const mg::Vec3D& dir, mg::Real radius)
{
    assert( isValid() );

    auto shaderMan = GLShaderManager::getInstance();
    auto shader = shaderMan->getShader("Color");
    if(!shader)
    {
        return;
    }

    QOpenGLBuffer vbo(QOpenGLBuffer::VertexBuffer);
    vbo.setUsagePattern(QOpenGLBuffer::StaticDraw);
    if (!vbo.create())
    {
        return;
    }

    std::vector<mg::Vec3D> pos;
    auto primType = static_cast<GLenum>(createCircle(30, pos));

    mg::Matrix4D tm;
    mg::matrix_aim_at(tm , base, base + dir, mg::axis_order_yzx);
    mg::matrix_set_x_basis_vector(tm, mg::matrix_get_x_basis_vector(tm) * radius);
    mg::matrix_set_y_basis_vector(tm, mg::matrix_get_y_basis_vector(tm) * radius);
    mg::matrix_set_z_basis_vector(tm, mg::matrix_get_z_basis_vector(tm) * radius);
    tm = getTransform() * tm;

    QOpenGLVertexArrayObject::Binder raiivao(&m_defaultVAO);
    vbo.bind();
    vbo.allocate(pos.data(), static_cast<GLsizei>(pos.size() * sizeof(pos[0])));

    shader->bind();
    shader->setUniformValue("color", 0.f, 1.f, 0.f);
    shader->setUniformValue("mvp", *reinterpret_cast<const GLMatrix4x4*>(tm.data()));

    shader->enableAttributeArray("position");
    shader->setAttributeBuffer("position", GL_FLOAT, 0, 3);

    auto gl = m_context->functions();
    gl->glDrawArrays(primType, 0, static_cast<GLsizei>(pos.size()));
}

void GLRenderer::box(const mg::Vec3D& base, const mg::Vec3D& dir, const mg::Vec3D& scale)
{
    assert( isValid() );

    auto shaderMan = GLShaderManager::getInstance();
    auto shader = shaderMan->getShader("Color");
    if(!shader)
    {
        return;
    }

    QOpenGLBuffer vbo(QOpenGLBuffer::VertexBuffer);
    vbo.setUsagePattern(QOpenGLBuffer::StaticDraw);
    if (!vbo.create())
    {
        return;
    }

    std::vector<mg::Vec3D> pos;
    std::vector<unsigned> idx;
    auto primType = static_cast<GLenum>(createBox(pos, idx));

    mg::Matrix4D tm;
    mg::matrix_aim_at(tm , base, base + dir, mg::axis_order_zxy);
    mg::matrix_set_x_basis_vector(tm, mg::matrix_get_x_basis_vector(tm) * scale[0]);
    mg::matrix_set_y_basis_vector(tm, mg::matrix_get_y_basis_vector(tm) * scale[1]);
    mg::matrix_set_z_basis_vector(tm, mg::matrix_get_z_basis_vector(tm) * scale[2]);
    tm = getTransform() * tm;

    QOpenGLVertexArrayObject::Binder raiiVAO(&m_defaultVAO);
    vbo.bind();
    vbo.allocate(pos.data(), static_cast<GLsizei>(pos.size() * sizeof(pos[0])));

    shader->bind();
    shader->setUniformValue("color", 0.f, 0.f, 1.f);
    shader->setUniformValue( "mvp", *reinterpret_cast<const GLMatrix4x4*>(tm.data()));

    shader->enableAttributeArray("position");
    shader->setAttributeBuffer("position", GL_FLOAT, 0, 3);

    auto gl = m_context->functions();
    gl->glDrawElements(primType, static_cast<GLsizei>(idx.size()), GL_UNSIGNED_INT, idx.data());
}

void GLRenderer::sphere(const mg::Vec3D& base, mg::Real radius)
{
    assert( isValid() );

    auto shaderMan = GLShaderManager::getInstance();
    auto shader = shaderMan->getShader("Color");
    if(!shader)
    {
        return;
    }

    QOpenGLBuffer vbo(QOpenGLBuffer::VertexBuffer);
    vbo.setUsagePattern(QOpenGLBuffer::StaticDraw);
    if (!vbo.create())
    {
        return;
    }

    std::vector<mg::Vec3D> pos;
    std::vector<unsigned> idx;
    auto primType = static_cast<GLenum>(createSphere(2, pos, idx));

    mg::Matrix4D tm;
    mg::matrix_uniform_scale(tm, radius);
    tm = getTransform() * tm;

    QOpenGLVertexArrayObject::Binder raiiVAO(&m_defaultVAO);
    vbo.bind();
    vbo.allocate(pos.data(), static_cast<GLsizei>(pos.size() * sizeof(pos[0])));

    shader->bind();
    shader->setUniformValue("color", 0.f, 0.f, 1.f);
    shader->setUniformValue( "mvp", *reinterpret_cast<const GLMatrix4x4*>(tm.data()));

    shader->enableAttributeArray("position");
    shader->setAttributeBuffer("position", GL_FLOAT, 0, 3);

    auto gl = m_context->functions();
    gl->glDrawElements(primType, static_cast<GLsizei>(idx.size()), GL_UNSIGNED_INT, idx.data());
}

void GLRenderer::cone(const mg::Vec3D& base, const mg::Vec3D& dir, mg::Real height, mg::Real radius)
{
    assert( isValid() );

    auto shaderMan = GLShaderManager::getInstance();
    auto shader = shaderMan->getShader("Color");
    if(!shader)
    {
        return;
    }

    QOpenGLBuffer vbo(QOpenGLBuffer::VertexBuffer);
    vbo.setUsagePattern(QOpenGLBuffer::StaticDraw);
    if (!vbo.create())
    {
        return;
    }

    std::vector<mg::Vec3D> pos;
    std::vector<unsigned> idx;
    auto primType = static_cast<GLenum>(createCone(6, pos, idx));

    mg::Matrix4D tm;
    mg::matrix_aim_at(tm , base, base + dir, mg::axis_order_yzx);
    mg::matrix_set_x_basis_vector(tm, mg::matrix_get_x_basis_vector(tm) * radius);
    mg::matrix_set_y_basis_vector(tm, mg::matrix_get_y_basis_vector(tm) * height);
    mg::matrix_set_z_basis_vector(tm, mg::matrix_get_z_basis_vector(tm) * radius);
    tm = getTransform() * tm;

    QOpenGLVertexArrayObject::Binder raiiVAO(&m_defaultVAO);
    vbo.bind();
    vbo.allocate(pos.data(), static_cast<GLsizei>(pos.size() * sizeof(pos[0])));

    shader->bind();
    shader->setUniformValue("color", 0.f, 0.f, 1.f);
    shader->setUniformValue( "mvp", *reinterpret_cast<const GLMatrix4x4*>(tm.data()));

    shader->enableAttributeArray("position");
    shader->setAttributeBuffer("position", GL_FLOAT, 0, 3);

    auto gl = m_context->functions();
    gl->glEnable(GL_PRIMITIVE_RESTART_FIXED_INDEX);
    gl->glDrawElements(primType, static_cast<GLsizei>(idx.size()), GL_UNSIGNED_INT, idx.data());
}
