#include "GLRenderer.h"
#include "Utils.h"
#include "Camera.h"
#include "GLUtils.h"

#include <QOpenGLContext>
#include <QOpenGLFunctions>
#include <QOpenGLFunctions_4_1_Core>
#include <QOpenGLBuffer>

#include <cassert>

namespace {

const unsigned IdxBreak = 0xffffffff;

static PrimitiveType createCircle(unsigned udiv, std::vector<mg::Vec3D>& o_vertData)
{
    o_vertData.reserve( o_vertData.size() + udiv );
    for (auto i = 0u; i < udiv; ++i)
    {
        mg::Real angle = i * mg::Constants::two_pi() / udiv;
        o_vertData.emplace_back(std::cos(angle), 0.f, std::sin(angle));
    }

    return PrimitiveType::LineLoop;
}

static PrimitiveType createBox(std::vector<mg::Vec3D>& o_vertData, std::vector<unsigned>& o_indices)
{
    const mg::Vec3D pos[] = {
        // front side
        {-0.5f, 0.5f, 0.5f},
        {-0.5f, -0.5f, 0.5f},
        {0.5f, -0.5f, 0.5f},
        {0.5f, 0.5f, 0.5f},
        // back side
        {-0.5f, 0.5f, -0.5f},
        {-0.5f, -0.5f, -0.5f},
        {0.5f, -0.5f, -0.5f},
        {0.5f, 0.5f, -0.5f},
    };

    const unsigned pindices[] = {
        0,1,2, 0,2,3,
        3,2,6, 6,7,3,
        7,6,5, 7,5,4,
        4,5,1, 4,1,0,
        0,3,4, 3,7,4,
        1,5,2, 2,5,6,
    };

    const mg::Vec3D normals[] = {
        {0.f, 0.f, 1.f},
        {1.f, 0.f, 0.f},
        {0.f, 0.f, -1.f},
        {-1.f, 0.f, 0.f},
        {0.f, 1.f, 0.f},
        {0.f, -1.f, 0.f},
    };

    o_vertData.reserve( o_vertData.size() + 2 * 6 * 4 );
    o_indices.reserve( o_indices.size() + 6 * 2 * 3);

    const auto idx = static_cast<unsigned>(o_vertData.size());
    for (int i = 0; i < mg::CountOf(pindices); ++i)
    {
        const unsigned pidx = pindices[i];
        const unsigned nidx = i / 6;
        o_vertData.push_back(pos[pidx]);
        o_vertData.push_back(normals[nidx]);
        o_indices.push_back(idx + i);
    }

    return PrimitiveType::Triangle;
}

static PrimitiveType createSphere(unsigned udiv, std::vector<mg::Vec3D>& o_vertData, std::vector<unsigned>& o_indices)
{
    // for reference see:
    // https://stackoverflow.com/questions/7687148/drawing-sphere-in-opengl-without-using-glusphere
    const unsigned n = udiv + 2;
    o_vertData.reserve( o_vertData.size() + (udiv + 2) * (1 + (udiv + 1) * 1) / 2); // aritmetic series 1+2+3+...
    o_indices.reserve(o_indices.size() + (udiv + 1) * (1 + udiv * 2) / 2); // aritmetic series: 1+3+5+ ...

    const mg::Vec3D octahedron[] = {
        {1.f, 0.f, 0.f},
        {0.f, 1.f, 0.f},
        {0.f, 0.f, 1.f},

        {1.f, 0.f, 0.f},
        {0.f, 0.f, 1.f},
        {0.f, -1.f, 0.f},

        {1.f, 0.f, 0.f},
        {0.f, -1.f, 0.f},
        {0.f, 0.f, -1.f},

        {1.f, 0.f, 0.f},
        {0.f, 0.f, -1.f},
        {0.f, 1.f, 0.f},

        {-1.f, 0.f, 0.f},
        {0.f, 0.f, 1.f},
        {0.f, 1.f, 0.f},

        {-1.f, 0.f, 0.f},
        {0.f, -1.f, 0.f},
        {0.f, 0.f, 1.f},

        {-1.f, 0.f, 0.f},
        {0.f, 0.f, -1.f},
        {0.f, -1.f, 0.f},

        {-1.f, 0.f, 0.f},
        {0.f, 1.f, 0.f},
        {0.f, 0.f, -1.f},
    };

    for (int k = 0; k < mg::CountOf(octahedron); k += 3)
    {
        o_vertData.push_back(octahedron[k]);
        for (auto i = 1u; i < n; ++i)
        {
            const auto idx1 = static_cast<unsigned>(o_vertData.size());
            // compute positions
            const auto t = static_cast<mg::Real>(i)/(n-1);
            const auto v1 = mg::lerp(octahedron[k], octahedron[k+1], t);
            const auto v2 = mg::lerp(octahedron[k], octahedron[k+2], t);
            for (auto j = 0u; j <= i; ++j)
            {
                const auto k = static_cast<mg::Real>(j)/i;
                o_vertData.push_back(mg::lerp(v1, v2, k));
            }
            // compute faces for the above positions
            const auto idx2 = static_cast<unsigned>(o_vertData.size());
            for (auto j = 1u; j <= i; ++j)
            {
                o_indices.push_back(idx2-j-1);
                o_indices.push_back(idx2-j);
                o_indices.push_back(idx1-j);

                if (idx2-j-1 <= idx1)
                {
                    break;
                }

                o_indices.push_back(idx1-j);
                o_indices.push_back(idx1-j-1);
                o_indices.push_back(idx2-j-1);
            }
        }
    }

    for (auto& p : o_vertData)
    {
        p.normalize();
    }

    return PrimitiveType::Triangle;
}

static PrimitiveType createCone(unsigned udiv, std::vector<mg::Vec3D>& o_vertData, std::vector<unsigned>& o_indices)
{

    o_vertData.reserve( o_vertData.size() + udiv * 3 + 1 );
    o_indices.reserve( o_indices.size() + udiv * 3 * 2);

    const mg::Vec3D apex(0.f, 1.f, 0.f);
    for (auto i = 0u; i < udiv; ++i)
    {
        o_vertData.emplace_back(apex);
        o_vertData.emplace_back(0.f, 0.f, 0.f); // reserve for normal
        const mg::Real angle = i * mg::Constants::two_pi() / udiv;
        o_vertData.emplace_back(std::cos(angle), 0.f, std::sin(angle));
        o_vertData.emplace_back(0.f, 0.f, 0.f); // reserve for normal
    }

    for (auto i = 0u; i < udiv ; ++i)
    {
        // the apex
        const auto pidx1 = 2*i;
        // vertex
        const auto pidx2 = 2*i+1;
        // next vertex
        const auto pidx3 = 2*((i+1)%(udiv)) + 1;

        o_indices.push_back(pidx2);
        o_indices.push_back(pidx1);
        o_indices.push_back(pidx3);

        // compute normals
        const auto normal = mg::cross(o_vertData[2*pidx1] - o_vertData[2*pidx2], o_vertData[2*pidx3] - o_vertData[2*pidx2]);
        o_vertData[2*pidx1+1] += normal;
        o_vertData[2*pidx2+1] += normal;
    }

    // center
    const auto idx = 2 * udiv;
    o_vertData.emplace_back(0.f, 0.f, 0.f);
    o_vertData.emplace_back(0.f, -1.f, 0.f);
    for (auto i = 0u; i < udiv; ++i)
    {
        const mg::Real angle = i * mg::Constants::two_pi() / udiv;
        o_vertData.emplace_back(std::cos(angle), 0.f, std::sin(angle));
        o_vertData.emplace_back(0.f, -1.f, 0.f);

        o_indices.push_back(idx);
        o_indices.push_back(idx + 1 + i);
        o_indices.push_back(idx + 1 + ((i+1)%udiv));
    }

    return PrimitiveType::Triangle;
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
        gl->glEnable(GL_CULL_FACE);
//        gl->glEnable(GL_PRIMITIVE_RESTART_FIXED_INDEX);
        gl->glFrontFace(GL_CCW);
//        gl->glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

        auto shaderMan = GLShaderManager::getInstance();
        shaderMan->loadShader("Constant",
                              "shaders/ConstantVert.glsl",
                              "shaders/ConstantFrag.glsl");

        shaderMan->loadShader("Phong",
                              "shaders/PhongVert.glsl",
                              "shaders/PhongFrag.glsl");

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

        m_vMatrix.zero();
        m_pMatrix.zero();
        m_vpMatrix.zero();
        m_transformStack.push(mg::Matrix4D().identity());
    }

    return isValid();
}

void GLRenderer::setCamera(Camera &cam)
{
    m_vMatrix = cam.getVMatrix();
    m_pMatrix = cam.getPMatrix();
    m_vpMatrix = cam.getVPMatrix();
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
    auto shader = shaderMan->getShader("Constant");
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
    shader->setUniformValue( "mvp", *reinterpret_cast<const GLMatrix4x4*>(m_vpMatrix.data()) );
    const auto& color = getColor();
    shader->setUniformValue("color", color[0], color[1], color[2], 1.f);

    shader->enableAttributeArray("position");
    shader->setAttributeBuffer("position", GL_FLOAT, 0, 3);

    GLenum primType = (closed)? GL_LINE_LOOP : GL_LINE_STRIP;

    auto gl = m_context->functions();
    gl->glDrawArrays(primType, 0, cnt);
}

void GLRenderer::circle(const mg::Vec3D& center, const mg::Vec3D& normal, mg::Real radius)
{
    assert( isValid() );

    auto shaderMan = GLShaderManager::getInstance();
    auto shader = shaderMan->getShader("Constant");
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
    mg::matrix_aim_at(tm , center, center + normal, mg::axis_order_yzx);
    mg::matrix_set_x_basis_vector(tm, mg::matrix_get_x_basis_vector(tm) * radius);
    mg::matrix_set_y_basis_vector(tm, mg::matrix_get_y_basis_vector(tm) * radius);
    mg::matrix_set_z_basis_vector(tm, mg::matrix_get_z_basis_vector(tm) * radius);
    tm = m_vpMatrix * tm;

    QOpenGLVertexArrayObject::Binder raiivao(&m_defaultVAO);
    vbo.bind();
    vbo.allocate(pos.data(), static_cast<GLsizei>(pos.size() * sizeof(pos[0])));

    shader->bind();
    shader->setUniformValue("mvp", *reinterpret_cast<const GLMatrix4x4*>(tm.data()));
    const auto& color = getColor();
    shader->setUniformValue("color", color[0], color[1], color[2], 1.f);

    shader->enableAttributeArray("position");
    shader->setAttributeBuffer("position", GL_FLOAT, 0, 3);

    auto gl = m_context->functions();
    gl->glDrawArrays(primType, 0, static_cast<GLsizei>(pos.size()));
}

void GLRenderer::box(const mg::Vec3D& center, const mg::Vec3D& zdir, const mg::Vec3D& scale)
{
    assert( isValid() );

    auto shaderMan = GLShaderManager::getInstance();
    auto shader = shaderMan->getShader("Phong");
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

    std::vector<mg::Vec3D> vertData;
    std::vector<unsigned> indices;
    auto primType = static_cast<GLenum>(createBox(vertData, indices));

    mg::Matrix4D tm;
    mg::matrix_aim_at(tm , center, center + zdir, mg::axis_order_zxy);
    mg::matrix_set_x_basis_vector(tm, mg::matrix_get_x_basis_vector(tm) * scale[0]);
    mg::matrix_set_y_basis_vector(tm, mg::matrix_get_y_basis_vector(tm) * scale[1]);
    mg::matrix_set_z_basis_vector(tm, mg::matrix_get_z_basis_vector(tm) * scale[2]);
    mg::Matrix4D mv = m_vMatrix * tm;
    mg::Matrix4D mvp = m_vpMatrix * tm;

    QOpenGLVertexArrayObject::Binder raiiVAO(&m_defaultVAO);
    vbo.bind();
    vbo.allocate(vertData.data(), static_cast<GLsizei>(vertData.size() * sizeof(vertData[0])));

    shader->bind();
    shader->setUniformValue( "mv", *reinterpret_cast<const GLMatrix4x4*>(mv.data()));
    shader->setUniformValue( "mvp", *reinterpret_cast<const GLMatrix4x4*>(mvp.data()));

    m_headLight.loadToShader(*shader);
    m_defaultMtl.loadToShader(*shader);

    shader->enableAttributeArray("position");
    shader->setAttributeBuffer("position", GL_FLOAT, 0, 3, 2 * sizeof(vertData[0]));
    shader->enableAttributeArray("normal");
    shader->setAttributeBuffer("normal", GL_FLOAT, sizeof(vertData[0]), 3, 2 * sizeof(vertData[0]));

    auto gl = m_context->functions();
    gl->glDrawElements(primType, static_cast<GLsizei>(indices.size()), GL_UNSIGNED_INT, indices.data());
}

void GLRenderer::sphere(const mg::Vec3D& center, mg::Real radius)
{
    assert( isValid() );

    auto shaderMan = GLShaderManager::getInstance();
    auto shader = shaderMan->getShader("Phong");
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

    std::vector<mg::Vec3D> vertData;
    std::vector<unsigned> indices;
    auto primType = static_cast<GLenum>(createSphere(2, vertData, indices));

    mg::Matrix4D tm;
    mg::matrix_uniform_scale(tm, radius);
    mg::matrix_set_translation(tm, center);
    mg::Matrix4D mv = m_vMatrix * tm;
    mg::Matrix4D mvp = m_vpMatrix * tm;

    QOpenGLVertexArrayObject::Binder raiiVAO(&m_defaultVAO);
    vbo.bind();
    vbo.allocate(vertData.data(), static_cast<GLsizei>(vertData.size() * sizeof(vertData[0])));

    shader->bind();
    shader->setUniformValue( "mv", *reinterpret_cast<const GLMatrix4x4*>(mv.data()));
    shader->setUniformValue( "mvp", *reinterpret_cast<const GLMatrix4x4*>(mvp.data()));

    m_headLight.loadToShader(*shader);
    m_defaultMtl.loadToShader(*shader);

    shader->enableAttributeArray("position");
    shader->setAttributeBuffer("position", GL_FLOAT, 0, 3, sizeof(vertData[0]));
    shader->enableAttributeArray("normal");
    shader->setAttributeBuffer("normal", GL_FLOAT,0, 3, sizeof(vertData[0]));

    auto gl = m_context->functions();
    gl->glDrawElements(primType, static_cast<GLsizei>(indices.size()), GL_UNSIGNED_INT, indices.data());
}

void GLRenderer::cone(const mg::Vec3D& center, const mg::Vec3D& updir, mg::Real height, mg::Real radius)
{
    assert( isValid() );

    auto shaderMan = GLShaderManager::getInstance();
    auto shader = shaderMan->getShader("Phong");
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

    std::vector<mg::Vec3D> vertData;
    std::vector<unsigned> indices;
    auto primType = static_cast<GLenum>(createCone(6, vertData, indices));

    mg::Matrix4D tm;
    mg::matrix_aim_at(tm , center, center + updir, mg::axis_order_yxz);
    mg::matrix_set_x_basis_vector(tm, mg::matrix_get_x_basis_vector(tm) * radius);
    mg::matrix_set_y_basis_vector(tm, mg::matrix_get_y_basis_vector(tm) * height);
    mg::matrix_set_z_basis_vector(tm, mg::matrix_get_z_basis_vector(tm) * radius);
    mg::Matrix4D mv = m_vMatrix * tm;
    mg::Matrix4D mvp = m_vpMatrix * tm;

    QOpenGLVertexArrayObject::Binder raiiVAO(&m_defaultVAO);
    vbo.bind();
    vbo.allocate(vertData.data(), static_cast<GLsizei>(vertData.size() * sizeof(vertData[0])));

    shader->bind();
    shader->setUniformValue( "mv", *reinterpret_cast<const GLMatrix4x4*>(mv.data()));
    shader->setUniformValue( "mvp", *reinterpret_cast<const GLMatrix4x4*>(mvp.data()));

    m_headLight.loadToShader(*shader);
    m_defaultMtl.loadToShader(*shader);

    shader->enableAttributeArray("position");
    shader->setAttributeBuffer("position", GL_FLOAT, 0, 3, 2 * sizeof(vertData[0]));
    shader->enableAttributeArray("normal");
    shader->setAttributeBuffer("normal", GL_FLOAT, sizeof(vertData[0]), 3, 2 * sizeof(vertData[0]));

    auto gl = m_context->functions();
    gl->glDrawElements(primType, static_cast<GLsizei>(indices.size()), GL_UNSIGNED_INT, indices.data());
}
