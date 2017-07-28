#ifndef GLRENDERER_H
#define GLRENDERER_H

#include "config.h"

#include <QOpenGLVertexArrayObject>

#include <stack>

class QOpenGLContext;


typedef GLfloat GLMatrix2x2[2][2];
typedef GLfloat GLMatrix3x3[3][3];
typedef GLfloat GLMatrix4x4[4][4];

enum class PrimitiveType : uint16_t
{
    Point = GL_POINTS,
    Line = GL_LINES,
    LineStrip = GL_LINE_STRIP,
    LineLoop = GL_LINE_LOOP,
    Triangle = GL_TRIANGLES,
    TriangleStrip = GL_TRIANGLE_STRIP,
    TriangleFan = GL_TRIANGLE_FAN,
    Quad = GL_QUADS,
    QuadStrip = GL_QUAD_STRIP,
    Polygon = GL_POLYGON,
};

enum class PickMode : uint8_t
{
    NonPickable,
    Pickable,
};

using PickName = uint32_t;

class GLRenderer
{
public:
    GLRenderer() = default;
    ~GLRenderer() = default;

    inline bool isValid() const { return m_context  != nullptr; }
    bool initialize();

    /// primitive draw
    void beginDrawable(PickMode mode = PickMode::NonPickable, PickName pickName = 0);
    void endDrawable();

    void line(const mg::Vec3D& start, const mg::Vec3D& end);
    void polyline(const mg::Vec3D pos[], int cnt, bool closed = false);
    void circle(const mg::Vec3D& base, const mg::Vec3D& dir, mg::Real radius);
    void box(const mg::Vec3D& base, const mg::Vec3D& dir, const mg::Vec3D& scale);
    void sphere(const mg::Vec3D& base, mg::Real radius);
    void cone(const mg::Vec3D& base, const mg::Vec3D& dir, mg::Real height, mg::Real radius);

    /// transform stack
    inline const mg::Matrix4D& getTransform() const { return m_transformStack.top(); }
    inline void setTransform(const mg::Matrix4D& tm) { m_transformStack.top() = tm; }
    inline void setToIdentityTransform() { m_transformStack.top().identity(); }
    inline void multiplyTransform(const mg::Matrix4D& tm) { m_transformStack.top() *= tm; }
    inline void pushTransform() { m_transformStack.push( m_transformStack.top() ); }
    inline void popTransform()
    {
        // there should always be a top transform available
        if (m_transformStack.size() < 2)
        {
            return;
        }
        m_transformStack.pop();
    }

protected:
    /// renderer's gl context
    QOpenGLContext* m_context = nullptr;
    /// transform stack
    std::stack<mg::Matrix4D> m_transformStack;
    /// dummy VAO used for immediate mode rendering
    QOpenGLVertexArrayObject m_defaultVAO;
};

#endif // GLRENDERER_H
