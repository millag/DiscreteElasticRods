#ifndef HAIRYOBJECT_H
#define HAIRYOBJECT_H

#include "RenderObject.h"
#include "Hair.h"

class HairyObject : public RenderObject
{
public:
    HairyObject(const Mesh* mesh = NULL, const ngl::Mat4& transform = ngl::Mat4(), int shaderId = -1):
        RenderObject(mesh, transform, shaderId)
    {}

    void setTransform(const ngl::Mat4& t)
    {
        RenderObject::setTransform(t);
        m_hair->updateTransform(t);
    }

    void attachHair(Hair* hair, int nMaxStrands)
    {
        assert(hair != NULL);
        m_hair = hair;
        m_hair->initialize(*this, nMaxStrands);
    }

private:
    Hair* m_hair;
};

#endif // HAIRYOBJECT_H
