#ifndef ISURVANT_H
#define ISURVANT_H
#include <ngl/Vec4.h>
#include "RenderObject.h"

class IServant
{

public:
    virtual ~IServant() { }
    virtual void findObjectsWithinDistance(const ngl::Vec4& pos, ngl::Real dist, std::vector<RenderObject*>& o_objects) = 0;

};

#endif // ISURVANT_H
