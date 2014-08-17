#ifndef ISURVANT_H
#define ISURVANT_H
#include "RenderObject.h"

class IServant
{

public:
    virtual ~IServant() { }
    virtual void findObjectsWithinDistance(const mg::Vec3D& pos, mg::Real dist, std::vector<RenderObject*>& o_objects) = 0;

};

#endif // ISURVANT_H
