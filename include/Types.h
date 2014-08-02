#ifndef TYPES_H
#define TYPES_H

#include "cml_config.h"
#include <cml/cml.h>

namespace mg
{
    typedef float Real;

    typedef cml::constantsf Constants;

    typedef cml::vector2f Vec2D;
    typedef cml::vector3f Vec3D;
    typedef cml::vector4f Vec4D;

    typedef cml::matrix22f_c Matrix2D;
    typedef cml::matrix33f_c Matrix3D;
    typedef cml::matrix44f_c Matrix4D;

    typedef cml::matrix<Real, cml::dynamic<>, cml::col_basis> MatrixND;

    typedef cml::quaternion<Real, cml::fixed<>, cml::scalar_first> Quaternion;
}



#endif // TYPES_H
