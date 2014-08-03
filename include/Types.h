#ifndef TYPES_H
#define TYPES_H

#include "cml_config.h"
#include <cml/cml.h>

namespace mg
{
    using  namespace cml;

    typedef float Real;

    typedef constantsf Constants;

    typedef vector2f Vec2D;
    typedef vector3f Vec3D;
    typedef vector4f Vec4D;

    typedef matrix22f_c Matrix2D;
    typedef matrix33f_c Matrix3D;
    typedef matrix44f_c Matrix4D;

    typedef matrix32f_r Matrix32D;

    typedef matrix<Real, dynamic<>, col_basis> MatrixND;

    typedef quaternion<Real, fixed<>, scalar_first> Quaternion;
}



#endif // TYPES_H
