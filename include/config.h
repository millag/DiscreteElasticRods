#ifndef CONFIG_H
#define CONFIG_H

#include "cml_config.h"
#include <cml/cml.h>

#define MULTI_THREAD
//#define DBUGG

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

    typedef matrix23f_c Matrix23D;

    typedef matrix<Real, dynamic<>, col_basis, col_major> MatrixND;

    typedef quaternion<Real, fixed<>, scalar_first> Quaternion;
}



#endif // CONFIG_H
