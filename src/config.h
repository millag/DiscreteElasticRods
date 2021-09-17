#pragma once

#include "cml_config.h"
#include <cml/cml.h>

#include <cmath>
#include <memory>
#include <cassert>

#define UNUSED_VALUE( name ) ( void )( name )
#define STRINGIFY( name ) #name

#ifndef MULTI_THREADING_ON
    #define MULTI_THREADING_ON 0
#endif

#if MULTI_THREADING_ON != 0
    #ifdef _WIN32
        #define OMP_PARALLEL_LOOP __pragma( omp parallel for )
    #else
        #define OMP_PARALLEL_LOOP _Pragma( STRINGIFY( omp parallel for ) )
    #endif
#else
    #define OMP_PARALLEL_LOOP
#endif

#define NON_COPYABLE( TypeName ) \
    TypeName( const TypeName& other ) = delete; \
    TypeName& operator=( const TypeName& other ) = delete;

namespace mg
{
using namespace cml;

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

const static unsigned SEC = 1000;
const static unsigned UPS = 30;
const static unsigned FPS = 30;
const static unsigned UPDATERATE = SEC / UPS;
const static unsigned REFRESHRATE = SEC / FPS;

const static Real ERR = 1e-6f;
const static Real THRESHODL = 0.00001f;
const static mg::Vec3D Origin(0,0,0);
const static mg::Vec3D Ox(1,0,0);
const static mg::Vec3D Oy(0,1,0);
const static mg::Vec3D Oz(0,0,1);

const static mg::Vec3D Gravity(0,-9.81f,0);
}
