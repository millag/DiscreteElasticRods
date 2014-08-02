#ifndef UTILS_H
#define UTILS_H

#include "Types.h"
#include "ngl/Vec4.h"

namespace mg
{
    template<class T>
    T clamp(T v, T a, T b)
    {
        assert(a < b);
        return std::max(a, std::min(v, b));
    }


    const static unsigned short SEC = 1000;
    const static unsigned short UPS = 30;
    const static unsigned short FPS = 30;
    const static unsigned short UPDATERATE = SEC / UPS;
    const static unsigned short REFRESHRATE = SEC / FPS;

    const static Real ERR = 1e-6f;
    const static ngl::Vec4 EX(1, 0, 0, 0);
    const static ngl::Vec4 EY(0, 1, 0, 0);
    const static ngl::Vec4 EZ(0, 0, 1, 0);
    const static ngl::Vec4 EW(0, 0, 0, 1);

    const static ngl::Vec4 G(0, -9.81, 0, 0);
    const static mg::Vec4D GRAVITY(0, -9.81, 0, 0);

    mg::Real getSign(mg::Real value);
    mg::Real randf(mg::Real min = 0.0, mg::Real max = 1.0);
    ngl::Vec4 genRandPointInBox(mg::Real bBoxMin = -1.0, mg::Real bBoxMax = 1.0);
    ngl::Vec4 genRandPointOnSphere(mg::Real radius = 1.0, const ngl::Vec4& center = ngl::Vec4());
    ngl::Vec4 genRandPointOnDisk(mg::Real radius = 1.0, const ngl::Vec4& center = ngl::Vec4());

    void truncate(ngl::Vec4& io_v, mg::Real maxLength);
    ngl::Vec4 faceforward(const ngl::Vec4& n, const ngl::Vec4& v);
    ngl::Vec4 reflect(const ngl::Vec4& v, const ngl::Vec4& n);
}

#endif // UTILS_H
