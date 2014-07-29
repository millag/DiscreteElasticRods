#ifndef UTILS_H
#define UTILS_H

#include "ngl/Vec4.h"

namespace utils
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

    const static ngl::Real ERR = 1e-6f;
    const static ngl::Vec4 EX(1, 0, 0, 0);
    const static ngl::Vec4 EY(0, 1, 0, 0);
    const static ngl::Vec4 EZ(0, 0, 1, 0);
    const static ngl::Vec4 EW(0, 0, 0, 1);

    const static ngl::Vec4 G(0, -9.81, 0, 0);

    ngl::Real getSign(ngl::Real value);
    ngl::Real randf(ngl::Real min = 0.0, ngl::Real max = 1.0);
    ngl::Vec4 genRandPointInBox(ngl::Real bBoxMin = -1.0, ngl::Real bBoxMax = 1.0);
    ngl::Vec4 genRandPointOnSphere(ngl::Real radius = 1.0, const ngl::Vec4& center = ngl::Vec4());
    ngl::Vec4 genRandPointOnDisk(ngl::Real radius = 1.0, const ngl::Vec4& center = ngl::Vec4());

    void truncate(ngl::Vec4& io_v, ngl::Real maxLength);
    ngl::Vec4 faceforward(const ngl::Vec4& n, const ngl::Vec4& v);
    ngl::Vec4 reflect(const ngl::Vec4& v, const ngl::Vec4& n);
}

#endif // UTILS_H
