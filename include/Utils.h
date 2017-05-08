#ifndef UTILS_H
#define UTILS_H

#include "config.h"

namespace mg
{
    const static unsigned short SEC = 1000;
    const static unsigned short UPS = 30;
    const static unsigned short FPS = 30;
    const static unsigned short UPDATERATE = SEC / UPS;
    const static unsigned short REFRESHRATE = SEC / FPS;

    const static Real ERR = 1e-6f;
    const static Real THRESHODL = 0.00001f;
    const static mg::Vec3D EX(1, 0, 0);
    const static mg::Vec3D EY(0, 1, 0);
    const static mg::Vec3D EZ(0, 0, 1);

    const static mg::Vec3D GRAVITY(0, -9.81f, 0);

    mg::Real getSign(mg::Real value);
    mg::Real randf(mg::Real min = 0.0, mg::Real max = 1.0);
    mg::Vec3D genRandPointInBox(mg::Real bBoxMin = -1.0, mg::Real bBoxMax = 1.0);
    mg::Vec3D genRandPointOnSphere(mg::Real radius = 1.0, const mg::Vec3D& center = mg::Vec3D());
    mg::Vec3D genRandPointOnDisk(mg::Real radius = 1.0, const mg::Vec3D& center = mg::Vec3D());

    void truncate(mg::Vec3D& io_v, mg::Real maxLength);
    mg::Vec3D faceforward(const mg::Vec3D& n, const mg::Vec3D& v);
    mg::Vec3D reflect(const mg::Vec3D& v, const mg::Vec3D& n);
}

#endif // UTILS_H
