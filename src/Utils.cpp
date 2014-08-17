#include "Utils.h"
#include <cstdlib>
#include <cassert>

mg::Real mg::randf(mg::Real min, mg::Real max) {
    return min + (max - min) * ((mg::Real)std::rand() / RAND_MAX);
}

mg::Vec3D mg::genRandPointInBox(mg::Real bBoxMin, mg::Real bBoxMax )
{
    return mg::Vec3D(randf(bBoxMin, bBoxMax), randf(bBoxMin, bBoxMax), randf(bBoxMin, bBoxMax));
}

// generate with uniform distribution - thank you http://mathworld.wolfram.com
mg::Vec3D mg::genRandPointOnSphere(mg::Real radius, const mg::Vec3D& center)
{
    mg::Real u = randf(-1, 1);
    mg::Real theta = randf(0, mg::Constants::two_pi());

    mg::Real x = std::sqrt(1 - u * u) * std::cos(theta) * radius;
    mg::Real y = std::sqrt(1 - u * u) * std::sin(theta) * radius;
    mg::Real z = u * radius;
    return mg::Vec3D(x, y, z) + center;
}

// generate with uniform distribution - thank you http://mathworld.wolfram.com
mg::Vec3D mg::genRandPointOnDisk(mg::Real radius, const mg::Vec3D& center)
{
    mg::Real theta = randf() * mg::Constants::two_pi();
    mg::Real r = std::sqrt(randf()) * radius;

    mg::Real x = std::cos(theta) * r;
    mg::Real y = std::sin(theta) * r;
    mg::Real z = 0;

    return  mg::Vec3D(x, y, z) + center;
}


mg::Real mg::getSign(mg::Real value)
{
    return (value < 0)? -1 : 1;
}

void mg::truncate(mg::Vec3D& io_v, mg::Real maxLength)
{
    mg::Real maxLengthSqr = maxLength * maxLength;
    if (maxLengthSqr > ERR  && io_v.length_squared() > maxLengthSqr)
    {
        io_v.normalize();
        io_v *= maxLength;
    }
}

mg::Vec3D mg::faceforward(const mg::Vec3D& n, const mg::Vec3D& v)
{
    mg::Vec3D normal = n;
    return (mg::dot(v, normal) > 0)? -normal : normal;
}

mg::Vec3D mg::reflect(const mg::Vec3D& v, const mg::Vec3D& n)
{
    assert(n.length_squared() > ERR);

    mg::Vec3D normal = faceforward(n, v);
    normal.normalize();

    return v - normal * 2 * mg::dot(v, normal);
}
