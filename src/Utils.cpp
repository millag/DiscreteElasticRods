#include "Utils.h"
#include <cstdlib>
#include "ngl/Util.h"


ngl::Real utils::randf(ngl::Real min, ngl::Real max) {
    return min + (max - min) * ((ngl::Real)std::rand() / RAND_MAX);
}

ngl::Vec4 utils::genRandPointInBox(ngl::Real bBoxMin, ngl::Real bBoxMax )
{
    return ngl::Vec4(randf(bBoxMin, bBoxMax), randf(bBoxMin, bBoxMax), randf(bBoxMin, bBoxMax));
}

// generate with uniform distribution - thank you http://mathworld.wolfram.com
ngl::Vec4 utils::genRandPointOnSphere(ngl::Real radius, const ngl::Vec4& center)
{
    ngl::Real u = randf(-1, 1);
    ngl::Real theta = randf(0, ngl::PI * 2);

    ngl::Real x = std::sqrt(1 - u * u) * std::cos(theta) * radius;
    ngl::Real y = std::sqrt(1 - u * u) * std::sin(theta) * radius;
    ngl::Real z = u * radius;
    return ngl::Vec4(x, y, z) + center;
}

// generate with uniform distribution - thank you http://mathworld.wolfram.com
ngl::Vec4 utils::genRandPointOnDisk(ngl::Real radius, const ngl::Vec4& center)
{
    ngl::Real theta = randf() * ngl::PI * 2;
    ngl::Real r = std::sqrt(randf()) * radius;

    ngl::Real x = std::cos(theta) * r;
    ngl::Real y = std::sin(theta) * r;
    ngl::Real z = 0;

    return  ngl::Vec4(x, y, z) + center;
}


ngl::Real utils::getSign(ngl::Real value)
{
    return (value < 0)? -1 : 1;
}

void utils::truncate(ngl::Vec4& io_v, ngl::Real maxLength)
{
    ngl::Real maxLengthSqr = maxLength * maxLength;
    if (maxLengthSqr > ERR  && io_v.lengthSquared() > maxLengthSqr)
    {
        io_v.normalize();
        io_v *= maxLength;
    }
}

ngl::Vec4 utils::faceforward(const ngl::Vec4& n, const ngl::Vec4& v)
{
    ngl::Vec4 normal = n;
    return (v.dot(normal) > 0)? normal.negate() : normal;
}

ngl::Vec4 utils::reflect(const ngl::Vec4& v, const ngl::Vec4& n)
{
    assert(n.lengthSquared() > ERR);

    ngl::Vec4 normal = faceforward(n, v);
    normal.normalize();

    return v - normal * 2 * v.dot(normal);
}
