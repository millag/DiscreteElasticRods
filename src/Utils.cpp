#include "Utils.h"
#include <cstdlib>
#include "ngl/Util.h"


mg::Real mg::randf(mg::Real min, mg::Real max) {
    return min + (max - min) * ((mg::Real)std::rand() / RAND_MAX);
}

ngl::Vec4 mg::genRandPointInBox(mg::Real bBoxMin, mg::Real bBoxMax )
{
    return ngl::Vec4(randf(bBoxMin, bBoxMax), randf(bBoxMin, bBoxMax), randf(bBoxMin, bBoxMax));
}

// generate with uniform distribution - thank you http://mathworld.wolfram.com
ngl::Vec4 mg::genRandPointOnSphere(mg::Real radius, const ngl::Vec4& center)
{
    mg::Real u = randf(-1, 1);
    mg::Real theta = randf(0, ngl::PI * 2);

    mg::Real x = std::sqrt(1 - u * u) * std::cos(theta) * radius;
    mg::Real y = std::sqrt(1 - u * u) * std::sin(theta) * radius;
    mg::Real z = u * radius;
    return ngl::Vec4(x, y, z) + center;
}

// generate with uniform distribution - thank you http://mathworld.wolfram.com
ngl::Vec4 mg::genRandPointOnDisk(mg::Real radius, const ngl::Vec4& center)
{
    mg::Real theta = randf() * ngl::PI * 2;
    mg::Real r = std::sqrt(randf()) * radius;

    mg::Real x = std::cos(theta) * r;
    mg::Real y = std::sin(theta) * r;
    mg::Real z = 0;

    return  ngl::Vec4(x, y, z) + center;
}


mg::Real mg::getSign(mg::Real value)
{
    return (value < 0)? -1 : 1;
}

void mg::truncate(ngl::Vec4& io_v, mg::Real maxLength)
{
    mg::Real maxLengthSqr = maxLength * maxLength;
    if (maxLengthSqr > ERR  && io_v.lengthSquared() > maxLengthSqr)
    {
        io_v.normalize();
        io_v *= maxLength;
    }
}

ngl::Vec4 mg::faceforward(const ngl::Vec4& n, const ngl::Vec4& v)
{
    ngl::Vec4 normal = n;
    return (v.dot(normal) > 0)? normal.negate() : normal;
}

ngl::Vec4 mg::reflect(const ngl::Vec4& v, const ngl::Vec4& n)
{
    assert(n.lengthSquared() > ERR);

    ngl::Vec4 normal = faceforward(n, v);
    normal.normalize();

    return v - normal * 2 * v.dot(normal);
}
