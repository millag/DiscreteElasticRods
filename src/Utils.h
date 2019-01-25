#pragma once

#include "config.h"

namespace mg
{

template <typename T, int N>
constexpr int CountOf(T const (&)[N]) noexcept
{
	return N;
}

constexpr inline mg::Real getSign(mg::Real value)
{
	return static_cast<mg::Real>( (value < 0)? -1 : 1 );
}

inline mg::Real randf(mg::Real min = 0.0, mg::Real max = 1.0)
{
	return min + (max - min) * ((mg::Real)std::rand() / RAND_MAX);
}

inline mg::Vec3D genRandPointInBox(mg::Real bBoxMin = -1.0, mg::Real bBoxMax = 1.0)
{
	return mg::Vec3D(randf(bBoxMin, bBoxMax), randf(bBoxMin, bBoxMax), randf(bBoxMin, bBoxMax));
}

/// Generate random point on sphere with uniform distribution
mg::Vec3D genRandPointOnSphere(mg::Real radius = 1.0, const mg::Vec3D& center = mg::Vec3D());

/// Generate random point on disk with uniform distribution
mg::Vec3D genRandPointOnDisk(mg::Real radius = 1.0, const mg::Vec3D& center = mg::Vec3D());

inline void truncate(mg::Vec3D& io_v, mg::Real maxLength)
{
	if ( io_v.length_squared() > maxLength * maxLength )
	{
		io_v.normalize();
		io_v *= maxLength;
	}
}

inline mg::Vec3D faceforward( const mg::Vec3D& v, const mg::Vec3D& n )
{
	return ( mg::dot( v, n ) > 0 )? -n : n;
}

inline mg::Vec3D reflect( const mg::Vec3D& v, const mg::Vec3D& n )
{
	mg::Vec3D normal = faceforward( v, n );
	normal.normalize();
	return v - normal * 2 * mg::dot( v, normal );
}

template <typename V, typename T>
inline V lerp( const V& v1, const V& v2, const T& t )
{
	return  v1 * (1 - t) + v2 * t;
}

inline bool almostEqual( mg::Real v1, mg::Real v2, mg::Real tolerance = mg::ERR )
{
	return std::fabs( v1 - v2 ) < tolerance;
}
}
