#pragma once

#include "config.h"

namespace mg
{

template <typename T, int N>
constexpr int CountOf(T const (&)[N]) noexcept
{
	return N;
}

mg::Real getSign(mg::Real value);

mg::Real randf(mg::Real min = 0.0, mg::Real max = 1.0);

mg::Vec3D genRandPointInBox(mg::Real bBoxMin = -1.0, mg::Real bBoxMax = 1.0);

mg::Vec3D genRandPointOnSphere(mg::Real radius = 1.0, const mg::Vec3D& center = mg::Vec3D());

mg::Vec3D genRandPointOnDisk(mg::Real radius = 1.0, const mg::Vec3D& center = mg::Vec3D());

void truncate(mg::Vec3D& io_v, mg::Real maxLength);

mg::Vec3D faceforward(const mg::Vec3D& n, const mg::Vec3D& v);

mg::Vec3D reflect(const mg::Vec3D& v, const mg::Vec3D& n);

}
