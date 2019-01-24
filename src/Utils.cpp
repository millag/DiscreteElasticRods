#include "Utils.h"

mg::Vec3D mg::genRandPointOnSphere(mg::Real radius, const mg::Vec3D& center)
{
	mg::Real u = randf(-1, 1);
	mg::Real theta = randf(0, mg::Constants::two_pi());

	mg::Real x = std::sqrt(1 - u * u) * std::cos(theta) * radius;
	mg::Real y = std::sqrt(1 - u * u) * std::sin(theta) * radius;
	mg::Real z = u * radius;
	return mg::Vec3D(x, y, z) + center;
}

mg::Vec3D mg::genRandPointOnDisk(mg::Real radius, const mg::Vec3D& center)
{
	mg::Real theta = randf() * mg::Constants::two_pi();
	mg::Real r = std::sqrt(randf()) * radius;

	mg::Real x = std::cos(theta) * r;
	mg::Real y = std::sin(theta) * r;
	mg::Real z = 0;

	return  mg::Vec3D(x, y, z) + center;
}
