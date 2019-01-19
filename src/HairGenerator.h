#pragma once

#include "Hair.h"

class HairGenerator
{
public:
	static void generateCurlyHair(const RenderObject* object, const std::vector<unsigned>& findices, Hair& o_hair);
	static void generateStraightHair(const RenderObject* object, const std::vector<unsigned>& findices, Hair& o_hair);

	static void generateHelicalRod(const HairParams& params,
	                        const mg::Vec3D& p, const mg::Vec3D& n, const mg::Vec3D& u,
	                        ElasticRod& o_rod);
	static void generateStraightRod(const HairParams& params,
	                         const mg::Vec3D& p, const mg::Vec3D& n, const mg::Vec3D& u,
	                         ElasticRod& o_rod);
};
