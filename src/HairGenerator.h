#pragma once

#include "Hair.h"

class HairGenerator
{
public:

	void generateHair( const RenderObject& skin,
	                   const std::vector<unsigned>& findices,
	                   Hair& o_hair );

	void generateCurlyHair( const RenderObject& skin,
	                        const std::vector<unsigned>& findices,
	                        Hair& o_hair );

	void generateHelicalRod( const HairParams& params,
	                         const mg::Vec3D& root,
	                         const mg::Vec3D& n,
	                         const mg::Vec3D& up,
	                         ElasticRod& o_rod );

	void generateStraightRod( const HairParams& params,
	                          const mg::Vec3D& root,
	                          const mg::Vec3D& n,
	                          const mg::Vec3D& up,
	                          ElasticRod& o_rod );
};
