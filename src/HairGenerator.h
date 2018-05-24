#ifndef HAIRGENERATOR_H
#define HAIRGENERATOR_H

#include "Hair.h"

class HairGenerator
{
public:
    HairGenerator();
    ~HairGenerator();

    static void generateCurlyHair(const RenderObject* object, const std::vector<unsigned>& findices, Hair& o_hair);
    static void generateStraightHair(const RenderObject* object, const std::vector<unsigned>& findices, Hair& o_hair);

    static void generateHelicalRod(const HairParams& params,
                            const mg::Vec3D& p, const mg::Vec3D& n, const mg::Vec3D& u,
                            ElasticRod& o_rod);
    static void generateStraightRod(const HairParams& params,
                             const mg::Vec3D& p, const mg::Vec3D& n, const mg::Vec3D& u,
                             ElasticRod& o_rod);
};

#endif // HAIRGENERATOR_H
