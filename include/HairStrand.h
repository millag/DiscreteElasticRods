#ifndef HAIRSTRAND_H
#define HAIRSTRAND_H

#include "Hair.h"

class HairGenerator
{
public:
    HairGenerator();
    ~HairGenerator();

    void generateCurlyHair(const RenderObject* object, Hair& o_hair) const;
    void generateStraightHair(const RenderObject* object, Hair& o_hair) const;

    void generateHelicalRod(const HairParams& params,
                            const mg::Vec3D& p, const mg::Vec3D& n, const mg::Vec3D& u,
                            ElasticRod& o_rod) const;
    void generateStraightRod(const HairParams& params,
                             const mg::Vec3D& p, const mg::Vec3D& n, const mg::Vec3D& u,
                             ElasticRod& o_rod) const;

};

#endif // HAIRSTRAND_H
