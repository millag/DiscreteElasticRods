#ifndef MESH_H
#define MESH_H

#include <vector>
#include <cassert>
#include "Types.h"

class Mesh
{
public:
    struct ShadingMode{ enum Enum { FLAT, GOURAUD }; };
    struct PrimitiveMode{ enum Enum { TRIANGLES = 3, LINES = 2 }; };

    Mesh(unsigned id, PrimitiveMode::Enum mode = PrimitiveMode::TRIANGLES);
    ~Mesh();

    inline unsigned getId() const { return m_id; }
    inline unsigned getNVertices() const { return (m_vertices.size()); }
    inline unsigned getNVerticesPerPrimitive() const { return m_mode; }
    unsigned getPrimitiveOffset(unsigned primitiveIdx) const;
    unsigned getNPrimitives() const;

public:
    static Mesh* createSphere(int id, unsigned divu = 20, unsigned divv = 10);
    static void calcNormals(Mesh& o_mesh, ShadingMode::Enum mode);

public:
    std::vector<mg::Vec3D> m_vertices;
    std::vector<mg::Vec3D> m_normals;
    std::vector<unsigned> m_vindices;

private:
    Mesh():m_id(-1) {}

private:
    unsigned m_id;
    PrimitiveMode::Enum m_mode;
};

#endif // MESH_H
