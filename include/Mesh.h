#ifndef MESH_H
#define MESH_H

#include <vector>
#include <cassert>
#include "config.h"

class Mesh
{
public:
    struct ShadingMode{ enum Enum { FLAT, GOURAUD }; };
    struct PrimitiveMode{ enum Enum { TRIANGLES = 3, LINES = 2 }; };

    Mesh(unsigned id = -1, PrimitiveMode::Enum mode = PrimitiveMode::TRIANGLES):m_id(id), m_mode(mode) { }
    ~Mesh() { }

    inline unsigned getId() const { return m_id; }
    inline void setId(unsigned id) { m_id = id; }

    inline unsigned getNVertices() const { return (m_vertices.size()); }
    inline unsigned getNVerticesPerPrimitive() const { return m_mode; }
    unsigned getPrimitiveOffset(unsigned primitiveIdx) const;
    unsigned getNPrimitives() const;

public:
    static Mesh* createSphere(int id, unsigned divu = 20, unsigned divv = 10);
    static void computeNormals(Mesh& o_mesh, ShadingMode::Enum mode);

public:
    std::vector<mg::Vec3D> m_vertices;
    std::vector<mg::Vec3D> m_normals;
    std::vector<unsigned> m_vindices;

private:
    unsigned m_id;
    PrimitiveMode::Enum m_mode;
};

#endif // MESH_H
