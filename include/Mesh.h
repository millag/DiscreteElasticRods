#ifndef MESH_H
#define MESH_H

#include <vector>
#include "ngl/Vec4.h"

class Mesh
{
public:
    enum SHADING_MODE { FLAT, GOURAUD };
    static void calcNormals(Mesh& o_mesh, Mesh::SHADING_MODE mode);


    Mesh(unsigned id): m_id(id) {}
    virtual ~Mesh() {}

    unsigned getId() const { return m_id; }
    unsigned getNFaces() const { return (m_vindices.size() / 3); }
    unsigned getNVertices() const { return (m_vertices.size()); }

    std::vector<ngl::Vec4> m_vertices;
    std::vector<unsigned> m_vindices;
    std::vector<ngl::Vec4> m_normals;
    std::vector<unsigned> m_nindices;

protected:
    const unsigned m_id;

private:
    Mesh():m_id(-1) {}
};

#endif // MESH_H
