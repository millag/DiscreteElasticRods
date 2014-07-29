#include "Mesh.h"

#include <limits>
#include <ngl/Util.h>

//TODO: mode is ignored for now - FLAT assumed
void Mesh::calcNormals(Mesh& o_mesh, Mesh::SHADING_MODE mode)
{
    o_mesh.m_normals.resize(o_mesh.m_vindices.size());

    unsigned i  = 0;
    typedef std::vector<unsigned>::const_iterator ULIter;
    for (ULIter it = o_mesh.m_vindices.begin(); it != o_mesh.m_vindices.end(); it += 3)
    {
        o_mesh.m_normals[3*i] = ngl::calcNormal( o_mesh.m_vertices[ *(it) ], o_mesh.m_vertices[ *(it + 2) ], o_mesh.m_vertices[ *(it + 1) ] );
        o_mesh.m_normals[3*i + 1] = o_mesh.m_normals[3*i];
        o_mesh.m_normals[3*i + 2] = o_mesh.m_normals[3*i];
        ++i;
    }
    o_mesh.m_nindices.resize(o_mesh.m_vindices.size());
    o_mesh.m_nindices.insert(o_mesh.m_nindices.begin(), o_mesh.m_vindices.begin(), o_mesh.m_vindices.end());
}
