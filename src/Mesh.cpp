#include "Mesh.h"
#include <limits>


unsigned Mesh::getNPrimitives() const
{
    switch (m_mode)
    {
    case PrimitiveMode::TRIANGLES:
        return (m_vindices.size() / 3);

    case PrimitiveMode::LINES:
        return (m_vindices.size() / 2);

    default:
        break;
    }
    return 0;
}

unsigned Mesh::getPrimitiveOffset(unsigned primitiveIdx) const
{
    assert( ((primitiveIdx + 1) * getNVerticesPerPrimitive()) <= m_vindices.size() );
    return primitiveIdx * getNVerticesPerPrimitive();
}

//TODO: mode is always ignored - GOURAUD assumed
void Mesh::computeNormals(Mesh& o_mesh, ShadingMode::Enum mode)
{
    o_mesh.m_normals.resize(o_mesh.m_vertices.size(), mg::Vec3D(0,0,0));

    mg::Vec3D dirx, dirz;
    for (unsigned j = 0; j < o_mesh.m_vindices.size(); j += 3)
    {
        for (unsigned i = 0; i < 3; ++i)
        {
            dirz = o_mesh.m_vertices[ o_mesh.m_vindices[j + (i + 1)%3] ] -  o_mesh.m_vertices[ o_mesh.m_vindices[j + i] ];
            dirx = o_mesh.m_vertices[ o_mesh.m_vindices[j + (i + 2)%3] ] -  o_mesh.m_vertices[ o_mesh.m_vindices[j + i] ];
            o_mesh.m_normals[ o_mesh.m_vindices[j + i] ] += mg::cross(dirz, dirx);
        }
    }

    typedef std::vector<mg::Vec3D>::iterator Iter;
    for (Iter it = o_mesh.m_normals.begin(); it != o_mesh.m_normals.end(); ++it)
    {
        it->normalize();
    }
}


Mesh* Mesh::createSphere(int id, unsigned divu, unsigned divv)
{
    Mesh* mesh = new Mesh(id);

    const mg::Real radius = 1.0;

    mesh->m_vertices.push_back(mg::Vec3D(0,radius,0));
    mg::Real y, r, x, z;
    for (unsigned i = 1; i < divv; i++)
    {
        y = std::cos( ((mg::Real)i / divv) * mg::Constants::pi() );
        r = std::sin( ((mg::Real)i / divv) * mg::Constants::pi() );

        for (unsigned j = 0; j < divu; j++)
        {
            x = std::cos( ((mg::Real)j / divu) * mg::Constants::two_pi()) * r;
            z = std::sin( ((mg::Real)j / divu) * mg::Constants::two_pi()) * r;
            mesh->m_vertices.push_back(mg::Vec3D(x,y,z));
        }
    }
    mesh->m_vertices.push_back(mg::Vec3D(0,-radius,0));

    for (unsigned j = 0; j < divu; j++)
    {
        unsigned idxp1 = 0;
        unsigned idxc1 = j % divu + 1;
        unsigned idxc2 = (j + 1) % divu + 1;

        mesh->m_vindices.push_back(idxp1);
        mesh->m_vindices.push_back(idxc2);
        mesh->m_vindices.push_back(idxc1);

        idxc1 = (divv - 2) * divu + j % divu + 1;
        idxc2 = (divv - 2) * divu + (j + 1) % divu + 1;
        idxp1 = (divv - 1) * divu + 1;

        mesh->m_vindices.push_back(idxp1);
        mesh->m_vindices.push_back(idxc1);
        mesh->m_vindices.push_back(idxc2);

    }

    for (unsigned i = 1; i < divv - 1; i++)
    {
        for (unsigned j = 0; j < divu; j++)
        {
            unsigned idxp1 = (i - 1) * divu + j % divu + 1;
            unsigned idxp2 = (i - 1) * divu + (j + 1) % divu + 1;
            unsigned idxc1 = i * divu + j % divu + 1;
            unsigned idxc2 = i * divu + (j + 1) % divu + 1;

            mesh->m_vindices.push_back(idxp1);
            mesh->m_vindices.push_back(idxp2);
            mesh->m_vindices.push_back(idxc1);

            mesh->m_vindices.push_back(idxc2);
            mesh->m_vindices.push_back(idxc1);
            mesh->m_vindices.push_back(idxp2);
        }
    }

    Mesh::computeNormals( *mesh,  Mesh::ShadingMode::GOURAUD );
    return mesh;
}
