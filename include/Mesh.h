#ifndef MESH_H
#define MESH_H

#include "config.h"

#include <vector>

class Mesh
{
public:
	enum ShadingMode
	{
		FLAT,
		GOURAUD,
	};

	enum Primitive
	{
		POINTS = 1,
		LINES = 2,
		TRIANGLES = 3,
	};

	static Mesh* createSphere( int id, unsigned udiv = 20u, unsigned vdiv = 10u );
	static void computeNormals( ShadingMode mode, Mesh& o_mesh );

	Mesh( unsigned id = -1, Primitive primitive = TRIANGLES );
	~Mesh();

	unsigned getId() const { return m_id; }
	void setId( unsigned id ) { m_id = id; }

	Primitive getPrimitiveType() const { return m_primitive; }
	unsigned getNPrimitives() const;
	unsigned getPrimitiveOffset( unsigned primitiveIdx ) const;
	unsigned getNVertices() const { return static_cast<unsigned>( m_vertices.size() ); }
	unsigned getNVerticesPerPrimitive() const { return static_cast<unsigned>( m_primitive ); }

	bool hasNormals() const { return m_normals.size() > 0; }

public:
	std::vector<mg::Vec3D> m_vertices;
	std::vector<mg::Vec3D> m_normals;
	std::vector<unsigned> m_vindices;

private:
	unsigned m_id = -1;
	Primitive m_primitive = TRIANGLES;
};

#endif // MESH_H
