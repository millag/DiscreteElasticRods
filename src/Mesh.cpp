#include "Mesh.h"

Mesh::Mesh( unsigned id, Primitive primitive ):
    m_id( id ),
    m_primitive( primitive )
{ }

Mesh::~Mesh() = default;

unsigned Mesh::getNPrimitives() const
{
	switch ( m_primitive )
	{
	case Primitive::TRIANGLES:
	case Primitive::LINES:
		return static_cast<unsigned>( m_vindices.size() / getNVerticesPerPrimitive() );
	case Primitive::POINTS:
		return static_cast<unsigned>( m_vertices.size() );
	default:
		break;
	}

	return 0;
}

unsigned Mesh::getPrimitiveOffset(unsigned primitiveIdx) const
{
	assert( ( ( primitiveIdx + 1 ) * getNVerticesPerPrimitive() ) <= m_vindices.size() );
	return primitiveIdx * getNVerticesPerPrimitive();
}

Mesh* Mesh::createSphere( int id, unsigned udiv, unsigned vdiv )
{
	Mesh* mesh = new Mesh( id );

//	north pole vertex
	mesh->m_vertices.emplace_back( 0.f, 1.f, 0.f );

	mg::Real y, r, x, z;
	for ( auto i = 1u; i < vdiv; ++i )
	{
		const auto vangle = ( static_cast<mg::Real>( i ) / vdiv ) * mg::Constants::pi();
		y = std::cos( vangle );
		r = std::sin( vangle );

		for ( auto j = 0u; j < udiv; ++j )
		{
			const auto uangle = ( static_cast<mg::Real>( j ) / udiv ) * mg::Constants::two_pi();
			x = std::cos( uangle ) * r;
			z = std::sin( uangle ) * r;
			mesh->m_vertices.emplace_back( x, y, z );
		}
	}

//	south pole vertex
	mesh->m_vertices.emplace_back( 0.f, -1.f, 0.f );

//	north pole indices
	for ( auto i = 0u; i < udiv; ++i )
	{
		mesh->m_vindices.push_back( 0 );
		mesh->m_vindices.push_back( 1 + ( i + 1 ) % udiv );
		mesh->m_vindices.push_back( 1 + i % udiv );
	}

	for ( auto i = 0u; i < vdiv - 2; ++i )
	{
		const auto topIdx = 1 + i * udiv;
		const auto bottomIdx = 1 + ( i + 1 ) * udiv;
		for ( auto j = 0u; j < udiv; ++j )
		{
			const auto currOffset = j % udiv;
			const auto nextOffset = ( j + 1 ) % udiv;

			mesh->m_vindices.push_back( topIdx + currOffset );
			mesh->m_vindices.push_back( topIdx + nextOffset );
			mesh->m_vindices.push_back( bottomIdx + currOffset );

			mesh->m_vindices.push_back( bottomIdx  + currOffset );
			mesh->m_vindices.push_back( topIdx + nextOffset );
			mesh->m_vindices.push_back( bottomIdx + nextOffset );
		}
	}

//	south pole indices
	const auto lastIdx = static_cast<unsigned>( mesh->m_vertices.size() - 1 );
	for ( auto i = 0u; i < udiv; ++i )
	{
		mesh->m_vindices.push_back( lastIdx );
		mesh->m_vindices.push_back( lastIdx - 1 - ( i + 1 ) % udiv );
		mesh->m_vindices.push_back( lastIdx - 1 - i % udiv );
	}

	Mesh::computeNormals( GOURAUD, *mesh );

	return mesh;
}

//TODO: mode is always ignored - GOURAUD assumed
void Mesh::computeNormals( ShadingMode mode, Mesh& o_mesh )
{
	UNUSED_VALUE( mode );

	o_mesh.m_normals.resize( o_mesh.m_vertices.size(), mg::Vec3D( 0.f, 0.f, 0.f ) );

	mg::Vec3D dirx, dirz;
	for ( auto j = 0u; j < o_mesh.m_vindices.size(); j += 3 )
	{
		for ( auto i = 0; i < 3; ++i )
		{
			dirz = o_mesh.m_vertices[ o_mesh.m_vindices[j + ( i + 1 ) % 3] ] -  o_mesh.m_vertices[o_mesh.m_vindices[j + i]];
			dirx = o_mesh.m_vertices[ o_mesh.m_vindices[j + ( i + 2 ) % 3] ] -  o_mesh.m_vertices[o_mesh.m_vindices[j + i]];
			o_mesh.m_normals[o_mesh.m_vindices[j + i]] += mg::cross( dirz, dirx );
		}
	}

	for ( auto& n : o_mesh.m_normals )
	{
		n.normalize();
	}
}
