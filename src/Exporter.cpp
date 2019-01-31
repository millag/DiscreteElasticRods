#include "Exporter.h"

#include <iostream>
#include <fstream>


bool Exporter::exportGeometry( const char *filename, const RenderObject& object )
{
	std::ofstream ofile( filename );
	if ( !ofile.is_open() )
	{
		std::cerr << "Can't open file " << filename << std::endl;
		ofile.close();
		return false;
	}

	auto mesh = object.getMesh();
	assert( mesh != nullptr );

	ofile << "# Total vertices " << object.getMesh()->m_vertices.size() << "\n";
	ofile << "# Total normals " << object.getMesh()->m_normals.size() << "\n";
	ofile << "# Total faces " << object.getMesh()->getNPrimitives() << "\n";

	for ( const auto& v : mesh->m_vertices )
	{
		const auto wv = mg::transform_point( object.getTransform(), v );
		ofile << "v " << wv[0] << " " <<  wv[1] << " " << wv[2] << "\n";
	}

	const auto normalTransform = mg::inverse( object.getTransform() ).transpose();
	for ( const auto& n : mesh->m_normals )
	{
		const auto wn = mg::transform_vector( normalTransform, n );
		ofile << "vn " << wn[0] << " " <<  wn[1] << " " << wn[2] << "\n";
	}

	for ( auto i = 0u; i < object.getMesh()->m_vindices.size(); i += object.getMesh()->getNVerticesPerPrimitive() )
	{
		ofile << "f";
		for ( auto j = 0u; j < object.getMesh()->getNVerticesPerPrimitive(); ++j )
		{
			ofile << " " << ( object.getMesh()->m_vindices[i+j] + 1 ) << "//" << ( object.getMesh()->m_vindices[i+j] + 1 );
		}
		ofile << "\n";
	}

	ofile.flush();
	ofile.close();
	return true;
}

bool Exporter::exportCurves( const char* filename, const std::vector<ElasticRod>& strands )
{
	std::ofstream ofile( filename );
	if ( !ofile.is_open() )
	{
		std::cerr << "Can't open file " << filename << std::endl;
		ofile.close();
		return false;
	}

	auto nv = 0ul;
	auto nvBegin = 0ul;
	auto nvEnd = 0ul;
	auto nCurve = 0ul;

	for ( const auto& strand : strands )
	{
		ofile << "g curve" << nCurve << "\n";
		nvBegin = nv;
		for ( const auto& pos : strand.m_ppos )
		{
			ofile << "v " << pos[0] << " " <<  pos[1] << " " << pos[2] << "\n";
			++nv;
		}
		nvEnd = nv;
		ofile << "g curve" << nCurve << "\n";
		++nCurve;

		if ( ( nvEnd - nvBegin ) == 0 )
		{
			continue;
		}

		ofile << "l";
//        for (unsigned long i = nvBegin + 1; i < nvEnd + 1; ++i)
//        {
//            ofile << " " << i;
//        }

		for ( auto i = nvEnd; i > nvBegin; --i )
		{
			ofile << " " << i;
		}
		ofile << "\n";
	}

	ofile << "# Total vertices " << nv << "\n";
	ofile << "# Total curves " << strands.size() << "\n";
	ofile.flush();
	ofile.close();

	return true;
}
