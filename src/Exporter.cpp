#include "Exporter.h"

#include <iostream>
#include <fstream>


bool Exporter::exportGeometry(const char *filename, const RenderObject& object)
{
	std::ofstream ofile (filename);
	if (!ofile.is_open())
	{
		std::cerr << "Can't open file " + std::string(filename) << std::endl;
		ofile.close();
		return false;
	}

	ofile << "# Total vertices " << object.getMesh()->m_vertices.size() << "\n";
	ofile << "# Total normals " << object.getMesh()->m_normals.size() << "\n";
	ofile << "# Total faces " << object.getMesh()->getNPrimitives() << "\n";

	mg::Vec3D v;
	mg::Matrix4D transform = object.getTransform();
	mg::Matrix4D transformInv = mg::inverse(object.getTransform()).transpose();

	typedef std::vector<mg::Vec3D>::const_iterator Iter;
	for (Iter it = object.getMesh()->m_vertices.begin(); it != object.getMesh()->m_vertices.end(); ++it)
	{
		v = mg::transform_point(transform, *it);
		ofile << "v " << v[0] << " " <<  v[1] << " " << v[2] << "\n";
	}

	for (Iter it = object.getMesh()->m_normals.begin(); it != object.getMesh()->m_normals.end(); ++it)
	{
		v = mg::transform_vector(transformInv, *it);
		ofile << "vn " << v[0] << " " <<  v[1] << " " << v[2] << "\n";
	}

	for (unsigned i = 0; i < object.getMesh()->m_vindices.size(); i+= object.getMesh()->getNVerticesPerPrimitive())
	{
		ofile << "f";
		for (unsigned j = 0; j < object.getMesh()->getNVerticesPerPrimitive(); ++j)
		{
			ofile << " " << (object.getMesh()->m_vindices[i + j] + 1) << "//" <<  (object.getMesh()->m_vindices[i + j] + 1);
		}
		ofile << "\n";
	}

	ofile.flush();
	// close on complete
	ofile.close();
	return true;
}

bool Exporter::exportCurves(const char *filename, const std::vector<ElasticRod *> &strands)
{
	std::ofstream ofile (filename);
	if (!ofile.is_open())
	{
		std::cerr << "Can't open file " + std::string(filename) << std::endl;
		ofile.close();
		return false;
	}

	unsigned long nv = 0;
	unsigned long nvBegin = 0;
	unsigned long nvEnd = 0;
	unsigned long nCurve = 0;

	typedef std::vector<ElasticRod*>::const_iterator Iter;
	typedef std::vector<mg::Vec3D>::const_iterator VIter;
	for (Iter it = strands.begin(); it != strands.end(); ++it)
	{
		ofile << "g curve" << nCurve << "\n";
		nvBegin = nv;
		for (VIter vit = (*it)->m_ppos.begin(); vit != (*it)->m_ppos.end(); ++vit)
		{
			ofile << "v " << (*vit)[0] << " " <<  (*vit)[1] << " " << (*vit)[2] << "\n";
			++nv;
		}
		nvEnd = nv;
		ofile << "g curve" << nCurve << "\n";
		nCurve++;

		if ((nvEnd - nvBegin) == 0)
		{
			continue;
		}
		ofile << "l";
//        for (unsigned long i = nvBegin + 1; i < nvEnd + 1; ++i)
//        {
//            ofile << " " << i;
//        }
		for (unsigned long i = nvEnd; i > nvBegin; --i)
		{
			ofile << " " << i;
		}
		ofile << "\n";
	}

	ofile << "# Total vertices " << nv << "\n";
	ofile << "# Total curves " << strands.size() << "\n";
	ofile.flush();
	// close on complete
	ofile.close();
	return true;
}
