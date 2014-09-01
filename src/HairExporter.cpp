#include "HairExporter.h"

#include <iostream>
#include <fstream>

void HairExporter::exportHair(const char *filename, const std::vector<ElasticRod *> &strands)
{
    std::ofstream ofile (filename);

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
        for (unsigned long i = nvBegin + 1; i < nvEnd + 1; ++i)
        {
            ofile << " " << i;
        }
//        for (unsigned long i = nvEnd; i > nvBegin; --i)
//        {
//            ofile << " " << i;
//        }
        ofile << "\n";
    }

    ofile << "# Total vertices " << nv << "\n";
    ofile << "# Total curves " << strands.size() << "\n";

    // close on complete
    ofile.close();

}
