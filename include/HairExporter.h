#ifndef HAIREXPORTER_H
#define HAIREXPORTER_H

#include "ElasticRod.h"

class HairExporter
{
public:
    HairExporter() { }
    ~HairExporter() { }

    void exportHair(const char *filename, const std::vector<ElasticRod*>& strands);
};

#endif // HAIREXPORTER_H
