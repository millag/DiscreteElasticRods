#pragma once

#include "ElasticRod.h"
#include "RenderObject.h"

class Exporter
{
public:
	bool exportGeometry(const char *filename, const RenderObject& object);
	bool exportCurves(const char *filename, const std::vector<ElasticRod*>& strands);
};
