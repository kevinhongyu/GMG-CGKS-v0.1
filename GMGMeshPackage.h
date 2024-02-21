#ifndef GMGMESHPACKAGE_H
#define GMGMESHPACKAGE_H
#pragma once
#include "SubFvZone.h"

class GMGMeshPackage
{
public:
	std::vector<SubFvZone*>& grids;
public:
	int nMGLevel;
public:
	GMGMeshPackage(std::vector<SubFvZone*>& mesh);
};

#endif;