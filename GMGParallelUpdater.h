#ifndef GMGPARALLELUPDATER_H
#define GMGPARALLELUPDATER_H
#pragma once
#include "FieldUpdate.h"
#include "GMGFieldPackage.h"
#include "ParallelUpdate.h"

class GMGParallelUpdater
{
public:
	std::vector<FieldUpdate> updaters;
	GMGFieldPackage& packs;
public:
	GMGParallelUpdater(GMGFieldPackage& given);
public:
	void Update(int level,int order=1);
public:
	//Revise Begin declare hihgorder updater
	//note after this first level will update in hipdater;
	ParallelUpdate<HighOrderDistribution> hupdater;
	//Revise end
public:
	void UpdateHighOrderDistribution();
};

#endif;