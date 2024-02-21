#ifndef GMGFIELDPACKAGE_H
#define GMGFIELDPACKAGE_H
#pragma once
#include "GMGMeshPackage.h"
#include "DataPackage.h"
#include "HighOrderGKSPackage.h"

class GMGFieldPackage
{
public:
	GMGMeshPackage& meshpack;
public:
	std::vector<DataPackage> packages;
	//Revise Begin For High Order avg_slope
	DataPackage avg_slope;
	//Revise End
	
	//Revise begin for high order distribution
	HighOrderGKSPackage pack;
	//Revise End
public:
	GMGFieldPackage(GMGMeshPackage& mesh);
public:
	std::vector<std::vector<Field<double>>> fields;
	std::vector<ElementField<HighOrderDistribution>> w_dis;
	std::vector<Field<double>> slopes;
	//Revise Begin,Add cf field
	Field<double> cf;
	//Revise End

};

#endif;