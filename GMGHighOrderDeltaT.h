#ifndef GMGHIGHORDERDELTAT_H
#define GMGHIGHORDERDELTAT_H
#pragma once
#include "HighOrderTransient.h"
#include "GMGFieldPackage.h"
#include "Block3d.h"

class GMGHighOrderDeltaT
{
public:
	std::vector<HighOrderTransient> DTs;
	GMGFieldPackage& fields;
public:
	GMGHighOrderDeltaT(
		double R_,
		double r_,
		double CFL,
		double Mu_,
		int K_, 
		GMGFieldPackage& fi
	);
};

#endif;