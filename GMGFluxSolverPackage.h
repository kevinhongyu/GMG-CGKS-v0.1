#ifndef GMGFLUXSOLVERPACKAGE_H
#define GMGFLUXSOLVERPACKAGE_H
#pragma once
#include "FirstOrderUpdate.h"
#include "GMGFieldPackage.h"
#include "GMGMeshPackage.h"
#include "GMGBoundaryTrait.h"
#include "GMGBoundaryPackage.h"
#include "HighOrderUpdate.h"

template <template<typename TauType,int N=1> typename FluxType,int TotalStage=1>
class GMGFluxSolverPackage
{
public:
	using BoundaryType = typename FluxSolverToBoundary<FluxType<NS, TotalStage>>::BoundaryType;
public:
	std::vector<FluxType<Euler,TotalStage>> fluxSolvers;
	//Revise Begin,Add ThirdOrder Flux Solver
	HighOrderUpdate<Euler, TotalStage> thirdorderfluxSolver;
public:
	GMGFluxSolverPackage(GMGFieldPackage& fields, GMGBoundaryPackage<BoundaryType>& Bpackage)
		:
		thirdorderfluxSolver(fields.packages[0],fields.avg_slope,fields.pack,fields.cf,fields.meshpack.grids[0],Bpackage.thirdorderPacks)
		//Revise End
	{
		for (int ilevel = 0; ilevel < fields.meshpack.nMGLevel; ilevel++)
		{
			this->fluxSolvers.push_back(FluxType<Euler, TotalStage>(
				fields.packages[ilevel],
				fields.meshpack.grids[ilevel],
				Bpackage.boundaryPacks[ilevel]
				));
		}
	}
public:
	/*FluxSolver.fluxsolver.Mu = Miu;
	FluxSolver.fluxsolver.r = gammac;
	FluxSolver.fluxsolver.R_gas = cas.Rgas;
	FluxSolver.fluxsolver.C2 = cas.C2;
	int K = cas.K;*/
	void SetSolverParameter(
		double Miu, double gammac, double R_gas, double C2, int K_
	,int Stage)
	{
		for (std::size_t ilevel = 0; ilevel < this->fluxSolvers.size(); ilevel++)
		{
			this->fluxSolvers[ilevel].fluxsolver.Mu = Miu;
			this->fluxSolvers[ilevel].fluxsolver.r = gammac;
			this->fluxSolvers[ilevel].fluxsolver.R_gas = R_gas;
			this->fluxSolvers[ilevel].fluxsolver.C2 = C2;
			this->fluxSolvers[ilevel].block.SetBlock3D(Stage);
		}
		//Revise Begin ,Add Third Order Parameter
		this->thirdorderfluxSolver.fluxsolver.Mu = Miu;
		this->thirdorderfluxSolver.fluxsolver.r = gammac;
		this->thirdorderfluxSolver.fluxsolver.R_gas = R_gas;
		this->thirdorderfluxSolver.fluxsolver.C2 = C2;
		this->thirdorderfluxSolver.block.SetBlock3D(Stage);
		//Revise End

	}
};

#endif;