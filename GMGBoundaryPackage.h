#ifndef GMGBOUNDARYPACKAGE_H
#define GMGBOUNDARYPACKAGE_H
#pragma once
#include "FirstOrderKFVSBoundaryPackage.h"
#include "GMGFieldPackage.h"
#include "HighOrderGKSBoundaryPackage.h"

template<typename BoundaryOrderType>
class GMGBoundaryPackage
{
public:
	std::vector<BoundaryOrderType> boundaryPacks;
	//Revise Begin if level==0&&order==3,add third order boundary manager
	//note only have dp_field not avg_slope.
	//we should construct it use 0 level's mesh pointer,0 level's dp_fields,0 level's highorder package
	HighOrderGKSBoundaryPackage<5> thirdorderPacks;
public:
	GMGFieldPackage& fieldpacks;
public:
	GMGBoundaryPackage(GMGFieldPackage& packs)
		:
		fieldpacks(packs),thirdorderPacks(packs.meshpack.grids[0],packs.packages[0],packs.pack)
		//Revise End
	{
		for (std::size_t ilevel = 0; ilevel < packs.meshpack.nMGLevel; ilevel++)
		{
			SubFvZone* mesh = fieldpacks.meshpack.grids[ilevel];
			DataPackage& pack = fieldpacks.packages[ilevel];
			this->boundaryPacks.push_back(BoundaryOrderType(mesh, pack));
		}
	}
public:
	void SetSuperSonic
	(
		const std::string& BoundaryName,
		double Rho,
		const Vector& U,
		double p,
		double gamma,
		int K
	)
	{
		for (std::size_t ilevel = 0; ilevel < this->boundaryPacks.size(); ilevel++)
		{
			this->boundaryPacks[ilevel].SetSuperSonic(BoundaryName, Rho, U, p, gamma, K);
		}
		//Revise Begin Add ThirdOrderBoundary,HighOrderPackage,Already have
		this->thirdorderPacks.SetSuperSonic(BoundaryName, Rho, U, p, gamma, K);
		//Revise End
	}

	void SetAdiabaticNonSlipWall(const std::string& BoundaryName)
	{
		for (std::size_t ilevel = 0; ilevel < this->boundaryPacks.size(); ilevel++)
		{
			this->boundaryPacks[ilevel].SetAidabaticNonSlipWall(BoundaryName);
		}
		//Revise Begin Add ThirdOrderBoundary
		this->thirdorderPacks.SetAidabaticNonSlipWall(BoundaryName);
		//Revise End
	}

	void SetEulerSlipWall(const std::string& BoundaryName)
	{
		for (std::size_t ilevel = 0; ilevel < this->boundaryPacks.size(); ilevel++)
		{
			this->boundaryPacks[ilevel].SetEulerSlipWall(BoundaryName);
		}
		//Revise Begin,Add third Order Boundary
		this->thirdorderPacks.SetEulerSlipWall(BoundaryName);
		//Revise End
	}

	void SetSubSonic
	(
		const std::string& BoundaryName,
		double Rho,
		const Vector& U,
		double p,
		double gamma,
		int K
	)
	{
		for (std::size_t ilevel = 0; ilevel < this->boundaryPacks.size(); ilevel++)
		{
			this->boundaryPacks[ilevel].SetSubSonic(BoundaryName, Rho, U, p, gamma, K);
		}
		//Revise Begin,Add third order Boundary
		this->thirdorderPacks.SetSubSonic(BoundaryName, Rho, U, p, gamma, K);
		//Revise End
	}

	void SetFree
	(
		const std::string& BoundaryName
	)
	{
		for (std::size_t ilevel = 0; ilevel < this->boundaryPacks.size(); ilevel++)
		{
			this->boundaryPacks[ilevel].SetFree(BoundaryName);
		}
		//Revise Begin,Add third order Boundary
		this->thirdorderPacks.SetFree(BoundaryName);
		//Revise End
	}

};


#endif;