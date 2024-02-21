#include "GMGParallelUpdater.h"

GMGParallelUpdater::GMGParallelUpdater(GMGFieldPackage& given)
	:
	packs(given)
{
	this->updaters.resize(packs.packages.size());
	for (std::size_t ilevel = 0; ilevel < packs.meshpack.nMGLevel; ilevel++)
	{
		for (int ivar = 0; ivar < 5; ivar++)
		{
			updaters[ilevel] + (&(packs.packages[ilevel].GetScalarField(ivar)));
		}
	}
	//Revise Begin First Add avg value
	for (int ivar = 0; ivar < 5; ivar++)
	{
		this->hupdater + (&this->packs.packages[0].GetScalarField(ivar));
	}
	//Second Add avg slope 
	for (int ivar = 0; ivar < 15; ivar++)
	{
		this->hupdater + (&this->packs.avg_slope.GetScalarField(ivar));
	}
	//Third Add H distribution
	for (int ivar = 0; ivar < 5; ivar++)
	{
		this->hupdater + (&this->packs.pack[ivar]);
	}
	//Revise End
}

void GMGParallelUpdater::Update(int level,int order)
{
	//Revise Begin
	// if order==3 and at level 0 use hupdater 
	if ((order == 3) && (level == 0))
	{
		this->hupdater.UpdateFields();
	}
	else
	{
		this->updaters[level].Update();
	}
	//Revise End
}

void GMGParallelUpdater::UpdateHighOrderDistribution()
{
	this->hupdater.UpdateElementFields();
}
