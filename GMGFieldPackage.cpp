#include "GMGFieldPackage.h"

GMGFieldPackage::GMGFieldPackage(GMGMeshPackage& mesh)
	:
	//Revise Begin,Add cf
	meshpack(mesh),cf(mesh.grids[0],1.0,"cf")
	//Revise End
{
	this->packages.resize(meshpack.grids.size());
	this->fields.resize(meshpack.grids.size());
	//First Get iLevel Mesh Then i pack Create Field
	for (std::size_t imesh = 0; imesh < meshpack.grids.size(); imesh++)
	{
		SubFvZone* sub = meshpack.grids[imesh];
		DataPackage& pack = this->packages[imesh];
		pack.p_blockzone = sub;
		this->fields[imesh].push_back(std::move(Field<double>(sub, 0.0, "Rho")));
		this->fields[imesh].push_back(std::move(Field<double>(sub, 0.0, "RhoU")));
		this->fields[imesh].push_back(std::move(Field<double>(sub, 0.0, "RhoV")));
		this->fields[imesh].push_back(std::move(Field<double>(sub, 0.0, "RhoW")));
		this->fields[imesh].push_back(std::move(Field<double>(sub, 0.0, "RhoE")));
		for (std::size_t ivar = 0; ivar < 5; ivar++)
		{
			pack.PackScalarField(fields[imesh][ivar], ivar);
		}
	}
	//Revise Begin Construct HighOrder WDis
	//First Put Them in vector
	for (int ivar = 0; ivar < 5; ivar++)
	{
		this->w_dis.push_back(ElementField<HighOrderDistribution>(meshpack.grids[0], HighOrderDistribution()));
	}
	//Second Put them in highorder gks data package
	for (int ivar = 0; ivar < 5; ivar++)
	{
		this->pack + w_dis[ivar];
	}
	//Revise End

	//Revise Begin Construct avg_slope
	//First put 15 blank field
	for (int ivar = 0; ivar < 15; ivar++)
	{
		this->slopes.push_back(Field<double>(this->meshpack.grids[0], 0.0, "Der"));
	}
	//Second Pack the fields
	for (int ivar = 0; ivar < 15; ivar++)
	{
		this->avg_slope.PackScalarField(this->slopes[ivar], ivar);
	}
	//Revise End

}
