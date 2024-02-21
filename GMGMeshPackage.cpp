#include "GMGMeshPackage.h"

GMGMeshPackage::GMGMeshPackage(std::vector<SubFvZone*>& mesh)
	:
	grids(mesh)
{
	this->nMGLevel = mesh.size();
}
