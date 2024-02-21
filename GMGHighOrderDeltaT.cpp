#include "GMGHighOrderDeltaT.h"

GMGHighOrderDeltaT::GMGHighOrderDeltaT
(
	double R_,
	double r_,
	double CFL,
	double Mu_,
	int K_,
	GMGFieldPackage& fi
)
	:
	fields(fi)
{
	;
}
