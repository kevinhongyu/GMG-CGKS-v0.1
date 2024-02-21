#include "GMGSolverCenter.h"

GMGSolverCenter::GMGSolverCenter
(
	std::vector<SubFvZone*>& grids, 
	double R_,
	double r_,
	double CFL,
	double Mu_,
	int K_,
	int Stage,
	double C2,
	int orderG
)
	:
	meshpacks(grids),fieldpacks(meshpacks),updaterpacks(fieldpacks),
	boundarypacks(fieldpacks),mgOperator(fieldpacks),dtpacks(R_,r_,CFL,Mu_,K_,
		fieldpacks),solverpacks(fieldpacks,boundarypacks),order(orderG)
{
	//this is original
	/*for (int ilevel = 0; ilevel < fieldpacks.meshpack.nMGLevel; ilevel++)
	{
		dtpacks.DTs.push_back(HighOrderTransient(R_, r_, CFL, Mu_, K_, solverpacks.fluxSolvers[ilevel].block,
			fieldpacks.meshpack.grids[ilevel], fieldpacks.packages[ilevel]));
	}*/
	//original end
	//Revise Begin,for level 0 third order
	if (order == 3)
	{
		dtpacks.DTs.push_back(HighOrderTransient(R_, r_, CFL, Mu_, K_, solverpacks.thirdorderfluxSolver.block,
			fieldpacks.meshpack.grids[0], fieldpacks.packages[0]));

		for (int ilevel = 1; ilevel < fieldpacks.meshpack.nMGLevel; ilevel++)
		{
			dtpacks.DTs.push_back(HighOrderTransient(R_, r_, CFL, Mu_, K_, solverpacks.fluxSolvers[ilevel].block,
				fieldpacks.meshpack.grids[ilevel], fieldpacks.packages[ilevel]));
		}
	}
	else
	{
		for (int ilevel = 0; ilevel < fieldpacks.meshpack.nMGLevel; ilevel++)
		{
			dtpacks.DTs.push_back(HighOrderTransient(R_, r_, CFL, Mu_, K_, solverpacks.fluxSolvers[ilevel].block,
				fieldpacks.meshpack.grids[ilevel], fieldpacks.packages[ilevel]));
		}
	}
	//revise end


	for (int ilevel = 0; ilevel < fieldpacks.meshpack.nMGLevel; ilevel++)
	{
		solverpacks.SetSolverParameter(Mu_, r_, R_, C2, K_, Stage);
	}
}

void GMGSolverCenter::InitFlowField(double* convar)
{
	for (int ivar = 0; ivar < 5; ivar++)
	{
		for (int ilevel = 0; ilevel < fieldpacks.meshpack.nMGLevel; ilevel++)
		{
			this->fieldpacks.packages[ilevel].
				GetScalarField(ivar).elementField.Initialize(convar[ivar]);
		}
	}
	//Revise Begin,if third order set avg_slope=0.0,cf=1.0
	if (this->order == 3)
	{
		for (int ivar = 0; ivar < 15; ivar++)
		{
			this->fieldpacks.avg_slope.GetScalarField(ivar).elementField.Initialize(0.0);
		}
		this->fieldpacks.cf.elementField.Initialize(1.0);
	}
	//Revise End
	
}

void GMGSolverCenter::ZeroResidual(int level)
{
	//this is original
	//this->solverpacks.fluxSolvers[level].ZeroResidual();
	//original end
	// 
	//Revise Begin,use highordersolver to zero residual,remember only level=0 will use. after this all level related function will note
	if ((this->order == 3)&&(level==0))
	{
		this->solverpacks.thirdorderfluxSolver.ZeroResidual();
	}
	else
	{
		this->solverpacks.fluxSolvers[level].ZeroResidual();
	}
	//Revise End
}

void GMGSolverCenter::Relaxation(int level,int step)
{
	////this is original version
	////1 Cal deltaT
	//this->dtpacks.DTs[level].GetDeltaTForGKSInvicid();
	////2 Save this level's element0
	//SaveOld(level);
	////3 Update Parallel InterFace
	//this->updaterpacks.updaters[level].Update();
	////4 SetBounddary Value
	//this->boundarypacks.boundaryPacks[level].SetOside();
	//this->boundarypacks.boundaryPacks[level].SetNside(this->dtpacks.DTs[level].K);
	////5 UpdateResidual  previous residual is not zero,but Here We should first load Res
	//int stage = 0;
	//LoadResiduals(level);
	////this->solverpacks.fluxSolvers[level].ZeroResidual();
	//this->solverpacks.fluxSolvers[level].UpdateResiduals(this->dtpacks.DTs[level].K, stage,0);
	//this->solverpacks.fluxSolvers[level].UpdateFlowField(stage);
	////original end

	//Revise if third order,use third order update
	if ((order == 3) && (level == 0))
	{
		//1 Cal delta T
		this->dtpacks.DTs[level].GetDeltaTForGKSInvicid();
		
		//std::cout << "dt = " << this->dtpacks.DTs[level].GetDeltaTForGKSInvicid() << std::endl;
		//std::cout << "dt = " << this->solverpacks.thirdorderfluxSolver.block.dt << std::endl;
		//2 Save Old convar field
		SaveOld(level);
		//3 set boundary value for ghost cell and ghost gradient
		this->boundarypacks.thirdorderPacks.SetGhostConvar();
		this->boundarypacks.thirdorderPacks.SetGhostGradient
		(
			this->fieldpacks.avg_slope.GetScalarField(0),
			this->fieldpacks.avg_slope.GetScalarField(5),
			this->fieldpacks.avg_slope.GetScalarField(10),
			this->fieldpacks.avg_slope.GetScalarField(1),
			this->fieldpacks.avg_slope.GetScalarField(6),
			this->fieldpacks.avg_slope.GetScalarField(11),
			this->fieldpacks.avg_slope.GetScalarField(2),
			this->fieldpacks.avg_slope.GetScalarField(7),
			this->fieldpacks.avg_slope.GetScalarField(12),
			this->fieldpacks.avg_slope.GetScalarField(3),
			this->fieldpacks.avg_slope.GetScalarField(8),
			this->fieldpacks.avg_slope.GetScalarField(13),
			this->fieldpacks.avg_slope.GetScalarField(4),
			this->fieldpacks.avg_slope.GetScalarField(9),
			this->fieldpacks.avg_slope.GetScalarField(14)
		);
		//4 update fields
		this->updaterpacks.Update(level, order);
		//5 reconstruct
		for (int ivar = 0; ivar < 5; ivar++)
		{
			hweno3_reconstruction_hybrid_with_limiter
			(
				2,
				this->fieldpacks.cf,
				this->fieldpacks.packages[level].GetScalarField(ivar),
				this->fieldpacks.pack[ivar],
				this->fieldpacks.avg_slope.GetScalarField(ivar),
				this->fieldpacks.avg_slope.GetScalarField(ivar + 5),
				this->fieldpacks.avg_slope.GetScalarField(ivar + 10)
			);
		}
		//6 update highorder pack
		this->updaterpacks.UpdateHighOrderDistribution();
		//7 Set Boundary ON Value
		this->boundarypacks.thirdorderPacks.SetOside();
		this->boundarypacks.thirdorderPacks.SetNside(this->dtpacks.DTs[level].K);
		//8 load residuals
		int stage = 0;
		LoadResiduals(level);
		//9 update residuals
		this->solverpacks.thirdorderfluxSolver.UpdateResidual(this->dtpacks.DTs[level].K, stage);
		//this->solverpacks.thirdorderfluxSolver.ResidualSmooth(0.2,10);
		//10 update flowfield
		this->solverpacks.thirdorderfluxSolver.UpdateFlowField(this->dtpacks.DTs[level].K, stage,step);

		//this->solverpacks.thirdorderfluxSolver.UpdateFlowFieldLocal();

		//this->solverpacks.thirdorderfluxSolver.UpdateFluxesAndCF(this->dtpacks.DTs[level].K, stage);

	}
	else
	{
		//1 Cal deltaT
		this->dtpacks.DTs[level].GetDeltaTForGKSInvicid();
		//this->dtpacks.DTs[level].block.dt = this->dtpacks.DTs[0].block.dt;
		//2 Save this level's element0
		//SaveOld(level);
		//3 Update Parallel InterFace
		this->updaterpacks.updaters[level].Update();
		//4 SetBounddary Value
		this->boundarypacks.boundaryPacks[level].SetOside();
		this->boundarypacks.boundaryPacks[level].SetNside(this->dtpacks.DTs[level].K);
		//5 UpdateResidual  previous residual is not zero,but Here We should first load Res
		int stage = 0;
		LoadResiduals(level);
		//this->solverpacks.fluxSolvers[level].ZeroResidual();
		this->solverpacks.fluxSolvers[level].UpdateResiduals(this->dtpacks.DTs[level].K, stage, 0);
		this->solverpacks.fluxSolvers[level].UpdateFlowField(stage);
		//rEVISE tRY lOCAL step
		//this->solverpacks.fluxSolvers[level].UpdateFlowFieldLocal();
	}
}

void GMGSolverCenter::SaveOld(int level)
{
	for (int ivar = 0; ivar < 5; ivar++)
	{
		this->fieldpacks.packages[level].GetScalarField(ivar).SaveOld();
	}
}

void GMGSolverCenter::LoadResiduals(int level)
{
	////this is original
	//for (int ivar = 0; ivar < 5; ivar++)
	//{
	//	for (std::size_t iCell = 0; iCell < this->meshpacks.grids[level]->interiorv_cellID.size()
	//		; iCell++)
	//	{
	//		this->solverpacks.fluxSolvers[level].ResFlux.v_value[iCell][0].x[ivar]
	//			= -this->solverpacks.fluxSolvers[level].RhsFlux.v_value[iCell][0].x[ivar];
	//	}
	//}
	////original end

	//Revise Begin, if third&&level==0, d will use "-" signal
	if ((order == 3) && (level == 0))
	{
		for (int ivar = 0; ivar < 5; ivar++)
		{
			for (std::size_t iCell = 0; iCell < this->meshpacks.grids[level]->interiorv_cellID.size()
				; iCell++)
			{
				this->solverpacks.thirdorderfluxSolver.ResFlux.v_value[iCell][0].x[ivar]
					= -this->solverpacks.thirdorderfluxSolver.RhsFlux.v_value[iCell][0].x[ivar];

				this->solverpacks.thirdorderfluxSolver.ResFlux.v_value[iCell][0].d[ivar]
					= -this->solverpacks.thirdorderfluxSolver.RhsFlux.v_value[iCell][0].d[ivar];

				this->solverpacks.thirdorderfluxSolver.ResFlux.v_value[iCell][0].d[ivar+5]
					= -this->solverpacks.thirdorderfluxSolver.RhsFlux.v_value[iCell][0].d[ivar+5];

				this->solverpacks.thirdorderfluxSolver.ResFlux.v_value[iCell][0].d[ivar + 10]
					= -this->solverpacks.thirdorderfluxSolver.RhsFlux.v_value[iCell][0].d[ivar + 10];
			}
		}
	}
	else
	{
		for (int ivar = 0; ivar < 5; ivar++)
		{
			for (std::size_t iCell = 0; iCell < this->meshpacks.grids[level]->interiorv_cellID.size()
				; iCell++)
			{
				this->solverpacks.fluxSolvers[level].ResFlux.v_value[iCell][0].x[ivar]
					= -this->solverpacks.fluxSolvers[level].RhsFlux.v_value[iCell][0].x[ivar];
			}
		}
	}
	//Revise end
}

void GMGSolverCenter::StoreRhsByResidual(int level)
{
	//this is first order original
	/*for (int ivar = 0; ivar < 5; ivar++)
	{
		for (std::size_t iCell = 0; iCell < this->meshpacks.grids[level]->interiorv_cellID.size()
			; iCell++)
		{
			this->solverpacks.fluxSolvers[level].RhsFlux.v_value[iCell][0].x[ivar]
				= this->solverpacks.fluxSolvers[level].ResFlux.v_value[iCell][0].x[ivar];
		}
	}*/
	//original end
	//Revisebegin
	if ((order == 3) && (level == 0))
	{
		for (int ivar = 0; ivar < 5; ivar++)
		{
			for (std::size_t iCell = 0; iCell < this->meshpacks.grids[level]->interiorv_cellID.size()
				; iCell++)
			{
				this->solverpacks.thirdorderfluxSolver.RhsFlux.v_value[iCell][0].x[ivar]
					= this->solverpacks.thirdorderfluxSolver.ResFlux.v_value[iCell][0].x[ivar];

				this->solverpacks.thirdorderfluxSolver.RhsFlux.v_value[iCell][0].d[ivar]
					= this->solverpacks.thirdorderfluxSolver.ResFlux.v_value[iCell][0].d[ivar];

				this->solverpacks.thirdorderfluxSolver.RhsFlux.v_value[iCell][0].d[ivar + 5]
					= this->solverpacks.thirdorderfluxSolver.ResFlux.v_value[iCell][0].d[ivar + 5];

				this->solverpacks.thirdorderfluxSolver.RhsFlux.v_value[iCell][0].d[ivar + 10]
					= this->solverpacks.thirdorderfluxSolver.ResFlux.v_value[iCell][0].d[ivar + 10];
			}
		}
	}
	else
	{
		for (int ivar = 0; ivar < 5; ivar++)
		{
			for (std::size_t iCell = 0; iCell < this->meshpacks.grids[level]->interiorv_cellID.size()
				; iCell++)
			{
				this->solverpacks.fluxSolvers[level].RhsFlux.v_value[iCell][0].x[ivar]
					= this->solverpacks.fluxSolvers[level].ResFlux.v_value[iCell][0].x[ivar];
			}
		}
	}
	//Revise end
}

void GMGSolverCenter::InitResiduals(int level)
{
	//this is original first order version
	/*for (int ivar = 0; ivar < 5; ivar++)
	{
		for (std::size_t iCell = 0; iCell < this->meshpacks.grids[level]->interiorv_cellID.size()
			; iCell++)
		{
			this->solverpacks.fluxSolvers[level].ResFlux.v_value[iCell][0].x[ivar]
				= -this->solverpacks.fluxSolvers[level].RhsFlux.v_value[iCell][0].x[ivar];
		}
	}*/
	//original end
	//Revise Begin
	if ((order == 3) && (level == 0))
	{
		for (int ivar = 0; ivar < 5; ivar++)
		{
			for (std::size_t iCell = 0; iCell < this->meshpacks.grids[level]->interiorv_cellID.size()
				; iCell++)
			{
				this->solverpacks.thirdorderfluxSolver.ResFlux.v_value[iCell][0].x[ivar]
					= -this->solverpacks.thirdorderfluxSolver.RhsFlux.v_value[iCell][0].x[ivar];

				this->solverpacks.thirdorderfluxSolver.ResFlux.v_value[iCell][0].d[ivar]
					= -this->solverpacks.thirdorderfluxSolver.RhsFlux.v_value[iCell][0].d[ivar];

				this->solverpacks.thirdorderfluxSolver.ResFlux.v_value[iCell][0].d[ivar+5]
					= -this->solverpacks.thirdorderfluxSolver.RhsFlux.v_value[iCell][0].d[ivar+5];

				this->solverpacks.thirdorderfluxSolver.ResFlux.v_value[iCell][0].d[ivar+10]
					= -this->solverpacks.thirdorderfluxSolver.RhsFlux.v_value[iCell][0].d[ivar+10];
			}
		}
	}
	else
	{
		for (int ivar = 0; ivar < 5; ivar++)
		{
			for (std::size_t iCell = 0; iCell < this->meshpacks.grids[level]->interiorv_cellID.size()
				; iCell++)
			{
				this->solverpacks.fluxSolvers[level].ResFlux.v_value[iCell][0].x[ivar]
					= -this->solverpacks.fluxSolvers[level].RhsFlux.v_value[iCell][0].x[ivar];
			}
		}
	}

}

void GMGSolverCenter::UpdateResiduals(int level)
{
	//this is original
	////1 Cal deltaT
	//this->dtpacks.DTs[level].GetDeltaTForGKSInvicid();
	////2 Update Parallel InterFace
	//this->updaterpacks.updaters[level].Update();
	////3 SetBounddary Value
	//this->boundarypacks.boundaryPacks[level].SetOside();
	//this->boundarypacks.boundaryPacks[level].SetNside(this->dtpacks.DTs[level].K);
	////4 UpdateResidual  previous residual is not zero,but Here We already Init res by -rhs
	//int stage = 0;
	//this->solverpacks.fluxSolvers[level].UpdateResiduals(this->dtpacks.DTs[level].K, stage,0);
	//original end
	//Revise begin
	if ((order == 3) && (level == 0))
	{
		//1 Cal delta T
		this->dtpacks.DTs[level].GetDeltaTForGKSInvicid();
		//2 Save Old convar field
		//SaveOld(level);
		//3 set boundary value for ghost cell and ghost gradient
		this->boundarypacks.thirdorderPacks.SetGhostConvar();
		this->boundarypacks.thirdorderPacks.SetGhostGradient
		(
			this->fieldpacks.avg_slope.GetScalarField(0),
			this->fieldpacks.avg_slope.GetScalarField(5),
			this->fieldpacks.avg_slope.GetScalarField(10),
			this->fieldpacks.avg_slope.GetScalarField(1),
			this->fieldpacks.avg_slope.GetScalarField(6),
			this->fieldpacks.avg_slope.GetScalarField(11),
			this->fieldpacks.avg_slope.GetScalarField(2),
			this->fieldpacks.avg_slope.GetScalarField(7),
			this->fieldpacks.avg_slope.GetScalarField(12),
			this->fieldpacks.avg_slope.GetScalarField(3),
			this->fieldpacks.avg_slope.GetScalarField(8),
			this->fieldpacks.avg_slope.GetScalarField(13),
			this->fieldpacks.avg_slope.GetScalarField(4),
			this->fieldpacks.avg_slope.GetScalarField(9),
			this->fieldpacks.avg_slope.GetScalarField(14)
		);
		//4 update fields
		this->updaterpacks.Update(level, order);
		//5 reconstruct
		for (int ivar = 0; ivar < 5; ivar++)
		{
			hweno3_reconstruction_hybrid_with_limiter
			(
				2,
				this->fieldpacks.cf,
				this->fieldpacks.packages[level].GetScalarField(ivar),
				this->fieldpacks.pack[ivar],
				this->fieldpacks.avg_slope.GetScalarField(ivar),
				this->fieldpacks.avg_slope.GetScalarField(ivar + 5),
				this->fieldpacks.avg_slope.GetScalarField(ivar + 10)
			);
		}
		//6 update highorder pack
		this->updaterpacks.UpdateHighOrderDistribution();
		//7 Set Boundary ON Value
		this->boundarypacks.thirdorderPacks.SetOside();
		this->boundarypacks.thirdorderPacks.SetNside(this->dtpacks.DTs[level].K);
		//8 load residuals
		int stage = 0;
		LoadResiduals(level);
		//9 update residuals
		this->solverpacks.thirdorderfluxSolver.UpdateResidual(this->dtpacks.DTs[level].K, stage);
	}
	else
	{
		//1 Cal deltaT
		this->dtpacks.DTs[level].GetDeltaTForGKSInvicid();
		//this->dtpacks.DTs[level].block.dt = this->dtpacks.DTs[0].block.dt;
		//2 Update Parallel InterFace
		this->updaterpacks.updaters[level].Update();
		//3 SetBounddary Value
		this->boundarypacks.boundaryPacks[level].SetOside();
		this->boundarypacks.boundaryPacks[level].SetNside(this->dtpacks.DTs[level].K);
		//4 UpdateResidual  previous residual is not zero,but Here We already Init res by -rhs
		int stage = 0;
		this->solverpacks.fluxSolvers[level].UpdateResiduals(this->dtpacks.DTs[level].K, stage, 0);
	}
	//Revise End
}

void GMGSolverCenter::RestrictAllQ(int level)
{
	this->mgOperator.RestrictAllQ(level);
}

void GMGSolverCenter::RestrictResidual(int level)
{
	//original version
	/*SubFvZone* fGrid = this->meshpacks.grids[level];
	SubFvZone* cGrid = this->meshpacks.grids[level + 1];
	ElementField<std::vector<Flux3d>>& fres = this->solverpacks.fluxSolvers[level].ResFlux;
	ElementField<std::vector<Flux3d>>& cres = this->solverpacks.fluxSolvers[level + 1].ResFlux;
	this->mgOperator.RestrictResidual(fGrid, fres, cGrid, cres);*/
	//original end
	//Revise Begin if(order&&level) fres=third, cres=flux[]
	SubFvZone* fGrid = this->meshpacks.grids[level];
	SubFvZone* cGrid = this->meshpacks.grids[level + 1];
	if ((order == 3) && (level == 0))
	{
		ElementField<std::vector<Flux3d>>& fres = this->solverpacks.thirdorderfluxSolver.ResFlux;
		ElementField<std::vector<Flux3d>>& cres = this->solverpacks.fluxSolvers[level + 1].ResFlux;
		this->mgOperator.RestrictResidual(fGrid, fres, cGrid, cres);
	}
	else
	{
		ElementField<std::vector<Flux3d>>& fres = this->solverpacks.fluxSolvers[level].ResFlux;
		ElementField<std::vector<Flux3d>>& cres = this->solverpacks.fluxSolvers[level + 1].ResFlux;
		this->mgOperator.RestrictResidual(fGrid, fres, cGrid, cres);
	}
	//Revise End

}

void GMGSolverCenter::PutCorrection(int level)
{
	this->mgOperator.PutCorrection(level);
}

void GMGSolverCenter::CorrectFineGrid(int level)
{
	this->mgOperator.CorrectFineGrid(level);
}

void GMGSolverCenter::PutCorrectionBack(int level)
{
	for (int ivar = 0; ivar < 5; ivar++)
	{
		this->fieldpacks.packages[level].GetScalarField(ivar).elementField0
			= this->fieldpacks.packages[level].GetScalarField(ivar).elementField0
			+ this->fieldpacks.packages[level].GetScalarField(ivar).elementField;
	}
}

void GMGSolverCenter::RecoverResidual(int level)
{
	//this is original
	/*for (int ivar = 0; ivar < 5; ivar++)
	{
		for (int iCell = 0; iCell < this->meshpacks.grids[level]->interiorv_cellID.size();
			iCell++)
		{
			this->solverpacks.fluxSolvers[level].ResFlux.v_value[iCell][0].x[ivar]
				= this->solverpacks.fluxSolvers[level].RhsFlux.v_value[iCell][0].x[ivar];
		}
	}*/
	//original end
	//Revise Begin
	if ((order == 3) && (level == 0))
	{
		for (int ivar = 0; ivar < 5; ivar++)
		{
			for (int iCell = 0; iCell < this->meshpacks.grids[level]->interiorv_cellID.size();
				iCell++)
			{
				this->solverpacks.thirdorderfluxSolver.ResFlux.v_value[iCell][0].x[ivar]
					= this->solverpacks.thirdorderfluxSolver.RhsFlux.v_value[iCell][0].x[ivar];

				this->solverpacks.thirdorderfluxSolver.ResFlux.v_value[iCell][0].d[ivar]
					= this->solverpacks.thirdorderfluxSolver.RhsFlux.v_value[iCell][0].d[ivar];

				this->solverpacks.thirdorderfluxSolver.ResFlux.v_value[iCell][0].d[ivar+5]
					= this->solverpacks.thirdorderfluxSolver.RhsFlux.v_value[iCell][0].d[ivar+5];

				this->solverpacks.thirdorderfluxSolver.ResFlux.v_value[iCell][0].d[ivar+10]
					= this->solverpacks.thirdorderfluxSolver.RhsFlux.v_value[iCell][0].d[ivar+10];
			}
		}
	}
	else
	{
		for (int ivar = 0; ivar < 5; ivar++)
		{
			for (int iCell = 0; iCell < this->meshpacks.grids[level]->interiorv_cellID.size();
				iCell++)
			{
				this->solverpacks.fluxSolvers[level].ResFlux.v_value[iCell][0].x[ivar]
					= this->solverpacks.fluxSolvers[level].RhsFlux.v_value[iCell][0].x[ivar];
			}
		}
	}
	//Revise End
}

void GMGSolverCenter::MultiGrid(int level,int step)
{
	boost::mpi::communicator word;
	StoreRhsByResidual(level);

	//! If grid is in the coarsest level.
	if (ISCoarsestGrid(level))
	{

		//Relaxation(level);
		SaveOld(level);
		for (int i = 0; i < 4; i++)
		{
			Relaxation(level);
		}
	}
	else
	{
		int cLevel = level + 1;


		//! Smooth the solutions before corrections.
		if (level > 0)
		{
			SaveOld(level);
			for (int i = 0; i < 2; i++)
			{
				Relaxation(level);
			}
			UpdateResiduals(level);
		}
		else
		{
		//Relaxation(level);
			SaveOld(level);
			for (int i = 0; i < 1; i++)
			{
				Relaxation(level,step);
			}
			UpdateResiduals(level);
		}
		//Relaxation(level);

		//! Restrict from fine to coarse grid for all q.
		//! Corresponding to nsmb: wsav = restr(w).
		//! Call the following function to update the Q value(cq) on the coarse grid.
		RestrictAllQ(level);



		//! Defect.
		//! The followings corresponding to nsmb d = RL-1(wsav) - restr(RL(w)-f).
		//! in which (*res) = - (*rhs) correspond to - f
		//! Note that the value of res is on fine grid, rhs ought be f£¬and this should be checked.
		//! Here should be ni and nj,rather than cni and cnj!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		//InitResiduals(level);

		//! After UpdateResiduals,(*res) = RL(w) - f.
		//UpdateResiduals(level);
		if (word.rank() == 0)
		{
			//std::cout << "Here11" << std::endl;
		}

		//! After RestrictDefect,(*cres) = - restr(RL(w) - f).
		RestrictResidual(level);

		//! After UpdateResiduals,(*cres) = RL-1(wsav) - restr(RL(w) - f).
		UpdateResiduals(cLevel);

		//! Solve on coarse grid.
		int MGFasType = 1;

		//! After MultiGrid(cgrid), the value of q on coarse grid has been updated.
		//! At this time,the value of q on coarse grid equals to w0 of NSMB.
		//! MGFasType = 1 V cycle
		//! MGFasType = 2 W cycle
		for (int i = 0; i < MGFasType; ++i)
		{
			MultiGrid(cLevel);
		}


		//! Put correction back to the coarse grid.
		//! Actually, PutCorrection has done the operation of (w0 - wsav).
		PutCorrection(cLevel);


		//! Correct variables in fine grid.
		//! Actually CorrectFineGrid has done the operation of w = w + prol(w0 - wsav).
		CorrectFineGrid(level);

		////! Smooth the solutions after corrections.
		//int n_post = GlobalDataBase::GetIntParaFromDB("n_post");
		for (int iter = 0; iter < 1; ++iter)
		{
			Relaxation(level,step);
		}

		//! The following function actually recover the value of q on coarse grid.
		PutCorrectionBack(cLevel);

	}
	//! Put RHS back into its places.
	//! Here, we just recover the initial value of res.
	RecoverResidual(level);
}

bool GMGSolverCenter::ISCoarsestGrid(int level)
{
	SubFvZone* grid = this->meshpacks.grids[level];
	if (grid->cGrid)
	{
		return false;
	}
	else
	{
		return true;
	}
}

void GMGSolverCenter::SolveOneStep(int step)
{

	//! Zero the residuals of the finest grid.
	ZeroResidual(0);

	MultiGrid(0,step);

}

void GMGSolverCenter::InterpolatFineGrid(int level)
{
	this->mgOperator.InterpolatFineGrid(level);
	this->updaterpacks.updaters[level].Update();
}

void GMGSolverCenter::MultiGridInitFlow(int level)
{


	ZeroResidual(level);

	StoreRhsByResidual(level);

	//! Relaxation.
	Relaxation(level);
}

void GMGSolverCenter::InitFlow()
{
	OutputResidual outputresidualcoarse(fieldpacks.packages[1], 1, "D:/");
	OutputResidual outputresidualcoarse1(fieldpacks.packages[2], 1, "D:/");
	for (int iLevel = this->meshpacks.nMGLevel - 1; iLevel >= 0; --iLevel)
	{
		int levelStep = 0;
		int totalLevelStep = 400 * (iLevel+1);



		while (levelStep < totalLevelStep)
		{
			++levelStep;
			//outputresidualcoarse.Output(levelStep);
			if (boost::mpi::communicator().rank() == 0)
			{
				std::cout << "Coarset Level = " << std::endl;
			}
			//outputresidualcoarse1.Output(levelStep);
			

		   MultiGridInitFlow(iLevel);
		   if (boost::mpi::communicator().rank() == 0)
		   {
			   std::cout << "ilevel init = " << iLevel << "level step = " << levelStep << std::endl;
		   }


			
		}

		//! Interpolate from coarse grid to fine grid.
		
		if (iLevel >= 1)
		{
			InterpolatFineGrid(iLevel - 1);
		}


		
	}
}
