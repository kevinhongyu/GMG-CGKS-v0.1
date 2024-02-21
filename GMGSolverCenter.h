#ifndef GMGSOLVERCENTER_H
#define GMGSOLVERCENTER_H
#pragma once
#include "GMGMeshPackage.h"
#include "GMGFieldPackage.h"
#include "GMGParallelUpdater.h"
#include "GMGBoundaryPackage.h"
#include "GMGOperator.h"
#include "GMGHighOrderDeltaT.h"
#include "GMGFluxSolverPackage.h"
#include "OutputResidual.h"
#include "OutputResidual.h"
#include <sstream>
#include <fstream>
#include <iomanip>

class GMGSolverCenter
{
public:
	GMGMeshPackage meshpacks;
	GMGFieldPackage fieldpacks;
	GMGParallelUpdater updaterpacks;
	GMGBoundaryPackage<FirstOrderKFVSBoundaryPackage> boundarypacks;
	GMGOperator mgOperator;
	GMGHighOrderDeltaT dtpacks;
	GMGFluxSolverPackage<FirstOrderUpdate, 1> solverpacks;
	int order{ 1 };
public:
	GMGSolverCenter
	(
		std::vector<SubFvZone*>& grids, 
		double R_,
		double r_,
		double CFL,
		double Mu_,
		int K_,
		int Stage,
		double C2,
		int orderG=1
	);
public:

	void InitFlowField(double* convar);

	void ZeroResidual(int level);

	//In this Function we use stored Res and W To update Res and W
	void Relaxation(int level,int step=0);

	void SaveOld(int level);

	//Exactly We load force term stored by rhs, use-rhs to res
	void LoadResiduals(int level);

	//Store Rhs by + Res
	void StoreRhsByResidual(int level);

	//exactly we init Res by -rhs
	void InitResiduals(int level);

	//We update residual but not update W. res -signal!!!
	void UpdateResiduals(int level);

	//Restrict W 
	void RestrictAllQ(int level);

	//restrict res to coarse mesh
	void RestrictResidual(int level);

	//use oldw=oldw-w store delta w
	void PutCorrection(int level);

	//correct fine grid zero order
	void CorrectFineGrid(int level);

	//use qold = qold+ qnow exactly its = qold, its useful?
	void PutCorrectionBack(int level);

	//use rhs recover res
	void RecoverResidual(int level);

	//MultiGrid Loop
	void MultiGrid(int level = 0,int step=0);

	//If is coarset
	bool ISCoarsestGrid(int level);

	void SolveOneStep(int step=0);

	void InterpolatFineGrid(int level);

	void MultiGridInitFlow(int level);

	void InitFlow();

	void WritePolyMesh(boost::mpi::communicator& word, Field<double>& f)
	{
		//this function is used to output the unstructured volume mesh


		std::stringstream name;
		name << "result" << ".plt" << std::endl;
		std::string s;
		name >> s;
		std::ofstream ofs(s);
		if (!ofs.is_open())
		{
			std::cout << "cannot write the file" << std::endl;
			exit(0);
		}
		int var = 0;
		ofs << "Variables=\"x\",\"y\",\"z\"," << std::endl;// 1-2
		ofs << "\"var\"" << std::endl;//3
		var++;
		FvZone* fvzone = (f.p_blockzone);
		int nbnodes = fvzone->Nodes.size();
		int nbelements = fvzone->interiorv_cellID.size() /*+ fvzone->nBFaces + fvzone->nIFaces*/;
		int nbfaces = fvzone->Faces.size();
		int facenodes = 0;
		for (int iface = 0; iface < fvzone->Faces.size(); iface++)
		{
			FvFace& face = fvzone->Faces[iface];
			facenodes += face.NodeIDs.GetLengh();
		}

		ofs << "ZONE T=\"Finite Element Data\"" << std::endl;
		ofs << "Nodes=" << nbnodes << " ,Faces=" << nbfaces << " ,Elements=" << nbelements;
		ofs << "ZONETYPE=FEPolyhedron DATAPACKING=BLOCK" << std::endl;
		ofs << "TotalNumFaceNodes=" << facenodes << ", NumConnectedBoundaryFaces=0, TotalNumBoundaryConnections=0" << std::endl;
		ofs << "VarLocation=(NODAL, NODAL, NODAL, ";
		ofs << "CellCentered)" << std::endl;
		//ofs << std::setiosflags(std::ios::scientific) << std::setiosflags(std::ios::right) << std::setprecision(15);

		ofs << "#node_location_x,y,z" << std::endl;


		for (std::size_t i = 0; i < fvzone->Nodes.size(); i++)
		{
			ofs << fvzone->Nodes[i].X << " ";
			if (i % 8 == 0)
			{
				ofs << std::endl;
			}
		}
		ofs << std::endl;

		for (std::size_t i = 0; i < fvzone->Nodes.size(); i++)
		{
			ofs << fvzone->Nodes[i].Y << " ";
			if (i % 8 == 0)
			{
				ofs << std::endl;
			}
		}
		ofs << std::endl;

		for (std::size_t i = 0; i < fvzone->Nodes.size(); i++)
		{
			ofs << fvzone->Nodes[i].Z << " ";
			if (i % 8 == 0)
			{
				ofs << std::endl;
			}
		}
		ofs << std::endl;
		ofs << std::endl;

		for (std::size_t i = 0; i < fvzone->interiorv_cellID.size(); i++)
		{
			ofs << f.elementField.v_value[i] << " ";
			if (i % 8 == 0)
			{
				ofs << std::endl;
			}
		}
		ofs << std::endl;
		ofs << std::endl;

		for (int iface = 0; iface < fvzone->Faces.size(); iface++)
		{
			FvFace& face = fvzone->Faces[iface];
			ofs << face.NodeIDs.GetLengh() << std::endl;
		}
		ofs << std::endl;
		ofs << std::endl;

		for (int iface = 0; iface < fvzone->Faces.size(); iface++)
		{
			FvFace& face = fvzone->Faces[iface];
			for (int j = 0; j < face.NodeIDs.GetLengh(); j++)
			{
				ofs << face.NodeIDs[j] + 1 << " ";
			}
			ofs << std::endl;
		}
		ofs << std::endl;
		ofs << std::endl;

		std::vector<int> owner;
		std::vector<int> neighbor;
		int ghostcellpara = fvzone->interiorv_cellID.size() + fvzone->nBFaces;
		for (std::size_t i = 0; i < fvzone->interiorv_faceID.size(); i++)
		{
			int FaceID = fvzone->interiorv_faceID.at(i);
			FvFace& face = fvzone->Faces.at(FaceID);
			int own = face.OSideCell;
			int nei = face.NSideCell;
			owner.push_back(own + 1);
			neighbor.push_back(nei + 1);
		}
		/*for (std::size_t i = 0; i < fvzone->BoundThreads.size(); i++)
		{
			FvThread* p_thread = fvzone->BoundThreads.at(i).get();
			for (std::size_t j = 0; j < p_thread->v_faceID.size(); j++)
			{
				int FaceID = p_thread->v_faceID[j];
				FvFace& face = fvzone->Faces.at(FaceID);
				int own = face.OSideCell;
				int nei = face.NSideCell;
				owner.push_back(own+1);
				neighbor.push_back(0);
			}
		}*/
#if defined(_BaseParallelMPI_)
		SubFvZone* sub = dynamic_cast<SubFvZone*>(fvzone);
		std::vector<std::vector<PairedElem>>& vv_pair = sub->smep_GhostElemPair.vv_pairedElems;
		for (std::size_t i = 0; i < vv_pair.size(); i++)
		{
			for (std::size_t j = 0; j < vv_pair[i].size(); j++)
			{
				int FaceID = vv_pair[i][j].faceID;
				FvFace& face = sub->Faces.at(FaceID);
				int EleID = vv_pair[i][j].localID.first;
				/*bool isown = vv_pair[i][j].bOwnFlag;
				if (isown == true)
				{
					owner.push_back(EleID);
					neighbor.push_back(ghostcellpara);
				}
				else
				{
					owner.push_back(ghostcellpara);
					neighbor.push_back(EleID);
				}*/
				//owner.push_back(EleID+1);
				//neighbor.push_back(0);
				face.OSideCell = EleID;
				ghostcellpara++;
			}
		}
#endif;
		for (std::size_t i = sub->interiorv_faceID.size(); i < sub->nFaces - sub->nBFaces; i++)
		{
			FvFace& face = sub->Faces[i];
			owner.push_back(face.OSideCell + 1);
			neighbor.push_back(0);
		}

		for (std::size_t i = 0; i < fvzone->BoundThreads.size(); i++)
		{
			FvThread* p_thread = fvzone->BoundThreads.at(i).get();
			for (std::size_t j = 0; j < p_thread->v_faceID.size(); j++)
			{
				int FaceID = p_thread->v_faceID[j];
				FvFace& face = fvzone->Faces.at(FaceID);
				int own = face.OSideCell;
				int nei = face.NSideCell;
				owner.push_back(own + 1);
				neighbor.push_back(0);
			}
		}


		int owncount = 0;
		for (auto x : owner)
		{
			ofs << x << " ";
			if (owncount % 8 == 0)
			{
				ofs << std::endl;
			}
			owncount++;
		}
		ofs << std::endl;

		int neicount = 0;
		for (auto x : neighbor)
		{
			ofs << x << " ";
			if (neicount % 8 == 0)
			{
				ofs << std::endl;
			}
			neicount++;
		}
		ofs << std::endl;

		/*for (int iface = 0; iface < fvzone->Faces.size(); iface++)
		{
			FvFace& face = fvzone->Faces[iface];
			ofs << face.OSideCell << " ";
		}
		ofs << std::endl;

		for (int iface = 0; iface < fvzone->Faces.size(); iface++)
		{
			FvFace& face = fvzone->Faces[iface];
			ofs << face.NSideCell << " ";
		}*/
		ofs << std::endl;

		ofs.close();
	}

	void SetFinestOrder(int FOrder)
	{
		this->order = FOrder;
	}

};

#endif;