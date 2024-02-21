#ifndef GMGOPERATOR_H
#define GMGOPERATOR_H
#pragma once
#include "GMGFieldPackage.h"
#include "Flux3d.h"
#include <boost/mpi.hpp>
#include "AppRes.h"

class GMGOperator
{
public:
    GMGFieldPackage& fields;
public:
    GMGOperator(GMGFieldPackage& f);
public:

    void RestrictAllQ(int level);

    void RestrictResidual(SubFvZone* fgrid_in, ElementField<std::vector<Flux3d>>& fres,
        SubFvZone* cgrid_in, ElementField<std::vector<Flux3d>>& cres);

    void PutCorrection(int level);

    void CorrectFineGrid(int level);

    void InterpolatFineGrid(int level);

    //this is called by fine grid to use coarse grid's delta w,but we need store old first!
    //When to saveold need to be consider!
    template<typename T>
    void PutCorrection(Field<T>& cq)
    {
        cq.elementField0 = cq.elementField0 - cq.elementField;
       /* if (boost::mpi::communicator().rank() == 0)
        {
            std::cout << "correction ele0 = " << cq.elementField0.v_value[0] << std::endl;
        }*/
    }


    template <typename T>
    void RestrictQ(SubFvZone* fgrid_in, Field<T>& fq, SubFvZone* cgrid_in, Field<T>& cq)
    {
        SubFvZone* fgrid = fgrid_in;
        SubFvZone* cgrid = cgrid_in;

        int      fnTotalCell = fgrid->interiorv_cellID.size();
        int      cnTotalCell = cgrid->interiorv_cellID.size();
        int* cell2coarsegridcell = fgrid->cell2coarsegridcell;

        for (int iCell = 0; iCell < cnTotalCell; ++iCell)
        {
            cq.elementField.v_value[iCell] = 0.0;
        }

        for (int iCell = 0; iCell < fnTotalCell; ++iCell)
        {
            //cq.elementField.v_value[cell2coarsegridcell[iCell]] += fvol[iCell] * fq[iCell];
            cq.elementField.v_value[cell2coarsegridcell[iCell]] += fq.elementField.v_value[iCell]
                * fgrid->Cells[iCell].Volumn;
           /* if (boost::mpi::communicator().rank() == 0)
            {
                if (cell2coarsegridcell[iCell] == 1755)
                {
                    std::cout << "Fine CellID = "<<iCell << std::endl;
                    std::cout << "Fine Cell Value = " << fq.elementField.v_value[iCell] << std::endl;
                    std::cout << "Fine Cell Volum = " << fgrid->Cells[iCell].Volumn << std::endl;
                    std::cout << "Coarse Cell Value = " << cq.elementField.v_value[cell2coarsegridcell[iCell]] << std::endl;
                }
            }*/
           
        }

        for (int iCell = 0; iCell < cnTotalCell; ++iCell)
        {
            //cq.elementField.v_value[iCell] /= cvol[iCell];
            cq.elementField.v_value[iCell] /= cgrid->Cells[iCell].Volumn;
           /* if (boost::mpi::communicator().rank() == 0)
            {
                if (iCell == 1755)
                {
                    std::cout << "Coarse Cell Volum = " << cgrid->Cells[iCell].Volumn << std::endl;
                    std::cout << "Coarse Cell Value = " << cq.elementField.v_value[iCell] << std::endl;
                }
            }*/
        }
    }

    void CorrectFineGrid(SubFvZone* fgrid_in, Field<double>& fq, SubFvZone* cgrid_in, Field<double>& cq);


    template < typename T >
    void InterpolatQ(SubFvZone* fgrid_in, Field<T>& fq, SubFvZone* cgrid_in, Field<T>& cq)
    {
        SubFvZone* fgrid = fgrid_in;

        int fnTotalCell = fgrid->interiorv_cellID.size();

        int* cell2coarsegridcell = fgrid->cell2coarsegridcell;
        int cc;
        for (int iCell = 0; iCell < fnTotalCell; ++iCell)
        {
            cc = cell2coarsegridcell[iCell];

            fq.elementField.v_value[iCell] = cq.elementField.v_value[cc];
        }
    }

};

template < typename T >
inline T MIN(const T& a, const T& b)
{
    return (a < b) ? a : b;
}

template < typename T >
inline T ABS(const T& a)
{
    return (a < 0) ? -a : a;
}

template < typename T >
inline T SIGN(const T& a, const T& b)
{
    int s = 1;
    if (b < 0) s = -1;
    return s * ABS(a);
}



#endif;