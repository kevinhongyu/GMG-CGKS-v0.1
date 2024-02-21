#include "GMGOperator.h"

GMGOperator::GMGOperator(GMGFieldPackage& f)
    :
    fields(f)
{
    ;
}

void GMGOperator::RestrictAllQ(int level)
{
    SubFvZone* grid = fields.meshpack.grids[level];
    SubFvZone* coarseGrid = grid->cGrid;
    if (coarseGrid)
    {
        for (std::size_t ivar = 0; ivar < 5; ivar++)
        {
            RestrictQ(grid, fields.packages[level].GetScalarField(ivar),
                coarseGrid, fields.packages[level + 1].GetScalarField(ivar));
        }
    }
}

void GMGOperator::RestrictResidual(SubFvZone* fgrid_in, ElementField<std::vector<Flux3d>>& fres,
    SubFvZone* cgrid_in, ElementField<std::vector<Flux3d>>& cres)
{

    SubFvZone* fgrid = fgrid_in;
    SubFvZone* cgrid = cgrid_in;

    int  fnTotalCell = fgrid->interiorv_cellID.size();
    int  cnTotalCell = cgrid->interiorv_cellID.size();
    int* cell2coarsegridcell = fgrid->cell2coarsegridcell;


    for (int ivar = 0; ivar < 5; ivar++)
    {
        for (int iCell = 0; iCell < cnTotalCell; ++iCell)
        {
            cres.v_value[iCell][0].x[ivar] = 0.0;
        }
    }

    
    for (int ivar = 0; ivar < 5; ivar++)
    {
        for (int iCell = 0; iCell < fnTotalCell; ++iCell)
        {
            cres.v_value[cell2coarsegridcell[iCell]][0].x[ivar] -= fres.v_value[iCell][0].x[ivar];
            /*if (boost::mpi::communicator().rank() == 0)
            {
                if (ivar == 0)
                {
                    if (cell2coarsegridcell[iCell] == 1755)
                    {
                        std::cout << "iCell = " << iCell << std::endl;
                        std::cout << "Fine Res = " << fres.v_value[iCell][0].x[ivar] << std::endl;
                    }
                }
            }*/
        }
        
    }

    //Revise
    //for (int ivar = 0; ivar < 5; ivar++)
    //{
    //    for (int iCell = 0; iCell < cnTotalCell; ++iCell)
    //    {
    //        cres.v_value[iCell][0].x[ivar] = 0.0;
    //    }
    //}

    //for (int ivar = 0; ivar < 5; ivar++)
    //{
    //    for (int iCell = 0; iCell < fnTotalCell; ++iCell)
    //    {
    //        cres.v_value[cell2coarsegridcell[iCell]][0].x[ivar] -= fres.v_value[iCell][0].x[ivar] * fgrid->Cells[iCell].Volumn;
    //    }
    //}

    //for (int ivar = 0; ivar < 5; ivar++)
    //{
    //    for (int iCell = 0; iCell < cnTotalCell; ++iCell)
    //    {
    //        //cq.elementField.v_value[iCell] /= cvol[iCell];
    //        cres.v_value[cell2coarsegridcell[iCell]][0].x[ivar] /= cgrid->Cells[iCell].Volumn;
    //    }
    //}



    
    

}

void GMGOperator::PutCorrection(int level)
{
    SubFvZone* grid = fields.meshpack.grids[level];
  
    for (std::size_t ivar = 0; ivar < 5; ivar++)
    {
        PutCorrection(fields.packages[level].GetScalarField(ivar));
    }
    
}

void GMGOperator::CorrectFineGrid(int level)
{


    SubFvZone* grid = fields.meshpack.grids[level];
    SubFvZone* coarseGrid = grid->cGrid;
    if (coarseGrid)
    {
        for (std::size_t ivar = 0; ivar < 5; ivar++)
        {
            CorrectFineGrid(grid, fields.packages[level].GetScalarField(ivar),
                coarseGrid, fields.packages[level + 1].GetScalarField(ivar));
        }
    }



}

void GMGOperator::InterpolatFineGrid(int level)
{
    SubFvZone* grid = fields.meshpack.grids[level];
    SubFvZone* coarseGrid = grid->cGrid;

     if (coarseGrid)
     {
         for (std::size_t ivar = 0; ivar < 5; ivar++)
         {
             InterpolatQ(grid, fields.packages[level].GetScalarField(ivar),
                 coarseGrid, fields.packages[level + 1].GetScalarField(ivar));
         }
     }
    
}

void GMGOperator::CorrectFineGrid(SubFvZone* fgrid_in, Field<double>& fq, SubFvZone* cgrid_in, Field<double>& cq)
{
    SubFvZone* fgrid = fgrid_in;
    SubFvZone* cgrid = cgrid_in;

    int fnTCell = fgrid->interiorv_cellID.size();
    int* cell2coarsegridcell = fgrid->cell2coarsegridcell;

    double mgCorrectionLimit = 0.01;

    for (int iCell = 0; iCell < fnTCell; ++iCell)
    {
        int coarseCell = cell2coarsegridcell[iCell];

        double fgridValue = fq.elementField.v_value[iCell];
        double correction = cq.elementField0.v_value[coarseCell];

        double DQ = MIN(fabs(correction), fabs(fgridValue) * mgCorrectionLimit);
       fq.elementField.v_value[iCell] -= SIGN(DQ, correction);
       //fq.elementField.v_value[iCell] -= correction;
       /* if (boost::mpi::communicator().rank() == 0)
        {
            std::cout << "Correction = " << correction << std::endl;
        }*/

    }
}
