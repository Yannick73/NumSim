#include "pressure_solver/cg.h"

CG::CG(std::shared_ptr<PartitionShell> partition, double epsilon, int maximumNumberOfIterations) :
       PressureSolver(partition, epsilon, maximumNumberOfIterations),
       ma_(discretization_->nCells(), "M*A")
{
    discretization_->makeCGFields();
}

void CG::init()
{
    #pragma omp collapse(2)
    for(int j = 0; j < discretization_->pjN(); j++)
    {
        for(int i = 0; i < discretization_->piN(); i++)
        {
            // Wiki notation: Ax = b, our notation D2pDx2_ij + D2pDy2_ij = RHS_ij
            const double D2pDx2 = discretization_->computeD2pDx2(i,j);
            const double D2pDy2 = discretization_->computeD2pDy2(i,j);
            const double rhs    = discretization_->rhs(i,j);
            // which makes the r_ij = RHS_ik - D2pDx2_ij - D2pDy2_ij
            discretization_->r(i,j) = rhs - D2pDx2 - D2pDy2;
            discretization_->a(i,j) = discretization_->r(i,j);

        }
    }
}

void CG::step()
{
    if(iteration_ == 0)
        init();

    partition_->exchangeRA();
    
    // calculate alpha
    // basically every entry is squared and the result summed
    double rTr_sum_local  = 0;
    double aTMa_sum_local = 0;
    #pragma omp collapse(2)
    for(int j = 0; j < discretization_->pjN(); j++)
    {
        for(int i = 0; i < discretization_->piN(); i++)
        {
            const double r = discretization_->r(i,j);;
            rTr_sum_local += r*r;
            // see comment in init
            const double D2aDx2 = discretization_->computeD2aDx2(i,j);
            const double D2aDy2 = discretization_->computeD2aDy2(i,j);
            ma_(i,j) = (D2aDx2 + D2aDy2) * discretization_->a(i,j);
            aTMa_sum_local += discretization_->a(i,j) * ma_(i,j);
        }
    }
    double rTr, aTMa;
    MPI_Allreduce(&rTr_sum_local,  &rTr,  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&aTMa_sum_local, &aTMa, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    const double alpha = rTr / aTMa;

    // iterate x (p in our case) and r, and rTr anew
    double rTr_sum_local_new = 0;
    #pragma omp collapse(2)
    for(int j = 0; j < discretization_->pjN(); j++)
    {
        for(int i = 0; i < discretization_->piN(); i++)
        {
            discretization_->p(i,j) += alpha * discretization_->a(i,j);
            discretization_->r(i,j) -= alpha * ma_(i,j);
            rTr_sum_local_new += discretization_->r(i,j)*discretization_->r(i,j);
        }
    }
    double rTr_new;
    MPI_Allreduce(&rTr_sum_local_new,  &rTr_new,  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    double beta = rTr_new / rTr;
    // update a (p in Wiki)
    #pragma omp collapse(2)
    for(int j = 0; j < discretization_->pjN(); j++)
    {
        for(int i = 0; i < discretization_->piN(); i++)
        {
            discretization_->a(i,j) = discretization_->r(i,j) + beta * discretization_->a(i,j);
        }
    }
}