#include "pressure_solver.h"

PressureSolver::PressureSolver(const std::shared_ptr<Discretization> discretization, 
                               double epsilon, int maximumNumberOfIterations) : 
                               discretization_(discretization),
                               epsilon_(epsilon),
                               maximumNumberOfIterations_(maximumNumberOfIterations)
{
    assert(discretization != nullptr);
    assert(epsilon > 0);
    assert(maximumNumberOfIterations > 0);
}

inline void PressureSolver::setBoundaryValues()
{
    const int iEnd = discretization_->iEnd();
    const int jEnd = discretization_->jEnd();
    // left and right edges take priority
    for(int j = -1; j < jEnd+1; j++)
    {
        // set the left edge
        discretization_->p(-1,j) = discretization_->p(0,j);
        // set the right edge
        discretization_->p(iEnd,j) = discretization_->p(iEnd-1,j);
    }
    for(int i = 0; i < iEnd; i++)
    {
        // set the lower edge
        discretization_->p(i,-1) = discretization_->p(i,0);
        // set the upper edge
        discretization_->p(i,jEnd) = discretization_->p(i,jEnd-1);
    }
}

bool PressureSolver::solve()
{
    // set the initial residuum to max, it is updated in each step
    double residuum = std::numeric_limits<double>::max();
    int i = 0;

    // runs, until either the residuum is small, or it hits the max no. of iterations
    while(residuum > epsilon_ && i < maximumNumberOfIterations_) {
        // each step updates the pressure and returns the current residuum
        residuum = step();
        // then the boundary is set again and the iteration count updated
        setBoundaryValues();
        //residuum = calculateResiduum();
        i++;
    }
  
    const bool converged = residuum <= epsilon_;
    #ifndef NDEBUG
    if(converged)
    {
        std::cout << "The pressure-solver converged after \t" << i << " steps\n";
    }
    else
    {
        std::cerr << "The pressure-solver did not converged after \t" << i  << 
        " steps, the residuum left is \t" << residuum << " > " << epsilon_ << std::endl;
    }
    #endif
    // probably not used, but could be helpful
    return converged;
}

double PressureSolver::calculateResiduum()
{
    double residuum = 0;
    const double dx2 = discretization_->dx2();
    const double dy2 = discretization_->dy2();

    for(int j = 0; j < discretization_->jEnd(); j++)
    {
        for(int i = 0; i < discretization_->iEnd(); i++)
        {
            const double rhs    = discretization_->rhs(i,j);            
            const double D2px2  = discretization_->computeD2pDx2(i,j);
            const double D2py2  = discretization_->computeD2pDy2(i,j);
            const double res_ij = rhs - D2px2 - D2py2;

            residuum += res_ij*res_ij;
        }
    }
    const long nTotalCells = discretization_->nTotalCells();

    return std::sqrt(residuum / nTotalCells);
}