#include "pressure_solver/pressure_solver.h"

PressureSolver::PressureSolver(std::shared_ptr<PartitionShell> partition, 
                               double epsilon, int maximumNumberOfIterations) : 
                               partition_(partition),
                               discretization_(partition->getDiscretization()),
                               epsilon2_(epsilon*epsilon),
                               maximumNumberOfIterations_(maximumNumberOfIterations),
                               rank_(partition->pi_.ownRankNo())
{
    assert(partition->getDiscretization() != nullptr);
    if(epsilon <= 0)
        throw std::out_of_range("epsilon must be strictly positive\n");
    if(maximumNumberOfIterations <= 0)
        throw std::out_of_range("maximumNumberOfIterations must be strictly positive!");
    #ifndef NDEBUG
    // shared-ptr takes care of later object destruction, when pressure-solver is destructed
    pressureDebug = std::make_shared<OutputWriterText>(discretization_, rank_);
    #endif
}

bool PressureSolver::solve()
{
    // set the initial residuum to max, it is updated in each step
    double residuum2;
    iteration_ = 0;

    // runs, until either the residuum is small, or it hits the max no. of iterations
    do
    {
        partition_->setBoundaryP();
        // each step updates the pressure and returns the current residuum
        step();
        residuum2 = calculateResiduum2();
        //std::cout << "Res2: " << residuum2 << std::endl;
    } while (residuum2 > epsilon2_ && ++iteration_ < maximumNumberOfIterations_);
    #ifdef SOLVER_STATISTICS
    solverStepStatistics_.push_back(iteration_);
    #endif
    partition_->setBoundaryP();

    const bool converged = residuum2 <= epsilon2_;
    #ifndef NDEBUG
    if(rank_ == 0)
    {
        if(converged)
            std::cout << "The pressure-solver converged after \t" << iteration_ << " steps\n";
        else
        {
            std::cerr << "The pressure-solver did not converged after \t" << iteration_  << 
            " steps, the residuum left is \t" << residuum2 << " > " << epsilon2_ << std::endl;
            throw std::runtime_error("Pressure-Solver did not converge!\n");
        }
    }
    #endif
    // probably not used, but could be helpful
    return converged;
}

double PressureSolver::calculateResiduum2()
{
    #ifdef TIMER
    timer_.setT0();
    #endif

    double partitionResiduum = 0;
    const double dx2 = discretization_->dx2();
    const double dy2 = discretization_->dy2();

    #pragma omp simd collapse(2) reduction(+:partitionResiduum)
    for(int j = 0; j < discretization_->pjN(); j++)
    {
        for(int i = 0; i < discretization_->piN(); i++)
        {
            const double rhs    = discretization_->rhs(i,j);            
            const double D2px2  = discretization_->computeD2pDx2(i,j);
            const double D2py2  = discretization_->computeD2pDy2(i,j);
            const double res_ij = rhs - D2px2 - D2py2;

            partitionResiduum += res_ij*res_ij;
        }
    }

    double overallResiduum;
    MPI_Allreduce(&partitionResiduum, &overallResiduum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    #ifdef TIMER
    timer_.addTimeSinceT0();
    #endif

    return overallResiduum / partition_->pi_.totalNoOfCellsGlobal();
}

void PressureSolver::printIterationStats()
{
    #ifdef SOLVER_STATISTICS
    long sum = 0;
    for(double val : solverStepStatistics_) { sum += val; }
    double mean = ((double)sum) / solverStepStatistics_.size();
    double std_dev = 0.0;
    for(double val : solverStepStatistics_) { std_dev += (val - mean) * (val - mean); }
    std_dev = std::sqrt(std_dev / solverStepStatistics_.size());
    std::cout << "Solver " << typeid(*this).name() << " stats: mean = " << mean 
              << ", std-dev = " << std_dev << ", overall steps = " << sum << std::endl;
    #else
    std::cerr << "To print pressure solver statistics, set DSOLVER_STATISTICS=1\n";
    #endif
}