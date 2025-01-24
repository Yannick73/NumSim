#pragma once

#include <memory>
#include <limits>
#include <cassert>
#include <cmath>
#include <exception>
#include <vector>
#include <mpi.h>
#include "timekeeper.h"
#include "discretization/partition_shell.h"
#include "discretization/discretization.h"
#include "output_writer/output_writer_text.h"

class PressureSolver {
public:
    //! Constructor with shared access to the grid-discretization and parameter
    PressureSolver(std::shared_ptr<PartitionShell> partition, double epsilon, int maximumNumberOfIterations);
    
    //! Calls the solving-step in a loop, returns, whether it sufficiently converged
    bool solve();

    //! Output a statistic over how many iterations the solver needed, only defined for !NDEBUG
    void printIterationStats();

    //! returns the last number of iterations, the solver used
    inline int getLastIterations() const { return iteration_; }

    #ifdef TIMER
    //! timer
    Timekeeper timer_;
    #endif

protected:
    //! current iteration
    int iteration_ = 0;
    
    //! The computation loop called each step
    //! Residuum is computated in step, so it may be done inline with the iteration loops (if possible)
    virtual void step() = 0;

    //! Check the residuum in each computation step
    //! calculated using the difference to the solution with euclidian norm
    //! returns the squared result so it may be added over the whole domain and then sqrt
    double calculateResiduum2();

    std::shared_ptr<PartitionShell> partition_;
    const std::shared_ptr<Discretization> discretization_;
    const double epsilon2_;
    const int maximumNumberOfIterations_;
    const int rank_;

    #ifndef NDEBUG
    std::shared_ptr<OutputWriterText> pressureDebug;
    #endif
    #ifdef SOLVER_STATISTICS
    //! overall step of simulation
    std::vector<int> solverStepStatistics_;
    #endif
};