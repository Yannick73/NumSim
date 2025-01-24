#pragma once

#include <memory>
#include <limits>
#include <cassert>
#include "discretization/discretization.h"

class PressureSolver {
public:
    //! Constructor with shared access to the grid-discretization and parameter
    PressureSolver(const std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations);
    
    //! Calls the solving-step in a loop, returns, whether it sufficiently converged
    bool solve();

protected:
    //! Always called at the begin of each solve-iteration
    void setBoundaryValues();

    //! The computation loop called each step
    //! Residuum is computated in step, so it may be done inline with the iteration loops (if possible)
    virtual double step() = 0;

    //! Check the residuum in each computation step
    //! calculated using the difference to the solution with euclidian norm
    double calculateResiduum();

    const std::shared_ptr<Discretization> discretization_;
    const double epsilon_;
    const int maximumNumberOfIterations_;
};