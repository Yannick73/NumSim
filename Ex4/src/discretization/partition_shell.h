#pragma once

#include <cassert>
#include <vector>
#include <set>
#include "settings.h"
#include "timekeeper.h"
#include "discretization/discretization.h"
#include "discretization/partition_information.h"
#include "boundary/boundary.h"
#include "boundary/dirichlet.h"

//! Class to encapsulate the discretization and its boundaries
class PartitionShell
{
public:
    PartitionShell(const std::shared_ptr<Discretization> discretization,
                   const PartitionInformation &pi) :
                   discretization_(discretization), pi_(pi) { }

    //! set applies all boundary conditions on fix boundary an on neighbours
    virtual void setBoundaryUVW() = 0;
    virtual void setBoundaryFGH() = 0;
    virtual void setBoundaryP() = 0;
    //! only exchange the pressure values with the respective neighbours w/o setting dirichlet
    //#pragma message("Missing exchange")
    virtual void exchangeP() = 0;
    //! also only used in pressure solver, but specifically only for CG
    //virtual void exchangeRA() = 0;
    //! used before paraview output
    virtual void exchangeUVW() = 0;

    //! getter for the discretization, to make access in pressure-solver simpler
    inline std::shared_ptr<Discretization> getDiscretization() const { return discretization_; };

    //! wrapper for the relevant public discretization methods
    inline double calculateVelocityDelta() const { return discretization_->calculateVelocityDelta(); }
    inline void calculateFGH(double deltaT) const { discretization_->calculateFGH(deltaT); }
    inline void calculateRHS(double deltaT) const { discretization_->calculateRHS(deltaT); }
    inline void calculateUVW(double deltaT) const { discretization_->calculateUVW(deltaT); }
    inline double calculateReynoldsDelta() const { return discretization_->calculateReynoldsDelta(); }

    inline void fixDirichletBoundary()
    {
        for(std::shared_ptr<Dirichlet> fixBoundary : fixBoundaries_)
            fixBoundary->setUVW();
    }

    //! information relevant to the partition
    const PartitionInformation &pi_;

    #ifdef TIMER
    //! neighbour communication timer
    Timekeeper timer_;
    #endif

protected:
    //! Internal discretization pointer with getter above, so it can be manipulated, but not re-allocated
    const std::shared_ptr<Discretization> discretization_;
    //! Every partition will always have some dirichlet-boundary condition to set
    std::vector<std::shared_ptr<Dirichlet>> fixBoundaries_;
};