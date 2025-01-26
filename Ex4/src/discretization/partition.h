#pragma once

#include <cassert>
#include <vector>
#include <set>
#include "settings.h"
#include "timekeeper.h"
#include "discretization/discretization.h"
#include "discretization/partition_shell.h"
#include "discretization/partition_information.h"
#include "boundary/boundary.h"
#include "boundary/dirichlet.h"
#include "boundary/neighbour_boundary.h"

//! Class to encapsulate the discretization and its boundaries
class Partition : public PartitionShell
{
public:
    Partition(const std::shared_ptr<Discretization> discretization,
              const Settings &settings,
              const PartitionInformation &pi);

    //! set applies all boundary conditions on fix boundary an on neighbours
    virtual void setBoundaryUV();
    virtual void setBoundaryFG();
    virtual void setBoundaryP();
    //! only exchange the pressure values with the respective neighbours w/o setting dirichlet
    virtual void exchangeP();
    //! also only used in pressure solver, but specifically only for CG
    virtual void exchangeRA();
    //! used before paraview output
    virtual void exchangeUV();

    //! getter for the discretization, to make access in pressure-solver simpler
    inline std::shared_ptr<Discretization> getDiscretization() const { return discretization_; };

    //! wrapper for the relevant public discretization methods
    inline double calculateVelocityDelta() const { return discretization_->calculateVelocityDelta(); }
    inline void calculateFG(double deltaT) const { discretization_->calculateFGH(deltaT); }
    inline void calculateRHS(double deltaT) const { discretization_->calculateRHS(deltaT); }
    inline void calculateUV(double deltaT) const { discretization_->calculateUVW(deltaT); }
    inline double calculateReynoldsDelta() const { return discretization_->calculateReynoldsDelta(); }

protected:    
    //! Keep a secondary access here, to only exchange pressure values in exchangeP method
    std::vector<std::shared_ptr<NeighbourBoundary>> neighbours_;
};