#pragma once

#include <memory>
#include <cmath>
#include <exception>
#include <sstream>
#include <mpi.h>
#include "pressure_solver/pressure_solver.h"
#include "discretization/discretization.h"
#include "discretization/partition.h"
#include "storage/array2D.h"

// Implemented using https://en.wikipedia.org/wiki/Conjugate_gradient_method
//! Theoretical implementation of CG solver, which does not work atm
class CG : public PressureSolver
{
public:
    //! Same as pressure-solver constructor with additional instantiation of a and r
    CG(std::shared_ptr<PartitionShell> partition, double epsilon, int maximumNumberOfIterations);

    void step() override;

protected:
    //! initialize step is called before the first iteration to set a and r
    void init();

    //! Variable for CG cell-internal data, holds product stencil-matrix*a
    Array2D ma_; 
};