#pragma once

#include <memory>
#include <cmath>
#include <mpi.h>
#include "pressure_solver/pressure_solver.h"
#include "discretization/discretization.h"
#include "discretization/partition.h"
#include "storage/field_variable.h"

class GaussSeidel : public PressureSolver {
public:
    //! Use the standard pressure-solver constructor
    using PressureSolver::PressureSolver;
    void step() override;
};