#pragma once

#include <cmath>
#include <exception>
#include <sstream>
#include <mpi.h>
#include "pressure_solver/pressure_solver.h"
#include "discretization/discretization.h"
#include "discretization/partition_shell.h"
#include "storage/field_variable.h"

class Checkerboard : public PressureSolver {
public:
    Checkerboard(std::shared_ptr<PartitionShell> partition, double epsilon, int maximumNumberOfIterations, double omega);
    void step() override;
protected:
    const double omega_;
};