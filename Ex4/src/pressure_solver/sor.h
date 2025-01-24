#pragma once

#include <memory>
#include <cmath>
#include <exception>
#include <sstream>
#include <mpi.h>
#include "pressure_solver/pressure_solver.h"
#include "discretization/discretization.h"
#include "discretization/partition.h"
#include "storage/field_variable.h"

class SOR : public PressureSolver
{
public:
    SOR(std::shared_ptr<PartitionShell> partition, double epsilon, int maximumNumberOfIterations, double omega);
    void step() override;
protected:
    const double omega_;
};