#pragma once

#include <memory>
#include <cmath>
#include <cassert>
#include "algorithms/pressure_solver.h"
#include "discretization/discretization.h"
#include "storage/field_variable.h"

class SOR : public PressureSolver {
public:
    SOR(const std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations, double omega);
    double step() override;
protected:
    const double omega_;
};