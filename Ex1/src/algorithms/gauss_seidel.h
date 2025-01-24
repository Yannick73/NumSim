#pragma once

#include <memory>
#include <cmath>
#include "algorithms/pressure_solver.h"
#include "discretization/discretization.h"
#include "storage/field_variable.h"

class GaussSeidel : public PressureSolver {
public:
    GaussSeidel(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations):
        PressureSolver(discretization, epsilon, maximumNumberOfIterations) {};
    double step() override;
};