#pragma once

#include <array>
#include <exception>
#include <limits>
#include "settings.h"
#include "output_writer/output_writer.h"
#include "output_writer/output_writer_paraview.h"
#include "output_writer/output_writer_text.h"
#include "discretization/discretization.h"
#include "discretization/central_differences.h"
#include "discretization/donor_cell.h"
#include "algorithms/pressure_solver.h"
#include "algorithms/gauss_seidel.h"
#include "algorithms/sor.h"

class Computation {
    private:
        //PressureSolver* pressureSolver;
        std::shared_ptr<Discretization> discretization;
        std::shared_ptr<PressureSolver> pressureSolver;

    public:
        void run(Settings settings);
};