#pragma once

#include <array>
#include <memory>
#include <vector>
#include <algorithm>
#include <exception>
#include <limits>
#include <utility>
#include <mpi.h>
#include "settings.h"
#include "timekeeper.h"
#include "output_writer/output_writer_paraview_parallel.h"
#include "discretization/discretization.h"
#include "discretization/partition_shell.h"
#include "discretization/async_partition.h"
#include "discretization/central_differences.h"
#include "discretization/donor_cell.h"
#include "boundary/boundary.h"
#include "boundary/dirichlet.h"
#include "pressure_solver/pressure_solver.h"
#include "pressure_solver/gauss_seidel.h"
#include "pressure_solver/sor.h"
#include "pressure_solver/checkerboard.h"

//! only a small class intended to facilitate calculation and statistics about dt
class DtCalculator
{
public:
    //! sets all calculator constants
    DtCalculator(const Settings &settings, std::shared_ptr<PartitionShell> partition) :
                 dtConstant_(std::min(partition->calculateReynoldsDelta(), settings.maximumDt)),
                 partition_(partition), rank_(partition->pi_.ownRankNo()),
                 endTime_(settings.endTime), dtMin_(settings.minimumDt), dtOut_(settings.outputDt) { }

    //! calculates the next dt and defines, whether the next time step is to be output
    std::pair<double, bool> calculate(double simulationTime);

    //! prints statistics for dt when DT_STATISTICS macro is defined
    void printDtStatistics();

    #ifdef TIMER
    //! timer
    Timekeeper timer_;
    #endif

private:
    //! dtConstant is the smaller of maxDt and reynoldsDt, which does not change throughout the sim
    const double dtConstant_;

    //! minimal deltaT for stability reasons, used only, if set to > 0.0
    const double dtMin_;

    //! output writer interval
    const double dtOut_;

    //! partition reference, to calculate the velocity delta
    const std::shared_ptr<PartitionShell> partition_;

    //! next timepoint to generate a paraview output file
    double nextParaviewTime_ = dtOut_;

    //! end of the simulation
    const double endTime_;

    //! rank of the partition for debugging
    const int rank_;

    #ifdef DT_STATISTICS
    std::vector<double> dt_times_;
    #endif
};

//! This may be optimized into a class in future, but for now define it as simple function
void runComputation(const Settings &settings, int rank, int nRanks);

//! Generates and returns the appropriate pressure-solver
std::shared_ptr<PressureSolver> newPressureSolver(std::shared_ptr<PartitionShell> partition, const Settings &settings, int nRanks);