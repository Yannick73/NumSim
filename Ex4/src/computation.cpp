#include "computation.h"


void runComputation(const Settings &settings, int rank, int nRanks)
{    
    std::array<double, 3> meshWidth{settings.physicalSize[0]/settings.nCells[0],
                                    settings.physicalSize[1]/settings.nCells[1],
                                    settings.physicalSize[2]/settings.nCells[2]};
    // generate partition information
    PartitionInformation pi(settings.nCells, meshWidth, rank, nRanks);
    
    std::shared_ptr<Discretization> discretization;
    if(settings.useDonorCell == true)
        discretization = std::make_shared<DonorCell>(pi, settings);
    else
        discretization = std::make_shared<CentralDifferences>(pi, settings);
    
    std::shared_ptr<PartitionShell> partition;
    partition = std::make_shared<AsyncPartition>(discretization, settings, pi);
    /*if(settings.useAsyncComm)
        partition = std::make_shared<AsyncPartition>(discretization, settings, pi);
    else
        partition = std::make_shared<Partition>(discretization, settings, pi);*/
        
    std::shared_ptr<PressureSolver> pressureSolver = newPressureSolver(partition, settings, nRanks);

    #pragma message("Output writer missing")
    #ifndef NDEBUG
    //OutputWriterText debugOut(discretization, rank);
    #endif
    //OutputWriterParaviewParallel paraviewOut(partition);

    DtCalculator dt(settings.maximumDt, settings.endTime, partition);

    double simulationTime = 0;
    // used for debugging
    int simTimestep = 0;

    if(rank == 0)
        std::cout << "\nSimulation setup, start loop\n\n";

    #ifdef TIMER
    const auto t0 = timestamp();
    #endif

    while(simulationTime < settings.endTime)
    {
         partition->setBoundaryUVW();

        std::pair<double, bool> dtValues = dt.calculate(simulationTime);
        double deltaT = dtValues.first;
        bool outputParaview = dtValues.second;
        
        simulationTime = simulationTime + deltaT;
        
        partition->calculateFGH(deltaT);

        partition->setBoundaryFGH();

        partition->calculateRHS(deltaT);

        pressureSolver->solve();

        partition->calculateUVW(deltaT);

        // necassary for correct output files
        #pragma message("Missing")
        //partition->exchangeUVW();

        /*if(outputParaview)
        {
            paraviewOut.writeFile(simulationTime);
            // diagnostic debug data
            if(rank == 0)
                std::cout << "Output for step " << simTimestep 
                    << " at sim-time " << simulationTime << " with dt " << deltaT 
                    << " last solver needed " << pressureSolver->getLastIterations() << " steps\n";
            #ifndef NDEBUG
            debugOut.writeFile(simulationTime);
            #endif
        }*/
        simTimestep++;
    }

    #ifdef DT_STATISTICS
    if(rank == 0)
        dt.printDtStatistics();
    #endif
    #ifdef SOLVER_STATISTICS
    if(rank == 0)
        pressureSolver->printIterationStats();
    #endif
    #ifdef TIMER
    double rankDtTimer = dt.timer_.getCummulatedTime();
    double rankNeighbourTimer =  partition->timer_.getCummulatedTime();
    double rankSolverTimer = pressureSolver->timer_.getCummulatedTime();
    
    double summedDtTimer;
    double summedNeighbourTimer;
    double summedSolverTimer;
    MPI_Allreduce(&rankDtTimer, &summedDtTimer, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&rankNeighbourTimer, &summedNeighbourTimer, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&rankSolverTimer, &summedSolverTimer, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    if(rank == 0)
    {
        std::stringstream timeInfoStr;
        timeInfoStr << "\nSim loop with " << nRanks << " ranks took " << std::setprecision(4) << getDurationS(t0) << "s to calculate\n";
        if(settings.useAsyncComm)
            timeInfoStr << "It used asynchronous communication\n";
        else
            timeInfoStr << "It used blocking communication\n";
        timeInfoStr << "Times are means of the cummulated times all ranks recorded:\n";
        timeInfoStr << "Dt timer: " << summedDtTimer/(double)nRanks << "s.\n";
        timeInfoStr << "Neighbour timer: " << summedNeighbourTimer/(double)nRanks << "s.\n";
        timeInfoStr << "Solver residuum timer: " << summedSolverTimer/(double)nRanks << "s.\n\n";
        std::cout << timeInfoStr.str();
    }
    #endif
}

std::shared_ptr<PressureSolver> newPressureSolver(std::shared_ptr<PartitionShell> partition, const Settings &settings, int nRanks)
{
    #pragma message("Missing pressure solvers")
    return std::make_shared<GaussSeidel>(partition, settings.epsilon, settings.maximumNumberOfIterations);
    // for a single rank, any solver works
    /*if(nRanks == 1)
    {
        if(settings.pressureSolver == "SOR")
            return std::make_shared<SOR>(partition, settings.epsilon, settings.maximumNumberOfIterations, settings.omega);
        else if(settings.pressureSolver == "GaussSeidel")
            return std::make_shared<GaussSeidel>(partition, settings.epsilon, settings.maximumNumberOfIterations);
        else if(settings.pressureSolver == "Checkerboard")
            return std::make_shared<Checkerboard>(partition, settings.epsilon, settings.maximumNumberOfIterations, settings.omega);
        else if(settings.pressureSolver == "CG")
            // Unfortunately, CG could not be implemented in time. Oh well!
            throw std::invalid_argument("CG is not fully implemented and does not work yet.\n");
        else
            throw std::invalid_argument("Invalid or non-implemented pressure solver: " + settings.pressureSolver + ", stop simulation\n.");
    }
    else
    {
        // for multiple ranks, as of now, only the checkerboards is hardcoded (although other algorithms work as well)
        return std::make_shared<Checkerboard>(partition, settings.epsilon, settings.maximumNumberOfIterations, settings.omega);
    }*/
}

std::pair<double, bool> DtCalculator::calculate(double simulationTime)
{
    #ifdef TIMER
    timer_.setT0();
    #endif

    double deltaLocal = std::min(dtConstant_, partition_->calculateVelocityDelta());
    double deltaT;
    MPI_Allreduce(&deltaLocal, &deltaT, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

    //! time to next full second (or sim end) for paraview output
    double deltaOut = nextParaviewTime_ - simulationTime;
    // this feature was overlooked, it may have caused the slowdown and the error
    bool outputParaview = false;
    if(deltaOut < deltaT)
    {
        outputParaview = true;
        //! next output is either at the exact next whole number, or the end of the sim
        nextParaviewTime_ = std::min(nextParaviewTime_+1.0, endTime_);
        deltaT = deltaOut;
    }

    #ifndef NDEBUG
    assert(deltaT > 0.);
    #endif
    #ifdef DT_STATISTICS
    dt_times_.push_back(deltaT);
    #endif
    #ifdef TIMER
    timer_.addTimeSinceT0();
    #endif

    return std::pair<double, bool>(deltaT, outputParaview);
}

void DtCalculator::printDtStatistics()
{
    #ifdef DT_STATISTICS
    double mean = 0.0;
    for(double val : dt_times_) { mean += val; }
    mean /= dt_times_.size();
    double std_dev = 0.0;
    for(double val : dt_times_) { std_dev += (val - mean) * (val - mean); }
    std_dev = std::sqrt(std_dev / dt_times_.size());
    std::cout << "dt statistics: mean = \t" << mean << " std-dev = \t" << std_dev << " #no-steps = \t" << dt_times_.size() << std::endl;
    #else
    std::cerr << "To print dt-statistics, set DDT_STATISTICS=1\n";
    #endif
}
