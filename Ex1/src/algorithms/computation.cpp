#include "computation.h"
#define DT_STATISTICS

void Computation::run(Settings settings)
{    
    std::array<double, 2> meshWidth{settings.physicalSize[0]/settings.nCells[0], settings.physicalSize[1]/settings.nCells[1]};
    
    if(settings.useDonorCell == true)
    {
        discretization = std::make_shared<DonorCell>(settings.nCells, meshWidth, settings);
    }
    else
    {
        discretization = std::make_shared<CentralDifferences>(settings.nCells, meshWidth, settings);
    }

    if(settings.pressureSolver == "SOR")
    {
        pressureSolver = std::make_shared<SOR>(discretization, settings.epsilon, settings.maximumNumberOfIterations, settings.omega);
    }
    else if(settings.pressureSolver == "GaussSeidel")
    {
        pressureSolver = std::make_shared<GaussSeidel>(discretization, settings.epsilon, settings.maximumNumberOfIterations);
    }
    else
    {
        throw std::runtime_error("Invalid or non-implemented pressure solver: " + settings.pressureSolver + ", stop simulation\n.");
    }


    #ifndef NDEBUG
    OutputWriterText debugOut(discretization);
    #endif
    OutputWriterParaview paraviewOut(discretization);

    #ifdef DT_STATISTICS
    std::vector<double> dt_times;
    #endif

    discretization->setBoundary(settings);

    // define for more convenient access
    const double dx2 = discretization->dx2();
    const double dy2 = discretization->dy2();
    const double reynoldsDelta = settings.tau * (settings.re/2) * (dx2*dy2) / (dx2+dy2);
    double velocityDelta = std::numeric_limits<double>::max();

    double simulationTime = 0;
    // used for debugging
    int simTimestep = 0;

    while(simulationTime < settings.endTime)
    {
        discretization->setBoundaryValues();

        if(!settings.disableAdaptiveDt)
        {
            velocityDelta = discretization->calculateVelocityDelta();
        }
        double deltaStop = settings.endTime - simulationTime;
        double deltaT = std::min(reynoldsDelta, std::min(settings.maximumDt, velocityDelta));
        // extra if for debugging
        if(deltaStop < deltaT)
        {
            #ifndef NDEBUG
            std::cout << "\nExecuting the last sim time step!\n\n";
            #endif
            deltaT = deltaStop;
        }
        if(deltaT <= 0.0)
        {
            deltaT = std::numeric_limits<float>::min();
            std::cerr << "Warning: deltaT was smaller then 0.0! Set it to float-eps: " << deltaT << std::endl;
        }
        #ifdef DT_STATISTICS
        dt_times.push_back(deltaT);
        #endif
        simulationTime = simulationTime + deltaT;
        
        discretization->calculateFG_RHS(deltaT);

        pressureSolver->solve();

        discretization->calculateUV(deltaT);

        #ifndef NDEBUG
        std::cout << "Write output files for step \t" << simTimestep << " at sim-time \t" << simulationTime << " with dt \t" << deltaT << std::endl;
        debugOut.writeFile(simulationTime);
        #endif
        paraviewOut.writeFile(simulationTime);
        simTimestep++;
    }

    #ifdef DT_STATISTICS
    double mean = 0.0;
    for(double val : dt_times) { mean += val; }
    mean /= dt_times.size();
    double std_dev = 0.0;
    for(double val : dt_times) { std_dev += (val - mean) * (val - mean); }
    std_dev = std::sqrt(std_dev / dt_times.size());
    std::cout << "dt statistics: mean = \t" << mean << " std-dev = \t" << std_dev << " #no-steps = \t" << dt_times.size() << std::endl;
    #endif
}
