#pragma once
#include <algorithm>
#include <array>
#include <cctype>
#include <iostream>
#include <string>
#include <cassert>

/** All settings that parametrize a simulation run.

 */
struct Settings 
{
    std::array<int, 3> nCells{20, 20, 20};          //< number of cells in x, y and z direction
    std::array<double, 3> physicalSize{2, 2, 2}; //< physical size of the domain
    double re = 1000;                   //< reynolds number
    double endTime = 10.0;              //< end time of the simulation
    double tau = 0.5;                   //< safety factor for time step width
    double maximumDt = 0.1;             //< maximum time step width
    double minimumDt = 0.0;             //< minimum time step width for stability
    double outputDt = 1.0;              //< time interval for output
    bool disableAdaptiveDt = false;     //< if set true, dt will not be dependent on the fluid velocities
    std::array<double, 3> g{0., 0., 0.};    //< external forces
    bool useDonorCell = false; //< if the donor cell scheme schould be used
    double alpha = 0.5;        //< factor for donor-cell scheme

    std::string boundaryTop = "Inflow"; //< boundary condition at the top of the domain
    std::string boundaryBottom = "Inflow"; //< boundary condition at the bottom of the domain
    std::string boundaryLeft = "Inflow"; //< boundary condition at the left of the domain
    std::string boundaryRight = "Inflow"; //< boundary condition at the right of the domain
    std::string boundaryFront = "Inflow"; //< boundary condition at the front of the domain
    std::string boundaryHind = "Inflow"; //< boundary condition at the hind of the domain

    double pressureTop = 0.; //< prescribed pressure at top of domain
    double pressureBottom = 0.; //< prescribed pressure at bottom of domain
    double pressureLeft = 0.; //< prescribed pressure at left of domain
    double pressureRight = 0.; //< prescribed pressure at right of domain
    double pressureFront = 0.; //< prescribed pressure at front of domain
    double pressureHind = 0.; //< prescribed pressure at hind of domain

    std::array<double, 3>
        dirichletBcBottom{0., 0., 0.}; //< prescribed values of u,v at bottom of domain
    std::array<double, 3>
        dirichletBcTop{1., 0., 0.}; //< prescribed values of u,v at top of domain
    std::array<double, 3>
        dirichletBcLeft{0., 0., 0.}; //< prescribed values of u,v at left of domain
    std::array<double, 3>
        dirichletBcRight{0., 0., 0.}; //< prescribed values of u,v at right of domain
    std::array<double, 3>
        dirichletBcFront{0., 0., 0.};
    std::array<double, 3>
        dirichletBcHind{0., 0., 0.};
    std::string pressureSolver =
        "Checkerboard";          //< which pressure solver to use, "GaussSeidel" or "SOR"
    double omega = 1.6; //< overrelaxation factor
    double epsilon = 1e-5; //< tolerance for the residual in the pressure solver
    int maximumNumberOfIterations =
        1e4; //< maximum number of iterations in the solver
    bool useAsyncComm = true; //< If asynchronous MPI communication is to be used

    //! parse a text file with settings, each line contains "<parameterName> =
    //! <value>"
    void loadFromFile(std::string filename);

    //! output all settings to console
    void printSettings();

private:
    void setParameter(const std::string &name, const std::string &value);
};
