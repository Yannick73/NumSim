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
    std::array<int, 2> nCells{20, 20};          //< number of cells in x and y direction
    std::array<double, 2> physicalSize{2, 2}; //< physical size of the domain
    double re = 1000;                   //< reynolds number
    double endTime = 10.0;              //< end time of the simulation
    double tau = 0.5;                   //< safety factor for time step width
    double maximumDt = 0.1;             //< maximum time step width
    bool disableAdaptiveDt = false;     //< if set true, dt will not be dependent on the fluid velocities
    std::array<double, 2> g{0., 0.};    //< external forces
    bool useDonorCell = true; //< if the donor cell scheme schould be used
    double alpha = 0.5;        //< factor for donor-cell scheme
    std::array<double, 2>
        dirichletBcBottom{0., 0.}; //< prescribed values of u,v at bottom of domain
    std::array<double, 2>
        dirichletBcTop{1., 0.}; //< prescribed values of u,v at top of domain
    std::array<double, 2>
        dirichletBcLeft{0., 0.}; //< prescribed values of u,v at left of domain
    std::array<double, 2>
        dirichletBcRight{0., 0.}; //< prescribed values of u,v at right of domain
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
