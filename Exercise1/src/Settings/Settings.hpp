#pragma once


#include <iostream>
#include <ifstream>
#include <exception>
#include <sstream>
#include <array>
#include <cstdlib>


/** All settings that parametrize a simulation run.
 *  Header given by exercise.
 */

// helper function to convert an array to a human readable string
template <typename T, std::size_t N>
std::string array2str(const std::array<T, N> &input);

struct Settings

{

  std::array<int,2> nCells;          //< number of cells in x and y direction

  std::array<double,2> physicalSize; //< physical size of the domain

  double re = 1000;                  //< reynolds number

  double endTime = 10.0;             //< end time of the simulation

  double tau = 0.5;                  //< safety factor for time step width

  double maximumDt = 0.1;            //< maximum time step width


  std::array<double,2> g{0., 0.};    //< external forces


  bool useDonorCell = false;         //< if the donor cell scheme schould be used

  double alpha = 0.5;                //< factor for donor-cell scheme


  std::array<double,2> dirichletBottom;  //< prescribed values of u,v at bottom of domain

  std::array<double,2> dirichletTop;     //< prescribed values of u,v at top of domain

  std::array<double,2> dirichletLeft;    //< prescribed values of u,v at left of domain

  std::array<double,2> dirichletRight;   //< prescribed values of u,v at right of domain


  std::string pressureSolver = "SOR";      //< which pressure solver to use, "GaussSeidel" or "SOR"

  double omega = 1.0;                //< overrelaxation factor

  double epsilon = 1e-5;             //< tolerance for the residual in the pressure solver

  int maximumNumberOfIterations = 1e5;    //< maximum number of iterations in the solver


  //! parse a text file with settings, each line contains "<parameterName> = <value>"

  void loadFromFile(std::string filename);


  //! output all settings to console

  void printSettings();

};
