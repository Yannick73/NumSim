#pragma once

#include <iostream>
#include <fstream>
#include <exception>
#include <sstream>
#include <array>
#include <vector>
#include <algorithm>

// helper function to convert an array to a human readable string
template <typename T, std::size_t N>
std::string array2str(const std::array<T, N> &input);
// returns the same string without whitespaces and in lowercase
std::string strclean(const std::string input);


/** All settings that parametrize a simulation run.
 *  Header given by exercise.
 */
struct Settings
{

public:

  std::array<int,2> nCells{0, 0};          //< number of cells in x and y direction

  std::array<double,2> physicalSize{0., 0.}; //< physical size of the domain

  double re = 1000;                  //< reynolds number

  double endTime = 10.0;             //< end time of the simulation

  double tau = 0.5;                  //< safety factor for time step width

  double maximumDt = 0.1;            //< maximum time step width


  std::array<double,2> g{0., 0.};    //< external forces


  bool useDonorCell = false;         //< if the donor cell scheme schould be used

  double alpha = 0.5;                //< factor for donor-cell scheme


  std::array<double,2> dirichletBottom{0., 0.};  //< prescribed values of u,v at bottom of domain

  std::array<double,2> dirichletTop{0., 0.};     //< prescribed values of u,v at top of domain

  std::array<double,2> dirichletLeft{0., 0.};    //< prescribed values of u,v at left of domain

  std::array<double,2> dirichletRight{0., 0.};   //< prescribed values of u,v at right of domain


  std::string pressureSolver = "SOR";      //< which pressure solver to use, "GaussSeidel" or "SOR"

  double omega = 1.0;                //< overrelaxation factor

  double epsilon = 1e-5;             //< tolerance for the residual in the pressure solver

  int maximumNumberOfIterations = 1e5;    //< maximum number of iterations in the solver


  //! parse a text file with settings, each line contains "<parameterName> = <value>"

  void loadFromFile(std::string filename);


  //! output all settings to console

  void printSettings();

private:
  // each line of the file, that contains data
  // removes comments starting with '#', whitespaces and casts everything to lower case
  std::vector<std::string> file_lines;

  // given the name searches the file_lines for anything starting with the name and returns the value
  std::string extract_string(const std::string name);

  /**
   * \brief   Set the attribute to the value of the name found in the input-file
   * \param   *attribute is the double value to be set
   * \param   name defines the value in the text file to be found
   * \param   input_file define the file to be searched for the name
   * */
  void load_double(double *attribute, const std::string name);

  //! Same with integer
  void load_integer(int *attribute, const std::string name);

  //! Boolean same as double, but the strings "1", "true" or "yes" are interpreted as true
  void load_boolean(bool *attribute, const std::string name);

  //! Directly saves the result of extract_string in the *attribute
  void load_string(std::string *attribute, const std::string name);

};
