#pragma once

#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "output_writer/output_writer.h"

/** Write *.txt files that are useful for debugging.
 *  All values are written to the file as they are stored in the field variables,
 *  no interpolation takes place.
 */
class OutputWriterText : 
  public OutputWriter
{
public:
  OutputWriterText(std::shared_ptr<Discretization> discretization, int rank) :
  OutputWriter(discretization), rank_(rank) {}
  //! use constructor of base class
  //using OutputWriter::OutputWriter;

  //! write current velocities to file, filename is output_<count>.txt
  void writeFile(double currentTime);

  //! write only current values of pressure to file, filename is pressure_<count>.txt
  void writePressureFile(int simStep, int pressureStep, int numMaxPressureSteps, char substep = '\0');

private:
  const int rank_;

};
