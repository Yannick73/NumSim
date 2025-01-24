#pragma once

#include <memory>
#include <mpi.h>
#include <vtkSmartPointer.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkImageData.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include "output_writer/output_writer.h"
#include "discretization/partition.h"
#include "discretization/partition_information.h"
#include "discretization/discretization.h"
#include "storage/field_variable.h"

/** Write *.vti files that can be viewed with ParaView.
 *  The mesh that can be visualized in ParaView corresponds to the mesh of the computational domain.
 *  All values are given for the nodes of the mesh, i.e., the corners of each cell.
 *  This means, values will be interpolated because the values are stored at positions given by the staggered grid.
 */
class OutputWriterParaviewParallel : 
  public OutputWriter
{
public:
  //! constructor
  OutputWriterParaviewParallel(const std::shared_ptr<PartitionShell> partition);

  //! write current velocities to file, filename is output_<count>.vti
  void writeFile(double currentTime);

private:

  //! gather u,v and p values from all ranks to rank 0 and store them in the global field variables
  void gatherData();

  const PartitionInformation &pi_;

  vtkSmartPointer<vtkXMLImageDataWriter> vtkWriter_;   //< vtk writer to write ImageData

  std::array<int,2> nCellsGlobal_;   //< global number of cells
  std::array<int,2> nPointsGlobal_;  //< global number of points

  FieldVariable uLocal_;    // field variable for u with global size, contains only the local values, other entries are 0
  FieldVariable vLocal_;    // field variable for v with global size, contains only the local values, other entries are 0
  FieldVariable pLocal_;    // field variable for p with global size, contains only the local values, other entries are 0

  FieldVariable uGlobal_;    // on rank 0: field variable for u that gathers values from all ranks, on other ranks: nullptr
  FieldVariable vGlobal_;    // on rank 0: field variable for v that gathers values from all ranks, on other ranks: nullptr
  FieldVariable pGlobal_;    // on rank 0: field variable for p that gathers values from all ranks, on other ranks: nullptr
};
