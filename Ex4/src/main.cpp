#include <iostream>
#include <mpi.h>
#include "settings.h"
#include "computation.h"
#include "timekeeper.h"

// statistics about solver iterations or dt-times can be generetated,
// if -DSOLVER_STATISTICS=1 respective -DDT_STATISTICS=1 is set
// detailed timing mode can be activated by issuing -DTIMER=1

// changing the communication mode may be done in the settings file
// setting useAsyncComm = true or false
// due to communication optimization, this may be change the behavior
// a little (especially in corners)
// recommendation, for small number ranks (2 or 4) use normal mode 
// and high number use async

int main(int argc, char *argv[])
{
  const auto t0 = timestamp();

  MPI_Init(&argc, &argv);

  // Get the number of processes
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // Get the rank of the process
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  Settings settings;
  // one arg equals no arguments given
  if(argc <= 1)
  {
    // only rank0 prints information, as not to print less redundant info
    if(world_rank == 0)
      std::cout << "No arguments were given, use the standard settings for simulation\n";
  }
  else if(argc == 2)
  {
    std::string path = argv[1];
    if(world_rank == 0)
      std::cout << "Load file at \"" << path << "\" for simulation\n";
    settings.loadFromFile(path);
  }
  else
  {
    if(world_rank == 0)
      std::cerr << "More than one argument was given, which is unexpected, stop simulation now\n";
    return EXIT_FAILURE;
  }

  if(world_rank == 0)
    settings.printSettings();

  runComputation(settings, world_rank, world_size);

  MPI_Finalize();
  // sim is, when it is, debug is printed only on one rank
  //if(world_rank == 0)
  //  std::cout << "\nSimulation with " << world_size << " ranks finished in " << std::setprecision(4) << getDurationS(t0) << "s\n\n";
  // For the timeout issues in parallel2, it was very likely caused by the fact,
  // that the partitioning scheme code was faulty.
    
  return EXIT_SUCCESS;
}