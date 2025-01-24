#include "settings.h"
#include "algorithms/computation.h"
#include <cstdlib>
#include <iostream>

int main(int argc, char *argv[]) {
  
  Settings settings;
  // one arg equals no arguments given
  if(argc <= 1) {
    std::cout << "No arguments were given, use the standard settings for simulation\n";
  }
  else if(argc == 2) {
    std::string path = argv[1];
    std::cout << "Load file at \"" << path << "\" for simulation\n";
    settings.loadFromFile(path);
  }
  else {
    std::cerr << "More then one argument was given, which is unexpected, stop simulation now\n";
    return EXIT_FAILURE;
  }
  settings.printSettings();

  Computation computation;
  computation.run(settings);

  std::cout << "Simulation finished\n";

  return EXIT_SUCCESS;
}