#include "settings.h"

#include <fstream>
#include <sstream>
#include <exception>
#include <iomanip>

void Settings::loadFromFile(std::string filename) {
  // Open file
  std::ifstream file(filename);

  // Check if file is open
  if (!file.is_open()) {
    std::stringstream str;
    str << "Could not open parameter file \"" << filename << "\".\n";
    throw std::runtime_error(str.str());
  }

  // Loop over lines of file
  std::string line;
  while (std::getline(file, line)) {
    // Remove whitespace at beginning of line
    line.erase(line.begin(),
               std::find_if(line.begin(), line.end(), [](unsigned char ch) {
                 return !std::isspace(ch);
               }));

    // If line is empty or starts with '#', skip it
    if (line.empty() || line[0] == '#')
      continue;

    // Find the position of '=' sign
    size_t equalPos = line.find('=');
    if (equalPos == std::string::npos)
      continue; // Skip if '=' not found

    // Parse parameter name
    std::string parameterName = line.substr(0, equalPos);
    // Remove trailing spaces from parameterName
    parameterName.erase(
        std::find_if(parameterName.rbegin(), parameterName.rend(),
                     [](unsigned char ch) { return !std::isspace(ch); })
            .base(),
        parameterName.end());

    // Parse value
    std::string value = line.substr(equalPos + 1);
    // Remove whitespace at beginning of value
    value.erase(value.begin(),
                std::find_if(value.begin(), value.end(), [](unsigned char ch) {
                  return !std::isspace(ch);
                }));

    // Remove comments at end of value, //! fails for # before first equal
    // sign!!
    size_t commentPos = value.find('#');
    if (commentPos != std::string::npos) {
      value = value.substr(0, commentPos);
    }

    // Remove whitespace at end of value
    value.erase(std::find_if(value.rbegin(), value.rend(),
                             [](unsigned char ch) { return !std::isspace(ch); })
                    .base(),
                value.end());

    // Parse actual value and set corresponding parameter
    setParameter(parameterName, value);
  }

  assert(physicalSize[0] > 0. && physicalSize[1] > 0.);
  assert(endTime > 0.);
  assert(re > 0.);
  assert(nCells[0] > 0 && nCells[1] > 0);
  assert(0. < tau && tau < 1.);
  assert(maximumDt > 0.);
  assert(omega > 0.); // if omega=0, p_ij would stay constant
  assert(epsilon > 0.);
  assert(maximumNumberOfIterations > 0);
}

// Helper function to set parameters
void Settings::setParameter(const std::string &name, const std::string &value) {
  if (name == "physicalSizeX") {
    sscanf(value.c_str(), "%lf", &physicalSize[0]);
  } else if (name == "physicalSizeY") {
    sscanf(value.c_str(), "%lf", &physicalSize[1]);
  } else if (name == "nCellsX") {
    sscanf(value.c_str(), "%d", &nCells[0]);
  } else if (name == "nCellsY") {
    sscanf(value.c_str(), "%d", &nCells[1]);
  } else if (name == "endTime") {
    endTime = std::stod(value);
  } else if (name == "re") {
    re = std::stod(value);
  } else if (name == "gX") {
    sscanf(value.c_str(), "%lf", &g[0]);
  } else if (name == "gY") {
    sscanf(value.c_str(), "%lf", &g[1]);
  } else if (name == "tau") {
    tau = std::stod(value);
  } else if (name == "maximumDt") {
    maximumDt = std::stod(value);
  } else if (name == "dirichletBottomX") {
    sscanf(value.c_str(), "%lf", &dirichletBcBottom[0]);
  } else if (name == "dirichletBottomY") {
    sscanf(value.c_str(), "%lf", &dirichletBcBottom[1]);
  } else if (name == "dirichletTopX") {
    sscanf(value.c_str(), "%lf", &dirichletBcTop[0]);
  } else if (name == "dirichletTopY") {
    sscanf(value.c_str(), "%lf", &dirichletBcTop[1]);
  } else if (name == "dirichletLeftX") {
    sscanf(value.c_str(), "%lf", &dirichletBcLeft[0]);
  } else if (name == "dirichletLeftY") {
    sscanf(value.c_str(), "%lf", &dirichletBcLeft[1]);
  } else if (name == "dirichletRightX") {
    sscanf(value.c_str(), "%lf", &dirichletBcRight[0]);
  } else if (name == "dirichletRightY") {
    sscanf(value.c_str(), "%lf", &dirichletBcRight[1]);
  } else if (name == "useDonorCell") {
    useDonorCell = (value == "true" || value == "1");
  } else if (name == "alpha") {
    alpha = std::stod(value);
  } else if (name == "pressureSolver") {
    pressureSolver = value;
  } else if (name == "omega") {
    omega = std::stod(value);
  } else if (name == "epsilon") {
    epsilon = std::stod(value);
  } else if (name == "maximumNumberOfIterations") {
    // double should be able to convert numbers with e, int is not
    maximumNumberOfIterations = (int)std::stod(value);
  } else if (name == "disableAdaptiveDt") {
    disableAdaptiveDt = (value == "true" || value == "1");
  } else {
    std::cout << "Unknown parameter: " << name << std::endl;
  }
}

void Settings::printSettings()

{

  std::cout << "Settings: " << std::endl

            << "  physicalSize: " << physicalSize[0] << " x " << physicalSize[1]
            << ", nCells: " << nCells[0] << " x " << nCells[1] << std::endl

            << "  endTime: " << endTime << " s, re: " << re << ", g: (" << g[0]
            << "," << g[1] << "), tau: " << tau << ", maximum dt: " << maximumDt
            << std::endl

            << "  dirichletBC: bottom: (" << dirichletBcBottom[0] << ","
            << dirichletBcBottom[1] << ")"

            << ", top: (" << dirichletBcTop[0] << "," << dirichletBcTop[1]
            << ")"

            << ", left: (" << dirichletBcLeft[0] << "," << dirichletBcLeft[1]
            << ")"

            << ", right: (" << dirichletBcRight[0] << "," << dirichletBcRight[1]
            << ")" << std::endl

            << "  useDonorCell: " << std::boolalpha << useDonorCell
            << ", alpha: " << alpha << std::endl

            << "  pressureSolver: " << pressureSolver << ", omega: " << omega
            << ", epsilon: " << epsilon
            << ", maximumNumberOfIterations: " << maximumNumberOfIterations
            << std::endl;
}
