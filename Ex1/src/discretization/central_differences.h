#pragma once

#include <array>
#include "discretization.h"
#include "settings.h"

class CentralDifferences : public Discretization {
public:
  //! construct the object with given number of cells in x and y direction

  CentralDifferences(std::array<int, 2> nCells,
                     std::array<double, 2> meshWidth,
                     Settings settings);

  double computeDuvDx(int i, int j) const override;
  double computeDuvDy(int i, int j) const override;
  double computeDu2Dx(int i, int j) const override;
  double computeDv2Dy(int i, int j) const override;
  /*double computeD2uDx2(int i, int j) const;
  double computeD2uDy2(int i, int j) const;
  double computeD2vDx2(int i, int j) const;
  double computeD2vDy2(int i, int j) const;*/
};
