#pragma once

#include <cmath>
#include <cassert>
#include "discretization/discretization.h"

class DonorCell : public Discretization

{

public:
  //! use the constructor of the base class
  DonorCell(std::array<int, 2> nCells, std::array<double, 2> meshWidth,
            Settings settings);

  //! compute the 1st derivative ∂ u^2 / ∂x
  double computeDu2Dx(int i, int j) const override;

  //! compute the 1st derivative ∂ v^2 / ∂x
  double computeDv2Dy(int i, int j) const override;

  //! compute the 1st derivative ∂ (uv) / ∂x
  double computeDuvDx(int i, int j) const override;

  //! compute the 1st derivative ∂ (uv) / ∂y
  double computeDuvDy(int i, int j) const override;

private:
  double alpha_;
};
