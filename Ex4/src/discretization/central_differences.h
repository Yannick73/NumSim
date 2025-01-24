#pragma once

#include "discretization/discretization.h"

class CentralDifferences : public Discretization
{
public:
  //! re-use the Discretization constructors
  using Discretization::Discretization;

  double computeDuvDx(int i, int j) const override;
  double computeDuvDy(int i, int j) const override;
  double computeDu2Dx(int i, int j) const override;
  double computeDv2Dy(int i, int j) const override;
};
