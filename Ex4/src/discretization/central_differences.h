#pragma once

#include "discretization/discretization.h"

class CentralDifferences : public Discretization
{
public:
  //! re-use the Discretization constructors
  using Discretization::Discretization;

  double computeDuvDx(int i, int j, int k) const override;
  double computeDuvDy(int i, int j, int k) const override;

  double computeDuwDx(int i, int j, int k) const override;
  double computeDuwDz(int i, int j, int k) const override;

  double computeDvwDy(int i, int j, int k) const override;
  double computeDvwDz(int i, int j, int k) const override;

  double computeDu2Dx(int i, int j, int k) const override;
  double computeDv2Dy(int i, int j, int k) const override;
  double computeDw2Dz(int i, int j, int k) const override;
};
