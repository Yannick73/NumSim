#pragma once

#include <cmath>
#include <cassert>
#include "discretization/discretization.h"
#include "discretization/partition_information.h"

class DonorCell : public Discretization
{
public:
  //! same constructor, but using partition information
  DonorCell(PartitionInformation &pi, const Settings &settings) : 
            Discretization(pi, settings),
            alpha_(settings.alpha) { assert(alpha_ > 0); }

  double computeDuvDx(int i, int j, int k) const override;
  double computeDuvDy(int i, int j, int k) const override;

  double computeDuwDx(int i, int j, int k) const override;
  double computeDuwDz(int i, int j, int k) const override;

  double computeDvwDy(int i, int j, int k) const override;
  double computeDvwDz(int i, int j, int k) const override;

  double computeDu2Dx(int i, int j, int k) const override;
  double computeDv2Dy(int i, int j, int k) const override;
  double computeDw2Dz(int i, int j, int k) const override;

private:
  double alpha_;
};
