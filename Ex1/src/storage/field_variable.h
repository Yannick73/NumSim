#pragma once

#include <cmath>
#include "storage/array2D.h"

class FieldVariable : public Array2D {
public:
  FieldVariable(std::array<int, 2> size, std::array<double, 2> origin,
                std::array<double, 2> meshWidth);
  double interpolateAt(double x, double y);

private:
  const std::array<double, 2> origin_;
  const std::array<double, 2> meshWidth_;
};