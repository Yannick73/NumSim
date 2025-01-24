#pragma once

#include <cmath>
#include <iostream>
#include "storage/array2D.h"

//! storage class holding and managing the data variable
class FieldVariable : public Array2D 
{
public:
  FieldVariable(std::array<int, 2> size, std::array<double, 2> origin,
                std::array<double, 2> meshWidth, std::string name);

  //! sets all the data to zero using the fill command
  inline void setToZero() { std::fill(data_.begin(), data_.end(), 0.0); };

  //! gets the underlying buffer definition
  inline double* data() { return data_.data(); };

  //! simple linear interpolation schemata based on the assumption,
  //! that the corresponding other dimension is already very close
  double horizontalInterpolation(double x, double y);
  double verticalInterpolation(double x, double y);

  //! more advanced bilinear interpolation schema, if the points are anywhere
  double bilinearInterpolation(double x, double y);

  inline const std::array<double, 2> getOrigin() { return origin_; }

private:
  //! physical dimensions
  const std::array<double, 2> origin_;
  const std::array<double, 2> meshWidth_;
};