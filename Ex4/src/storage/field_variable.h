#pragma once

#include <cmath>
#include <iostream>
#include "storage/array3D.h"

//! storage class holding and managing the data variable
class FieldVariable : public Array3D 
{
public:
  FieldVariable(std::array<int, 3> size, std::array<double, 3> origin,
                std::array<double, 3> meshWidth, std::string name);

  //! sets all the data to zero using the fill command
  inline void setToZero() { std::fill(data_.begin(), data_.end(), 0.0); };

  //! gets the underlying buffer definition
  inline double* data() { return data_.data(); };

  //! simple linear interpolation schemata based on the assumption,
  //! that the corresponding other dimension is already very close
  //double horizontalInterpolation(double x, double y);
  //double verticalInterpolation(double x, double y);

  //! more advanced bilinear interpolation schema, if the points are anywhere
  //double bilinearInterpolation(double x, double y);
  // for u x lies on the correct x grid, but y and z does not
  double yzInterpolation (double x, double y, double z);
  double xzInterpolation (double x, double y, double z);
  double xyInterpolation (double x, double y, double z);
  double midInterpolation(double x, double y, double z);

  inline const std::array<double, 3> getOrigin() { return origin_; }

private:
  //! physical dimensions
  const std::array<double, 3> origin_;
  const std::array<double, 3> meshWidth_;
};