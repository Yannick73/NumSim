#include "storage/array2D.h"

Array2D::Array2D(std::array<int, 2> size) : size_(size), name_("unnamed")
{
  assert(size[0] > 0);
  assert(size[1] > 0);
  // allocate data, initialize to 0
  data_.resize(size_[0] * size_[1], 0.0);
}
Array2D::Array2D(std::array<int, 2> size, std::string name) : size_(size), name_(name)
{
  assert(size[0] > 0);
  assert(size[1] > 0);
  // allocate data, initialize to 0
  data_.resize(size_[0] * size_[1], 0.0);
}