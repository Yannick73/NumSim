#include "storage/array2D.h"


// size[0] is the number of rows (length y), size[1] is the number of columns
// (length x)
Array2D::Array2D(std::array<int, 2> size) : size_(size)
{
  assert(size[0] > 0);
  assert(size[1] > 0);
  // allocate data, initialize to 0
  data_.resize(size_[0] * size_[1], 0.0);
}