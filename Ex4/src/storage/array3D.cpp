#include "storage/array3D.h"

// size[0] is the number of rows (length y), size[1] is the number of columns
// (length x)
Array3D::Array3D(std::array<int, 3> size) : 
size_(size), name_("unnamed"), size0Xsize1_(size[0]*size[1])
{
  assert(size[0] > 0);
  assert(size[1] > 0);
  assert(size[2] > 0);
  std::size_t total_size = 
    ((std::size_t)size_[0]) * ((std::size_t)size_[1]) * ((std::size_t)size_[2]);
  // allocate data, initialize to 0
  data_.resize(total_size, 0.0);
}
Array3D::Array3D(std::array<int, 3> size, std::string name) : 
size_(size), name_(name), size0Xsize1_(size[0]*size[1])
{
  assert(size[0] > 0);
  assert(size[1] > 0);
  assert(size[2] > 0);
  std::size_t total_size = 
    ((std::size_t)size_[0]) * ((std::size_t)size_[1]) * ((std::size_t)size_[2]);
  // allocate data, initialize to 0
  data_.resize(total_size, 0.0);
}