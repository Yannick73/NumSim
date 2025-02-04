#include "storage/array3D.h"

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

inline double &Array3D::operator()(int i, int j, int k)
{
  const std::size_t index = compute_index(i, j, k);
  // make assertion conditional on DEBUG mode to optimize further
  #ifndef NDEBUG
  if(i < 0 || i >= size_[0] || j < 0 || j >= size_[1] || k < 0 || k >= size_[2])
  {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::stringstream str;
    str << "Out-of-bound access on " << name_ << "(i,j,k): (" << i << ',' << j << ',' << k
        << "), size: (" << size_[0] << ',' << size_[1] << ',' << size_[2] << ") in R:" << rank << "\n";
    throw std::out_of_range(str.str());
  }
  #endif

  return data_[index];
}

inline double Array3D::operator()(int i, int j, int k) const
{
  //std::size_t index = compute_index(i, j, k);
  const std::size_t index = compute_index(i, j, k);
  // assert that indices are in range
  #ifndef NDEBUG
  if(i < 0 || i >= size_[0] || j < 0 || j >= size_[1] || k < 0 || k >= size_[2])
  {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::stringstream str;
    str << "Out-of-bound access on " << name_ << "(i,j,k): (" << i << ',' << j << ',' << k
        << "), size: (" << size_[0] << ',' << size_[1] << ',' << size_[2] << ") in R:" << rank << "\n";
    throw std::out_of_range(str.str());
  }
  #endif

  return data_[index];
}