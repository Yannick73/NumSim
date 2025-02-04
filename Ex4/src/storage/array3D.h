#pragma once

#include <cassert>
#include <array>
#include <vector>
#include <sstream>
#include <iostream>
#include <exception>
#include <limits>
#include <mpi.h>

/** This class represents a 2D array of double values.
 *  Internally they are stored consecutively in memory.
 *  The entries can be accessed by two indices i,j.
 */
class Array3D {
public:
  //! constructor
  Array3D(std::array<int, 3> size);

  //! optional name argument for debugging convenience
  Array3D(std::array<int, 3> size, std::string name);

  //! get the size
  inline std::array<int, 3> size() const { return size_; }

  //! get the length of the underlying data
  inline std::size_t length() const { return size_[0]*size_[1]*size_[2]; }

  //! get index (to avoid it beeing different between both access functions)
  inline std::size_t compute_index(std::size_t i, std::size_t j, std::size_t k)
  const
  {
    return k * size0Xsize1_ + j * size_[0] + i;
  }

  // often used items are inlined for better optimisation
  //! access the value at coordinate (i,j), declared not const, i.e. the value
  //! can be changed
  inline double &operator()(int i, int j, int k)
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

  //! get the value at coordinate (i,j), declared const, i.e. it is not possible
  //! to change the value
  inline double operator()(int i, int j, int k) const
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

protected:
  std::vector<double> data_;      //< storage array values, in row-major order
  const std::array<int, 3> size_; //< width, height of the domain
  std::string name_;              //< name used for debugging
  const std::size_t size0Xsize1_; //< precomputed size of the cube area for efficiency
};