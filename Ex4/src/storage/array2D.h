#pragma once

#include <cassert>
#include <array>
#include <vector>
#include <sstream>
#include <iostream>
#include <exception>
#include <mpi.h>

/** This class represents a 2D array of double values.
 *  Internally they are stored consecutively in memory.
 *  The entries can be accessed by two indices i,j.
 */
class Array2D {
public:
  //! constructor
  Array2D(std::array<int, 2> size);

  //! optional name argument for debugging convenience
  Array2D(std::array<int, 2> size, std::string name);

  //! get the size
  inline std::array<int, 2> size() const { return size_; }

  //! gets the underlying buffer definition
  inline double* data() { return data_.data(); };

  //! gets the total length
  inline std::size_t length() const { return data_.size(); }

  // often used items are inlined for better optimisation
  //! access the value at coordinate (i,j), declared not const, i.e. the value
  //! can be changed
  inline double &operator()(int i, int j)
  {
    const int index = j * size_[0] + i;
    // make assertion conditional on DEBUG mode to optimize further
    #ifndef NDEBUG
    if(i < 0 || i >= size_[0] || j < 0 || j >= size_[1])
    {
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      std::stringstream str;
      str << "Out-of-bound access on " << name_ << "(i,j): (" << i << ',' << j
          << "), size: (" << size_[0] << ',' << size_[1] << ") in R:" << rank << "\n";
      throw std::out_of_range(str.str());
    }
    #endif

    return data_[index];
  }

  //! get the value at coordinate (i,j), declared const, i.e. it is not possible
  //! to change the value
  inline double operator()(int i, int j) const
  {
    const int index = j * size_[0] + i;
    // assert that indices are in range
    #ifndef NDEBUG
    if(i < 0 || i >= size_[0] || j < 0 || j >= size_[1])
    {
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      std::stringstream str;
      str << "Out-of-bound access on " << name_ << "(i,j): (" << i << ',' << j
          << "), size: (" << size_[0] << ',' << size_[1] << ") in R:" << rank << "\n";
      throw std::out_of_range(str.str());
    }
    #endif

    return data_[index];
  }

protected:
  std::vector<double> data_;      //< storage array values, in row-major order
  const std::array<int, 2> size_; //< width, height of the domain
  std::string name_;
};