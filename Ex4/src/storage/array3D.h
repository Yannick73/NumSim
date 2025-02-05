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

  //! optional name argument for debugging convenience
  Array3D(std::array<int, 3> size, std::string name = "unnamed");

  //! get the size
  inline std::array<int, 3> size() const { return size_; }

  //! get the length of the underlying data
  inline std::size_t length() const { return data_.size(); }

  //! gets the underlying buffer definition
  inline double* data() { return data_.data(); };

  //! get index (to avoid it beeing different between both access functions)
  inline std::size_t compute_index(std::size_t i, std::size_t j, std::size_t k)
  const
  {
    return k * size0Xsize1_ + j * size_[0] + i;
  }

  // often used items are inlined for better optimisation
  //! access the value at coordinate (i,j), declared not const, i.e. the value
  //! can be changed
  virtual double &operator()(int i, int j, int k);

  //! get the value at coordinate (i,j), declared const, i.e. it is not possible
  //! to change the value
  virtual double operator()(int i, int j, int k) const;

  inline void rename(std::string name) { name_ = name; }

protected:
  std::vector<double> data_;      //< storage array values, in row-major order
  const std::array<int, 3> size_; //< width, height of the domain
  std::string name_;              //< name used for debugging
  const std::size_t size0Xsize1_; //< precomputed size of the cube area for efficiency
};