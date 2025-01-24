#pragma once

#include <cassert>
#include <array>
#include <vector>

/** This class represents a 2D array of double values.
 *  Internally they are stored consecutively in memory.
 *  The entries can be accessed by two indices i,j.
 */
class Array2D {
public:
  //! constructor
  Array2D(std::array<int, 2> size);

  //! get the size
  inline std::array<int, 2> size() const { return size_; }

  // often used items are inlined for better optimisation
  //! access the value at coordinate (i,j), declared not const, i.e. the value
  //! can be changed
  inline double &operator()(int i, int j)
  {
    const int index = j * size_[0] + i;
    // assert that indices are in range
    assert(0 <= i && i < size_[0]);
    assert(0 <= j && j < size_[1]);

    return data_[index];
  }

  //! get the value at coordinate (i,j), declared const, i.e. it is not possible
  //! to change the value
  inline double operator()(int i, int j) const
  {

    const int index = j * size_[0] + i;
    // assert that indices are in range
    assert(0 <= i && i < size_[0]);
    assert(0 <= j && j < size_[1]);

    return data_[index];
  }

protected:
  std::vector<double> data_;      //< storage array values, in row-major order
  const std::array<int, 2> size_; //< width, height of the domain
};