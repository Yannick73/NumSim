#pragma once

#include <array>
#include "storage/field_variable.h"
#include "settings.h"

class StaggeredGrid {
public:
  StaggeredGrid(std::array<int, 2> nCells, std::array<double, 2> meshWidth);
  inline const std::array<double, 2> meshWidth() const { return meshWidth_; }
  inline const std::array<int, 2> nCells() const { return nCells_; }
  inline const long nTotalCells() const { return nCells_[0]*nCells_[1]; }
  // to access the variables
  inline FieldVariable &u() { return u_; }
  inline FieldVariable &v() { return v_; }
  inline FieldVariable &p() { return p_; }
  inline FieldVariable &f() { return f_; }
  inline FieldVariable &g() { return g_; }
  inline FieldVariable &rhs() { return rhs_; }

  // end variables for u and v respectively
  inline int uiEnd() const { return nCells_[0]-1; }
  inline int ujEnd() const { return nCells_[1]; }
  inline int viEnd() const { return nCells_[0]; }
  inline int vjEnd() const { return nCells_[1]-1; }
  inline int iEnd()  const { return nCells_[0]; }
  inline int jEnd()  const { return nCells_[1]; }

  // inlining to make access more efficient
  inline double  u(int i, int j) const { return u_(i+2,j+2); }
  inline double &u(int i, int j)       { return u_(i+2,j+2); }
  inline double  v(int i, int j) const { return v_(i+2,j+2); }
  inline double &v(int i, int j)       { return v_(i+2,j+2); }
  inline double  p(int i, int j) const { return p_(i+2,j+2); }
  inline double &p(int i, int j)       { return p_(i+2,j+2); }
  inline double  f(int i, int j) const { return f_(i+2,j+2); }
  inline double &f(int i, int j)       { return f_(i+2,j+2); }
  inline double  g(int i, int j) const { return g_(i+2,j+2); }
  inline double &g(int i, int j)       { return g_(i+2,j+2); }
  inline double  rhs(int i, int j) const { return rhs_(i+2,j+2); }
  inline double &rhs(int i, int j)       { return rhs_(i+2,j+2); }
  void setBoundary(Settings settings);
  void setBoundaryValues();

protected:
  std::array<double, 2> meshWidth_;
  std::array<int, 2> nCells_;
  FieldVariable u_, v_, p_, rhs_, f_, g_;
  std::array<double, 2> dirichletBcBottom;
  std::array<double, 2> dirichletBcTop;
  std::array<double, 2> dirichletBcLeft;
  std::array<double, 2> dirichletBcRight;
};
