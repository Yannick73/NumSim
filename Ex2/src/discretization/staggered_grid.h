#pragma once

#include <array>
#include <memory>
#include "discretization/partition_information.h"
#include "storage/field_variable.h"
#include "settings.h"

class StaggeredGrid {
public:
  //! normal constructor with number of cells and width of each cell
  StaggeredGrid(std::array<int, 2> nCells, std::array<double, 2> meshWidth);

  //! constructor using the partition information (which includes info about boundaries)
  StaggeredGrid(PartitionInformation &pi);

  void makeCGFields();


  // cell width and number
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
  std::shared_ptr<FieldVariable> r() { return r_; }
  std::shared_ptr<FieldVariable> a() { return a_; }

  // begin of each field, necassary for possible ghost layer
  inline int ui0() const { return ui0_; }
  inline int uj0() const { return uj0_; }
  inline int vi0() const { return vi0_; }
  inline int vj0() const { return vj0_; }
  inline int pi0() const { return pi0_; }
  inline int pj0() const { return pj0_; }
  // end variables for u and v respectively
  inline int uiN() const { return nCells_[0] - 1 + uGhost_; }
  inline int ujN() const { return nCells_[1]; }
  inline int viN() const { return nCells_[0]; }
  inline int vjN() const { return nCells_[1] - 1 + vGhost_; }
  inline int piN() const { return nCells_[0]; }
  inline int pjN() const { return nCells_[1]; }

  // inlining to make access more efficient
  inline double  u(int i, int j) const { return u_(i+ui0_,j+uj0_); }
  inline double &u(int i, int j)       { return u_(i+ui0_,j+uj0_); }
  inline double  v(int i, int j) const { return v_(i+vi0_,j+vj0_); }
  inline double &v(int i, int j)       { return v_(i+vi0_,j+vj0_); }
  inline double  p(int i, int j) const { return p_(i+pi0_,j+pj0_); }
  inline double &p(int i, int j)       { return p_(i+pi0_,j+pj0_); }
  inline double  f(int i, int j) const { return f_(i+ui0_,j+uj0_); }
  inline double &f(int i, int j)       { return f_(i+ui0_,j+uj0_); }
  inline double  g(int i, int j) const { return g_(i+vi0_,j+vj0_); }
  inline double &g(int i, int j)       { return g_(i+vi0_,j+vj0_); }
  inline double  rhs(int i, int j) const { return rhs_(i,j); }
  inline double &rhs(int i, int j)       { return rhs_(i,j); }
  // CG solver variables use same indexing schema as p
  inline double  r(int i, int j) const { return (*r_)(i+pi0_,j+pj0_); }
  inline double &r(int i, int j)       { return (*r_)(i+pi0_,j+pj0_); }
  inline double  a(int i, int j) const { return (*a_)(i+pi0_,j+pj0_); }
  inline double &a(int i, int j)       { return (*a_)(i+pi0_,j+pj0_); }

protected:
  //! normal fluid-sim variables
  std::array<double, 2> meshWidth_;
  std::array<int, 2> nCells_;
  FieldVariable u_, v_, p_, rhs_, f_, g_;
  //! CG solver variables, only defined when using said CG, same indexing as p
  std::shared_ptr<FieldVariable> r_, a_;
  /*std::array<double, 2> dirichletBcBottom;
  std::array<double, 2> dirichletBcTop;
  std::array<double, 2> dirichletBcLeft;
  std::array<double, 2> dirichletBcRight;*/

private:
  // field begin offsets
  const int ui0_;
  const int uj0_;
  const int vi0_;
  const int vj0_;
  const int pi0_;
  const int pj0_;

  const int uGhost_;
  const int vGhost_;
  
  bool cgFieldsMade_ = false;
};
