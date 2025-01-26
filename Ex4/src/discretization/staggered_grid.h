#pragma once

#include <array>
#include <memory>
#include "discretization/partition_information.h"
#include "storage/field_variable.h"
#include "settings.h"

class StaggeredGrid {
public:
  //! normal constructor with number of cells and width of each cell
  // StaggeredGrid(std::array<int, 3> nCells, std::array<double, 3> meshWidth);

  //! constructor using the partition information (which includes info about boundaries)
  StaggeredGrid(PartitionInformation &pi);

  void makeCGFields();


  // cell width and number
  inline const std::array<double, 3> meshWidth() const { return meshWidth_; }
  inline const std::array<int, 3> nCells() const { return nCells_; }
  inline const std::size_t nTotalCells() const { return nCells_[0]*nCells_[1]*nCells_[2]; }
  
  // to access the variables
  inline FieldVariable &u() { return u_; }
  inline FieldVariable &v() { return v_; }
  inline FieldVariable &w() { return w_; }
  inline FieldVariable &p() { return p_; }
  inline FieldVariable &f() { return f_; }
  inline FieldVariable &g() { return g_; }
  inline FieldVariable &h() { return h_; }
  inline FieldVariable &rhs() { return rhs_; }
  // CG variables are only initialized when required
  std::shared_ptr<FieldVariable> r() { return r_; }
  std::shared_ptr<FieldVariable> a() { return a_; }

  // begin of each field, necassary for possible ghost layer
  inline int ui0() const { return ui0_; }
  inline int uj0() const { return uj0_; }
  inline int uk0() const { return uk0_; }
  inline int vi0() const { return vi0_; }
  inline int vj0() const { return vj0_; }
  inline int vk0() const { return vk0_; }
  inline int wi0() const { return wi0_; }
  inline int wj0() const { return wj0_; }
  inline int wk0() const { return wk0_; }
  inline int pi0() const { return pi0_; }
  inline int pj0() const { return pj0_; }
  inline int pk0() const { return pk0_; }
  // end variables for u and v respectively
  inline int uiN() const { return nCells_[0] - 1 + uGhost_; }
  inline int ujN() const { return nCells_[1]; }
  inline int ukN() const { return nCells_[2]; }
  inline int viN() const { return nCells_[0]; }
  inline int vjN() const { return nCells_[1] - 1 + vGhost_; }
  inline int vkN() const { return nCells_[2]; }
  inline int wiN() const { return nCells_[0]; }
  inline int wjN() const { return nCells_[1]; }
  inline int wkN() const { return nCells_[2] - 1 + wGhost_; }
  inline int piN() const { return nCells_[0]; }
  inline int pjN() const { return nCells_[1]; }
  inline int pkN() const { return nCells_[2]; }

  // inlining to make access more efficient
  inline double  u(int i, int j, int k) const { return u_(i+ui0_,j+uj0_,k+uk0_); }
  inline double &u(int i, int j, int k)       { return u_(i+ui0_,j+uj0_,k+uk0_); }
  inline double  v(int i, int j, int k) const { return v_(i+vi0_,j+vj0_,k+vk0_); }
  inline double &v(int i, int j, int k)       { return v_(i+vi0_,j+vj0_,k+vk0_); }
  inline double  w(int i, int j, int k) const { return w_(i+wi0_,j+wj0_,k+wk0_); }
  inline double &w(int i, int j, int k)       { return w_(i+wi0_,j+wj0_,k+wk0_); }
  inline double  p(int i, int j, int k) const { return p_(i+pi0_,j+pj0_,k+pk0_); }
  inline double &p(int i, int j, int k)       { return p_(i+pi0_,j+pj0_,k+pk0_); }
  inline double  f(int i, int j, int k) const { return f_(i+ui0_,j+uj0_,k+uk0_); }
  inline double &f(int i, int j, int k)       { return f_(i+ui0_,j+uj0_,k+uk0_); }
  inline double  g(int i, int j, int k) const { return g_(i+vi0_,j+vj0_,k+vk0_); }
  inline double &g(int i, int j, int k)       { return g_(i+vi0_,j+vj0_,k+vk0_); }
  inline double  h(int i, int j, int k) const { return h_(i+wi0_,j+wj0_,k+wk0_); }
  inline double &h(int i, int j, int k)       { return h_(i+wi0_,j+wj0_,k+wk0_); }
  inline double  rhs(int i, int j, int k) const { return rhs_(i,j,k); }
  inline double &rhs(int i, int j, int k)       { return rhs_(i,j,k); }
  // CG solver variables use same indexing schema as p
  inline double  r(int i, int j, int k) const { return (*r_)(i+pi0_,j+pj0_,k+pk0_); }
  inline double &r(int i, int j, int k)       { return (*r_)(i+pi0_,j+pj0_,k+pk0_); }
  inline double  a(int i, int j, int k) const { return (*a_)(i+pi0_,j+pj0_,k+pk0_); }
  inline double &a(int i, int j, int k)       { return (*a_)(i+pi0_,j+pj0_,k+pk0_); }

protected:
  //! normal fluid-sim variables
  std::array<double, 3> meshWidth_;
  std::array<int, 3> nCells_;
  FieldVariable u_, v_, w_, p_, rhs_, f_, g_, h_;
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
  const int uk0_;
  const int vi0_;
  const int vj0_;
  const int vk0_;
  const int wi0_;
  const int wj0_;
  const int wk0_;
  const int pi0_;
  const int pj0_;
  const int pk0_;

  const int uGhost_;
  const int vGhost_;
  const int wGhost_;
  
  bool cgFieldsMade_ = false;
};
