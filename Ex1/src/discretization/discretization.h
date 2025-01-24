#pragma once

#include <array>
#include <cassert>
#include "discretization/staggered_grid.h"
#include "settings.h"

class Discretization : public StaggeredGrid {

public:
  //! construct the object with given number of cells in x and y direction
  Discretization(std::array<int, 2> nCells, std::array<double, 2> meshWidth, Settings settings);

  //! calculate the deltaT depending on the velocity and mesh-width, includes boundaries
  double calculateVelocityDelta() const;

  //! calculate the preliminary velocities (f, g) and the rhs for pressure solver
  void calculateFG_RHS(double deltaT);

  //! calculate the final velocities using f, g and p
  void calculateUV(double deltaT);

  //! return the meshwidth in x-direction
  inline double dx() const { return meshWidth_[0]; };

  //! return the meshwidth in y-direction
  inline double dy() const { return meshWidth_[1]; };

  //! return the meshwidth squared in x-direction
  inline double dx2() const { return dx()*dx(); };

  //! return the meshwidth squared in y-direction
  inline double dy2() const { return dy()*dy(); };

  //! compute the 1st derivative ∂p / ∂x
  double computeDpDx(int i, int j) const;

  //! compute the 1st derivative ∂p / ∂y
  double computeDpDy(int i, int j) const;

  //! compute the 2st derivative ∂^2p / ∂x^2
  double computeD2pDx2(int i, int j) const;

  //! compute the 2st derivative ∂^2p / ∂y^2
  double computeD2pDy2(int i, int j) const;

  //! compute the 2nd derivative ∂^2 u / ∂x^2
  double computeD2uDx2(int i, int j) const;

  //! compute the 2nd derivative ∂^2 u / ∂y^2
  double computeD2uDy2(int i, int j) const;

  //! compute the 2nd derivative ∂^2 v / ∂x^2
  double computeD2vDx2(int i, int j) const;

  //! compute the 2nd derivative ∂^2 v / ∂y^2
  double computeD2vDy2(int i, int j) const;

  //! compute the 1st derivative ∂ u^2 / ∂x
  virtual double computeDu2Dx(int i, int j) const = 0;

  //! compute the 1st derivative ∂ v^2 / ∂y
  virtual double computeDv2Dy(int i, int j) const = 0;

  //! compute the 1st derivative ∂ (uv) / ∂x
  virtual double computeDuvDx(int i, int j) const = 0;

  //! compute the 1st derivative ∂ (uv) / ∂y
  virtual double computeDuvDy(int i, int j) const = 0;

private:
  const Settings settings_;
};
