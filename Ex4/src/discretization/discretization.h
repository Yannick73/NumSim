#pragma once

#include <array>
#include <cassert>
#include "discretization/partition_information.h"
#include "discretization/staggered_grid.h"
#include "settings.h"

class Discretization : public StaggeredGrid
{
public:
  //! construct the object with given number of cells in x and y direction
  Discretization(std::array<int, 2> nCells, std::array<double, 2> meshWidth, const Settings &settings);

  //! construct the object with given number of cells in x and y direction
  Discretization(PartitionInformation pi, const Settings &settings);

  //! calculate the deltaT depending on the velocity and mesh-width, includes boundaries
  double calculateVelocityDelta() const;

  //! calculate the deltaT depending on the reynolds constant and mesh-width
  double calculateReynoldsDelta() const;

  //! calculate the preliminary velocities (f, g) and the rhs for pressure solver
  void calculateFG(double deltaT);

  //! calculate the preliminary velocities (f, g) and the rhs for pressure solver
  void calculateRHS(double deltaT);

  //! calculate the final velocities using f, g and p
  void calculateUV(double deltaT);

  //! return the meshwidth in x-direction
  inline double dx() const { return meshWidth_[0]; };

  //! return the meshwidth in y-direction
  inline double dy() const { return meshWidth_[1]; };

  //! return the meshwidth squared in x-direction
  inline double dx2() const { return meshWidth2_[0]; };

  //! return the meshwidth squared in y-direction
  inline double dy2() const { return meshWidth2_[1]; };

  //! compute the 1st derivative ∂p / ∂x
  inline double computeDpDx(int i, int j) const { return (p(i+1,j) - p(i,j)) / dx(); };

  //! compute the 1st derivative ∂p / ∂y
  inline double computeDpDy(int i, int j) const { return (p(i,j+1) - p(i,j)) / dy(); };

  //! compute the 2st derivative ∂^2p / ∂x^2
  inline double computeD2pDx2(int i, int j) const { return (p(i-1,j) - 2*p(i,j) + p(i+1,j)) / dx2(); };

  //! compute the 2st derivative ∂^2p / ∂y^2
  inline double computeD2pDy2(int i, int j) const { return (p(i,j-1) - 2*p(i,j) + p(i,j+1)) / dy2(); };

  //! compute the 2st derivative ∂^2a / ∂x^2
  inline double computeD2aDx2(int i, int j) const { return (a(i-1,j) - 2*a(i,j) + a(i+1,j)) / dx2(); };

  //! compute the 2st derivative ∂^2a / ∂y^2
  inline double computeD2aDy2(int i, int j) const { return (a(i,j-1) - 2*a(i,j) + a(i,j+1)) / dy2(); };

  //! compute the 2nd derivative ∂^2 u / ∂x^2
  inline double computeD2uDx2(int i, int j) const { return (u(i+1,j) - 2*u(i,j) + u(i-1,j)) / dx2(); };

  //! compute the 2nd derivative ∂^2 u / ∂y^2
  inline double computeD2uDy2(int i, int j) const { return (u(i,j+1) - 2*u(i,j) + u(i,j-1)) / dy2(); };

  //! compute the 2nd derivative ∂^2 v / ∂x^2
  inline double computeD2vDx2(int i, int j) const { return (v(i+1,j) - 2*v(i,j) + v(i-1,j)) / dx2(); };

  //! compute the 2nd derivative ∂^2 v / ∂y^2
  inline double computeD2vDy2(int i, int j) const { return (v(i,j+1) - 2*v(i,j) + v(i,j-1)) / dy2(); };

  //! compute the 1st derivative ∂ u^2 / ∂x
  virtual double computeDu2Dx(int i, int j) const = 0;

  //! compute the 1st derivative ∂ v^2 / ∂y
  virtual double computeDv2Dy(int i, int j) const = 0;

  //! compute the 1st derivative ∂ (uv) / ∂x
  virtual double computeDuvDx(int i, int j) const = 0;

  //! compute the 1st derivative ∂ (uv) / ∂y
  virtual double computeDuvDy(int i, int j) const = 0;

private:
  //! Settings used
  const Settings &settings_;
  //! mesh width squared used for second derivates, or more precise dx2()
  const std::array<double, 2> meshWidth2_;
};
