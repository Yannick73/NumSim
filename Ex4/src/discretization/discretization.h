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
  Discretization(PartitionInformation &pi, const Settings &settings);

  //! calculate the deltaT depending on the velocity and mesh-width, includes boundaries
  double calculateVelocityDelta() const;

  //! calculate the deltaT depending on the reynolds constant and mesh-width
  double calculateReynoldsDelta() const;

  //! calculate the preliminary velocities (f, g) and the rhs for pressure solver
  void calculateFGH(double deltaT);

  //! calculate the preliminary velocities (f, g) and the rhs for pressure solver
  void calculateRHS(double deltaT);

  //! calculate the final velocities using f, g and p
  void calculateUVW(double deltaT);

  //! return the meshwidth in x-direction
  inline double dx() const { return meshWidth_[0]; };

  //! return the meshwidth in y-direction
  inline double dy() const { return meshWidth_[1]; };

  //! return the meshwidth in z-direction
  inline double dz() const { return meshWidth_[2]; };

  //! return the meshwidth squared in x-direction
  inline double dx2() const { return meshWidth2_[0]; };

  //! return the meshwidth squared in y-direction
  inline double dy2() const { return meshWidth2_[1]; };

  //! return the meshwidth squared in z-direction
  inline double dz2() const { return meshWidth2_[2]; };

  //! compute the 1st derivative ∂p / ∂x
  inline double computeDpDx(int i, int j, int k) const { return (p(i+1,j,k) - p(i,j,k)) / dx(); };

  //! compute the 1st derivative ∂p / ∂y
  inline double computeDpDy(int i, int j, int k) const { return (p(i,j+1,k) - p(i,j,k)) / dy(); };

  //! compute the 1st derivative ∂p / ∂z
  inline double computeDpDz(int i, int j, int k) const { return (p(i,j,k+1) - p(i,j,k)) / dz(); };

  //! compute the 2st derivative ∂^2p / ∂x^2
  inline double computeD2pDx2(int i, int j, int k) const { return (p(i+1,j,k) - 2*p(i,j,k) + p(i-1,j,k)) / dx2(); };

  //! compute the 2st derivative ∂^2p / ∂y^2
  inline double computeD2pDy2(int i, int j, int k) const { return (p(i,j+1,k) - 2*p(i,j,k) + p(i,j-1,k)) / dy2(); };

  //! compute the 2st derivative ∂^2p / ∂z^2
  inline double computeD2pDz2(int i, int j, int k) const { return (p(i,j,k+1) - 2*p(i,j,k) + p(i,j,k-1)) / dz2(); };

  //! compute the 2nd derivative ∂^2 u / ∂x^2
  inline double computeD2uDx2(int i, int j, int k) const { return (u(i+1,j,k) - 2*u(i,j,k) + u(i-1,j,k)) / dx2(); };

  //! compute the 2nd derivative ∂^2 u / ∂y^2
  inline double computeD2uDy2(int i, int j, int k) const { return (u(i,j+1,k) - 2*u(i,j,k) + u(i,j-1,k)) / dy2(); };

  //! compute the 2nd derivative ∂^2 u / ∂z^2
  inline double computeD2uDz2(int i, int j, int k) const { return (u(i,j,k+1) - 2*u(i,j,k) + u(i,j,k-1)) / dz2(); };

  //! compute the 2nd derivative ∂^2 v / ∂x^2
  inline double computeD2vDx2(int i, int j, int k) const { return (v(i+1,j,k) - 2*v(i,j,k) + v(i-1,j,k)) / dx2(); };

  //! compute the 2nd derivative ∂^2 v / ∂y^2
  inline double computeD2vDy2(int i, int j, int k) const { return (v(i,j+1,k) - 2*v(i,j,k) + v(i,j-1,k)) / dy2(); };

  //! compute the 2nd derivative ∂^2 v / ∂z^2
  inline double computeD2vDz2(int i, int j, int k) const { return (v(i,j,k+1) - 2*v(i,j,k) + v(i,j,k-1)) / dz2(); };

  //! compute the 2nd derivative ∂^2 w / ∂x^2
  inline double computeD2wDx2(int i, int j, int k) const { return (w(i+1,j,k) - 2*w(i,j,k) + w(i-1,j,k)) / dx2(); };

  //! compute the 2nd derivative ∂^2 w / ∂y^2
  inline double computeD2wDy2(int i, int j, int k) const { return (w(i,j+1,k) - 2*w(i,j,k) + w(i,j-1,k)) / dy2(); };

  //! compute the 2nd derivative ∂^2 w / ∂z^2
  inline double computeD2wDz2(int i, int j, int k) const { return (w(i,j,k+1) - 2*w(i,j,k) + w(i,j,k-1)) / dz2(); };

  //! compute the 1st derivative ∂ u^2 / ∂x
  virtual double computeDu2Dx(int i, int j, int k) const = 0;

  //! compute the 1st derivative ∂ v^2 / ∂y
  virtual double computeDv2Dy(int i, int j, int k) const = 0;

  //! compute the 1st derivative ∂ w^2 / ∂z
  virtual double computeDw2Dz(int i, int j, int k) const = 0;

  //! compute the 1st derivative ∂ (uv) / ∂x
  virtual double computeDuvDx(int i, int j, int k) const = 0;

  //! compute the 1st derivative ∂ (uv) / ∂y
  virtual double computeDuvDy(int i, int j, int k) const = 0;

  //! compute the 1st derivative ∂ (uw) / ∂x
  virtual double computeDuwDx(int i, int j, int k) const = 0;

  //! compute the 1st derivative ∂ (uw) / ∂z
  virtual double computeDuwDz(int i, int j, int k) const = 0;

  //! compute the 1st derivative ∂ (vw) / ∂y
  virtual double computeDvwDy(int i, int j, int k) const = 0;

  //! compute the 1st derivative ∂ (vw) / ∂z
  virtual double computeDvwDz(int i, int j, int k) const = 0;

private:
  //! Settings used
  const Settings &settings_;
  //! mesh width squared used for second derivates, or more precise dx2()
  const std::array<double, 3> meshWidth2_;
};
