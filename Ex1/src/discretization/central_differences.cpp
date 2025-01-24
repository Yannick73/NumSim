#include "discretization/central_differences.h"

CentralDifferences::CentralDifferences(std::array<int, 2> nCells,
                                       std::array<double, 2> meshWidth,
                                       Settings settings)
    : Discretization(nCells, meshWidth, settings) {}

//! compute the 1st derivative ∂ u^2 / ∂x
double CentralDifferences::computeDu2Dx(int i, int j) const
{
  const double uRight = (u(i,j)   + u(i+1,j)) / 2.;
  const double uLeft  = (u(i-1,j) + u(i,j))   / 2.;
  return (uRight*uRight - uLeft*uLeft) / dx();
}

//! compute the 1st derivative ∂ v^2 / ∂x
double CentralDifferences::computeDv2Dy(int i, int j) const 
{
  const double vTop    = (v(i,j)   + v(i,j+1)) / 2.;
  const double vBottom = (v(i,j-1) + v(i,j))   / 2.;
  return (vTop*vTop - vBottom*vBottom) / dy();
}

//! compute the 1st derivative ∂ (uv) / ∂x
double CentralDifferences::computeDuvDx(int i, int j) const 
{
  const double vRight = (v(i+1,j)   + v(i,j))   / 2.; // v at top right corner of cell
  const double vLeft  = (v(i,j)     + v(i-1,j)) / 2.; // v at top left corner of cell
  const double uRight = (u(i,j+1)   + u(i,j))   / 2.; // u at top right corner of cell
  const double uLeft  = (u(i-1,j+1) + u(i-1,j)) / 2.; // u at top left corner of cell
  return (vRight*uRight - vLeft*uLeft) / dx();
}

//! compute the 1st derivative ∂ (uv) / ∂y
double CentralDifferences::computeDuvDy(int i, int j) const 
{
  const double vTop    = (v(i,j)   + v(i+1,j))   / 2.; // v at top right corner of cell
  const double vBottom = (v(i,j-1) + v(i+1,j-1)) / 2.; // v at bottom right corner of cell
  const double uTop    = (u(i,j)   + u(i,j+1))   / 2.; // u at top right corner of cell
  const double uBottom = (u(i,j-1) + u(i,j))     / 2.; // u at bottom right corner of cell
  return (vTop*uTop - vBottom*uBottom) / dy();
}