#include "discretization/central_differences.h"

//! compute the 1st derivative ∂ u^2 / ∂x
double CentralDifferences::computeDu2Dx(int i, int j, int k) const
{
  const double uRight = (u(i,j,k)   + u(i+1,j,k)) / 2.;
  const double uLeft  = (u(i-1,j,k) + u(i,j,k))   / 2.;
  return (uRight*uRight - uLeft*uLeft) / dx();
}

//! compute the 1st derivative ∂ v^2 / ∂y
double CentralDifferences::computeDv2Dy(int i, int j, int k) const 
{
  const double vTop    = (v(i,j,k)   + v(i,j+1,k)) / 2.;
  const double vBottom = (v(i,j-1,k) + v(i,j,k))   / 2.;
  return (vTop*vTop - vBottom*vBottom) / dy();
}

//! compute the 1st derivative ∂ w^2 / ∂z
double CentralDifferences::computeDw2Dz(int i, int j, int k) const 
{
  const double wFront = (w(i,j,k)   + w(i,j,k+1)) / 2.;
  const double vBack  = (w(i,j,k-1) + w(i,j,k))   / 2.;
  return (wFront*wFront - vBack*vBack) / dz();
}

// Assumption: In the uv plane, the discretization should remain unchanged
//! compute the 1st derivative ∂ (uv) / ∂x
double CentralDifferences::computeDuvDx(int i, int j, int k) const 
{
  const double vRight = (v(i+1,j,k)   + v(i,j,k))   / 2.; // v at top right corner of cell
  const double vLeft  = (v(i,j,k)     + v(i-1,j,k)) / 2.; // v at top left corner of cell
  const double uRight = (u(i,j+1,k)   + u(i,j,k))   / 2.; // u at top right corner of cell
  const double uLeft  = (u(i-1,j+1,k) + u(i-1,j,k)) / 2.; // u at top left corner of cell
  return (vRight*uRight - vLeft*uLeft) / dx();
}

//! compute the 1st derivative ∂ (uv) / ∂y
double CentralDifferences::computeDuvDy(int i, int j, int k) const 
{
  const double vTop    = (v(i,j,k)   + v(i+1,j,k))   / 2.; // v at top right corner of cell
  const double vBottom = (v(i,j-1,k) + v(i+1,j-1,k)) / 2.; // v at bottom right corner of cell
  const double uTop    = (u(i,j,k)   + u(i,j+1,k))   / 2.; // u at top right corner of cell
  const double uBottom = (u(i,j-1,k) + u(i,j,k))     / 2.; // u at bottom right corner of cell
  return (vTop*uTop - vBottom*uBottom) / dy();
}


//! compute the 1st derivative ∂ (uw) / ∂x
double CentralDifferences::computeDuwDx(int i, int j, int k) const 
{
  const double uLeft  = (u(i-1,j,k)   + u(i,  j,k)) / 2.;
  const double wLeft  = (w(i-1,j,k+1) + w(i-1,j,k)) / 2.;
  const double uRight = (u(i+1,j,k)   + u(i,  j,k)) / 2.;
  const double wRight = (w(i,  j,k+1) + w(i,  j,k)) / 2.; 
  return (uLeft*wLeft - uRight*wRight) / dx();
}

//! compute the 1st derivative ∂ (uw) / ∂z
double CentralDifferences::computeDuwDz(int i, int j, int k) const 
{
  // back to front z increases
  const double wBack  = (w(i,  j,k-1) + w(i,j,k))   / 2.;
  const double uBack  = (u(i+1,j,k-1) + u(i,j,k-1)) / 2.;
  const double wFront = (w(i,  j,k+1) + w(i,j,k))   / 2.;
  const double uFront = (u(i+1,j,k)   + u(i,j,k))   / 2.;
  return (wBack*uBack - wFront*uFront) / dz();
}

//! compute the 1st derivative ∂ (vw) / ∂y
double CentralDifferences::computeDvwDy(int i, int j, int k) const 
{
  const double vTop    = (v(i,j+1,k)   + v(i,j,  k)) / 2.;
  const double wTop    = (w(i,j,  k+1) + w(i,j,  k)) / 2.;
  const double vBottom = (v(i,j,  k-1) + w(i,j,  k)) / 2.;
  const double wBottom = (w(i,j-1,k+1) + w(i,j-1,k)) / 2.;
  return (vBottom*wBottom - vTop*wTop) / dy();
}

//! compute the 1st derivative ∂ (vw) / ∂z
double CentralDifferences::computeDvwDz(int i, int j, int k) const
{
  const double wBack  = (w(i,j,  k-1) + w(i,j,k))   / 2.;
  const double vBack  = (v(i,j+1,k-1) + v(i,j,k-1)) / 2.;
  const double wFront = (w(i,j,  k+1) + w(i,j,k))   / 2.;
  const double vFront = (v(i,j+1,k)   + v(i,j,k))   / 2.;
  return (vBack*wBack - vFront*wFront) / dz();
}