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
  const double wBack  = (w(i,j,k-1) + w(i,j,k))   / 2.;
  return (wFront*wFront - wBack*wBack) / dz();
}

// Assumption: In the uv plane, the discretization should remain unchanged
//! compute the 1st derivative ∂ (uv) / ∂x on the v points
double CentralDifferences::computeDuvDx(int i, int j, int k) const 
{
  const double vRight = (v(i+1,j,  k) + v(i,  j,k)) / 2.; // v i+1/2,j,k
  const double uRight = (u(i,  j+1,k) + u(i,  j,k)) / 2.; // u i,j+1/2,k
  const double vLeft  = (v(i,  j,  k) + v(i-1,j,k)) / 2.; // v i-1/2,j,k
  const double uLeft  = (u(i-1,j+1,k) + u(i-1,j,k)) / 2.; // u i-1,j+1/2,k
  return (vRight*uRight - vLeft*uLeft) / dx();
}

//! compute the 1st derivative ∂ (uv) / ∂y on the u points
double CentralDifferences::computeDuvDy(int i, int j, int k) const 
{
  const double vTop    = (v(i,j,  k) + v(i+1,j,  k)) / 2.; // v i+1/2,j,k
  const double uTop    = (u(i,j,  k) + u(i,  j+1,k)) / 2.; // u i,j+1/2,k
  const double vBottom = (v(i,j-1,k) + v(i+1,j-1,k)) / 2.; // v v+1/2,j-1,k
  const double uBottom = (u(i,j-1,k) + u(i,  j,  k)) / 2.; // u i,j-1/2,k
  return (vTop*uTop - vBottom*uBottom) / dy();
}


//! compute the 1st derivative ∂ (uw) / ∂x on the w points
double CentralDifferences::computeDuwDx(int i, int j, int k) const 
{
  const double uLeft  = (u(i-1,j,k+1) + u(i-1,j,k)) / 2.; // u i-1,j,k+1/2
  const double wLeft  = (w(i-1,j,k)   + w(i,  j,k)) / 2.; // w i-1/2,j,k
  const double uRight = (u(i,  j,k+1) + u(i,  j,k)) / 2.; // u i,j,k+1/2
  const double wRight = (w(i+1,j,k)   + w(i,  j,k)) / 2.; // w i+1/2,k (!+1/2)
  return (uRight*wRight - uLeft*wLeft) / dx();  // wrong order!
}

//! compute the 1st derivative ∂ (uw) / ∂z on the u points
double CentralDifferences::computeDuwDz(int i, int j, int k) const 
{
  // back to front z increases
  const double wBack  = (w(i+1,j,k-1) + w(i,j,k-1)) / 2.; // w i+1/2,j,k-1
  const double uBack  = (u(i,  j,k-1) + u(i,j,k))   / 2.; // u i,j,k-1/2
  const double wFront = (w(i+1,j,k)   + w(i,j,k))   / 2.; // w i+1/2,j,k
  const double uFront = (u(i,  j,k+1) + u(i,j,k))   / 2.; // u i,j,k+1/2 (i+1/2 to k+1/2)
  return (wFront*uFront - wBack*uBack) / dz();  // wrong order!
}

//! compute the 1st derivative ∂ (vw) / ∂y on the w points
double CentralDifferences::computeDvwDy(int i, int j, int k) const 
{
  const double vTop    = (v(i,j,  k+1) + v(i,j,  k)) / 2.; // v i,j,k+1/2
  const double wTop    = (w(i,j+1,  k) + w(i,j,  k)) / 2.; // w i,j+1/2,k (1/2 from i to j)
  const double vBottom = (v(i,j-1,k+1) + v(i,j-1,k)) / 2.; // v i,j-1,k+1/2
  //2nd was w!!!
  const double wBottom = (w(i,j-1,k)   + w(i,j,  k)) / 2.; // w i,j-1/2,k
  return (vTop*wTop - vBottom*wBottom) / dy();  // wrong order!
}

//! compute the 1st derivative ∂ (vw) / ∂z on the v points
double CentralDifferences::computeDvwDz(int i, int j, int k) const
{
  const double wBack  = (w(i,j+1,k-1) + w(i,j,k-1)) / 2.; // w i,j+1/2,k-1
  const double vBack  = (v(i,j,  k-1) + v(i,j,k))   / 2.; // v i,j,j-1/2
  const double wFront = (w(i,j+1,k)   + w(i,j,k))   / 2.; // w i,j+1/2,k
  const double vFront = (v(i,j,  k+1) + v(i,j,k))   / 2.; // v i,j,k+1/2
  return (vFront*wFront - vBack*wBack) / dz();  // wrong order!
}