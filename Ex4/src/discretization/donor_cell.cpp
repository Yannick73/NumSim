#include "discretization/donor_cell.h"

//! compute the 1st derivative ∂ u^2 / ∂x
double DonorCell::computeDu2Dx(int i, int j, int k) const
{
    const double uRight = (u(i,j,k)   + u(i+1,j,k)) / 2.;
    const double uLeft  = (u(i-1,j,k) + u(i,j,k))   / 2.;
    const double centralDiff = (uRight*uRight - uLeft*uLeft) / dx();

    const double uDonorLeft  = std::abs(uLeft)  * ((u(i-1,j,k) - u(i,j,k))   / 2.);
    const double uDonorRight = std::abs(uRight) * ((u(i,j,k)   - u(i+1,j,k)) / 2.);
    return centralDiff + alpha_ * (uDonorRight - uDonorLeft) / dx();
}

//! compute the 1st derivative ∂ v^2 / ∂x
double DonorCell::computeDv2Dy(int i, int j, int k) const 
{
    const double vTop    = (v(i,j,k)   + v(i,j+1,k)) / 2.;
    const double vBottom = (v(i,j-1,k) + v(i,j,k))   / 2.;
    const double centralDiff = (vTop*vTop - vBottom*vBottom) / dy();

    const double vDonorTop    = std::abs(vTop)    * ((v(i,j,k)   - v(i,j+1,k)) / 2.);
    const double vDonorBottom = std::abs(vBottom) * ((v(i,j-1,k) - v(i,j,k))   / 2.);
    return centralDiff + alpha_ * (vDonorTop - vDonorBottom) / dy();
}

//! compute the 1st derivative ∂ v^2 / ∂x
double DonorCell::computeDw2Dz(int i, int j, int k) const 
{
    const double wFront = (w(i,j,k)   + w(i,j,k+1)) / 2.;
    const double wBack  = (w(i,j,k-1) + w(i,j,k))   / 2.;
    const double centralDiff = (wFront*wFront - wBack*wBack) / dz();

    const double wDonorFront = std::abs(wFront) * (w(i,j,k) - w(i,j,k+1)) / 2.;
    const double wDonorBack  = std::abs(wBack)  * (w(i,j,k-1) - w(i,j,k)) / 2.;
    return centralDiff + alpha_ * (wDonorFront - wDonorBack) / dz();
}

//! compute the 1st derivative ∂ (uv) / ∂x on the v points
double DonorCell::computeDuvDx(int i, int j, int k) const 
{
    const double vRight = (v(i+1,j,k)   + v(i,j,k))   / 2.; // v at top right corner of cell
    const double vLeft  = (v(i,j,k)     + v(i-1,j,k)) / 2.; // v at top left corner of cell
    const double uRight = (u(i,j+1,k)   + u(i,j,k))   / 2.; // u at top right corner of cell
    const double uLeft  = (u(i-1,j+1,k) + u(i-1,j,k)) / 2.; // u at top left corner of cell
    const double centralDiff  = (vRight*uRight - vLeft*uLeft) / dx();

    const double uvDonorRight = std::abs(uRight) * ((v(i,j,k)   - v(i+1,j,k)) / 2.);
    const double uvDonorLeft  = std::abs(uLeft)  * ((v(i-1,j,k) - v(i,j,k))   / 2.);
    return centralDiff + alpha_ * (uvDonorRight - uvDonorLeft) / dx();
}

//! compute the 1st derivative ∂ (uv) / ∂y
double DonorCell::computeDuvDy(int i, int j, int k) const 
{
    const double vTop    = (v(i,j,k)   + v(i+1,j,k))   / 2.; // v at top right corner of cell
    const double vBottom = (v(i,j-1,k) + v(i+1,j-1,k)) / 2.; // v at bottom right corner of cell
    const double uTop    = (u(i,j,k)   + u(i,j+1,k))   / 2.; // u at top right corner of cell
    const double uBottom = (u(i,j-1,k) + u(i,j,k))     / 2.; // u at bottom right corner of cell
    const double centralDiff    = (vTop*uTop - vBottom*uBottom) / dy();

    const double uvDonorTop     = std::abs(vTop)    * ((u(i, j,k)     - u(i, j + 1,k)) / 2.);
    const double uvDonorBottom  = std::abs(vBottom) * ((u(i, j - 1,k) - u(i, j,k))     / 2.);
    return centralDiff + alpha_ * (uvDonorTop - uvDonorBottom) / dy();
}

//! compute the 1st derivative ∂ (uw) / ∂x
double DonorCell::computeDuwDx(int i, int j, int k) const 
{
  const double uLeft  = (u(i-1,j,k+1) + u(i-1,j,k)) / 2.; // u i-1,j,k+1/2
  const double wLeft  = (w(i-1,j,k)   + w(i,  j,k)) / 2.; // w i-1/2,j,k
  const double uRight = (u(i,  j,k+1) + u(i,  j,k)) / 2.; // u i,j,k+1/2
  const double wRight = (w(i+1,j,k)   + w(i,  j,k)) / 2.; // w i+1/2,k
  const double centralDiff = (uRight*wRight - uLeft*wLeft) / dx();

  const double uwDonorRight = std::abs(uRight) * (w(i,j,k) - w(i+1,j,k)) / 2.; 
  const double uwDonorLeft  = std::abs(uLeft)  * (w(i-1,j,k) - w(i,j,k)) / 2.;
  return centralDiff + alpha_ * (uwDonorRight - uwDonorLeft) / dx();
}

//! compute the 1st derivative ∂ (uw) / ∂z on u points
double DonorCell::computeDuwDz(int i, int j, int k) const 
{
  // Dunno, what was broken, but now it is fixed (reduced error by 10)
  // back to front z increases
  const double wBack  = (w(i+1,j,k-1) + w(i,j,k-1)) / 2.; // w i+1/2,j,k-1
  const double uBack  = (u(i,  j,k-1) + u(i,j,k))   / 2.; // u i,j,k-1/2
  const double wFront = (w(i+1,j,k)   + w(i,j,k))   / 2.; // w i+1/2,j,k
  const double uFront = (u(i,  j,k+1) + u(i,j,k))   / 2.; // u i,j,k+1/2 (i+1/2 to k+1/2)
  const double centralDiff = (wFront*uFront - wBack*uBack) / dz();  // wrong order!
  
  const double uwDonorFront = std::abs(wFront) * (u(i,  j,k) - u(i,j,k+1))   / 2.;
  const double uwDonorBack = std::abs(wBack) * (u(i,  j,k-1) - u(i,j,k))   / 2.;
  return centralDiff + alpha_ * (uwDonorFront - uwDonorBack) / dz();
}

//! compute the 1st derivative ∂ (vw) / ∂y on the w points
double DonorCell::computeDvwDy(int i, int j, int k) const 
{
  const double vTop    = (v(i,  j,  k+1) + v(i,j,  k)) / 2.; // v i,j,k+1/2
  const double wTop    = (w(i,  j+1,k)   + w(i,j,  k)) / 2.; // w i,j+1/2,k (i instead of j)
  const double vBottom = (v(i,  j-1,k+1) + v(i,j-1,k)) / 2.; // v i,j-1,k+1/2
  const double wBottom = (w(i,  j-1,k)   + w(i,j,  k)) / 2.; // w i,j-1/2,k
  const double centralDiff = (vTop*wTop - vBottom*wBottom) / dy();

  const double vwDonorTop    = std::abs(vTop)    * (w(i,j,  k) - w(i,j+1,k)) / 2.;
  const double vwDonorBottom = std::abs(vBottom) * (w(i,j-1,k) - w(i,j,k)) / 2.;
  return centralDiff + alpha_ * (vwDonorTop - vwDonorBottom) / dy();
}

//! compute the 1st derivative ∂ (vw) / ∂z on the v points
double DonorCell::computeDvwDz(int i, int j, int k) const
{
  const double wBack  = (w(i,j+1,k-1) + w(i,j,k-1)) / 2.; // w i,j+1/2,k-1
  const double vBack  = (v(i,j,  k-1) + v(i,j,k))   / 2.; // v i,j,k-1/2
  const double wFront = (w(i,j+1,k)   + w(i,j,k))   / 2.; // w i,j+1/2,k
  const double vFront = (v(i,j,  k+1) + v(i,j,k))   / 2.; // v i,j,k+1/2
  const double centralDiff = (vFront*wFront - vBack*wBack) / dz();

  const double vwDonorFront = std::abs(wFront) * (v(i,j,k)   - v(i,j,k+1)) / 2.;
  const double vwDonorBack  = std::abs(wBack)  * (v(i,j,k-1) - v(i,j,k))   / 2.;
  return centralDiff + alpha_ * (vwDonorFront - vwDonorBack) / dz();
}