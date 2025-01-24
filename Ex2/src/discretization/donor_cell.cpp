#include "discretization/donor_cell.h"

//! compute the 1st derivative ∂ u^2 / ∂x
double DonorCell::computeDu2Dx(int i, int j) const 
{
    const double uRight = (u(i,j)   + u(i+1,j)) / 2.;
    const double uLeft  = (u(i-1,j) + u(i,j))   / 2.;
    const double centralDiff = (uRight*uRight - uLeft*uLeft) / dx();

    const double uDonorLeft  = std::abs(uLeft)  * ((u(i-1,j) - u(i,j))   / 2.);
    const double uDonorRight = std::abs(uRight) * ((u(i,j)   - u(i+1,j)) / 2.);
    return centralDiff + alpha_ * (uDonorRight - uDonorLeft) / dx();
}

//! compute the 1st derivative ∂ v^2 / ∂x
double DonorCell::computeDv2Dy(int i, int j) const 
{
    const double vTop    = (v(i,j)   + v(i,j+1)) / 2.;
    const double vBottom = (v(i,j-1) + v(i,j))   / 2.;
    const double centralDiff = (vTop*vTop - vBottom*vBottom) / dy();

    const double vDonorTop    = std::abs(vTop)    * ((v(i,j)   - v(i,j+1)) / 2.);
    const double vDonorBottom = std::abs(vBottom) * ((v(i,j-1) - v(i,j))   / 2.);
    return centralDiff + alpha_ * (vDonorTop - vDonorBottom) / dy();
}

//! compute the 1st derivative ∂ (uv) / ∂x
double DonorCell::computeDuvDx(int i, int j) const 
{
    const double vRight = (v(i+1,j)   + v(i,j))   / 2.; // v at top right corner of cell
    const double vLeft  = (v(i,j)     + v(i-1,j)) / 2.; // v at top left corner of cell
    const double uRight = (u(i,j+1)   + u(i,j))   / 2.; // u at top right corner of cell
    const double uLeft  = (u(i-1,j+1) + u(i-1,j)) / 2.; // u at top left corner of cell
    const double centralDiff  = (vRight*uRight - vLeft*uLeft) / dx();

    const double uvDonorRight = std::abs(uRight) * ((v(i,j)   - v(i+1,j)) / 2.);
    const double uvDonorLeft  = std::abs(uLeft)  * ((v(i-1,j) - v(i,j))   / 2.);
    return centralDiff + alpha_ * (uvDonorRight - uvDonorLeft) / dx();
}

//! compute the 1st derivative ∂ (uv) / ∂y
double DonorCell::computeDuvDy(int i, int j) const 
{
    const double vTop    = (v(i,j)   + v(i+1,j))   / 2.; // v at top right corner of cell
    const double vBottom = (v(i,j-1) + v(i+1,j-1)) / 2.; // v at bottom right corner of cell
    const double uTop    = (u(i,j)   + u(i,j+1))   / 2.; // u at top right corner of cell
    const double uBottom = (u(i,j-1) + u(i,j))     / 2.; // u at bottom right corner of cell
    const double centralDiff    = (vTop*uTop - vBottom*uBottom) / dy();

    const double uvDonorTop     = std::abs(vTop)    * ((u(i, j)     - u(i, j + 1)) / 2.);
    const double uvDonorBottom  = std::abs(vBottom) * ((u(i, j - 1) - u(i, j))     / 2.);
    return centralDiff + alpha_ * (uvDonorTop - uvDonorBottom) / dy();

}
