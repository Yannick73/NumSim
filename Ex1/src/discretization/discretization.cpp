#include "discretization/discretization.h"

Discretization::Discretization(std::array<int, 2> nCells,
                               std::array<double, 2> meshWidth,
                               Settings settings) : 
                               StaggeredGrid(nCells, meshWidth),
                               settings_(settings)
{
  assert(nCells[0]    > 0 && nCells[1]    > 0);
  assert(meshWidth[0] > 0 && meshWidth[1] > 0);
}

//! calculate the deltaT depending on the velocity and mesh-width, includes boundaries
double Discretization::calculateVelocityDelta() const
{
  double max_u = std::numeric_limits<double>::epsilon();
  double max_v = std::numeric_limits<double>::epsilon();
  // find the minimal delta t due to the fluid velocity
  // velocityDelta may be deactivated in settings by setting "disableAdaptiveDt=true"
  for(int j = -1; j < ujEnd()+1; j++)
  {
      for(int i = -1; i < uiEnd()+1; i++)
      {
          max_u = std::max(max_u, u(i,j));
      }
  }
  for(int j = -1; j < vjEnd()+1; j++)
  {
      for(int i = -1; i < viEnd()+1; i++)
      {
          max_v = std::max(max_v, v(i,j));
      }
  }
  //std::cout << "max-vel \t" << max_vel << std::endl;
  double dTx = dx() / max_u;
  double dTy = dy() / max_v;
  return settings_.tau * std::min(dTx, dTy);
}

//! calculate the preliminary velocities (f, g) and the rhs for pressure solver
void Discretization::calculateFG_RHS(double deltaT)
{
  // calculate f
  for(int j = 0; j < ujEnd(); j++)
  {
      for(int i = 0; i < uiEnd(); i++)
      {
          const double D2uDx2 = computeD2uDx2(i,j);
          const double D2uDy2 = computeD2uDy2(i,j);
          const double Du2Dx  = computeDu2Dx(i,j);
          const double DuvDy  = computeDuvDy(i,j);
          f(i,j) = u(i,j) + deltaT*((D2uDx2+D2uDy2)/settings_.re - Du2Dx - DuvDy + settings_.g[0]);
      }
  }

  // calculate g
  for(int j = 0; j < vjEnd(); j++)
  {
      for(int i = 0; i < viEnd(); i++)
      {
          const double D2vDx2 = computeD2vDx2(i,j);
          const double D2vDy2 = computeD2vDy2(i,j);
          const double Dv2Dy  = computeDv2Dy(i,j);
          const double DuvDx  = computeDuvDx(i,j);
          g(i,j) = v(i,j) + deltaT*((D2vDx2+D2vDy2)/settings_.re - Dv2Dy - DuvDx + settings_.g[1]);
      }
  }

  // using f and g, calculate rhs
  for(int j = 0; j < jEnd(); j++)
  {
      for(int i = 0; i < iEnd(); i++)
      {
          const double DfDx = (f(i,j)-f(i-1,j))/dx();
          const double DgDy = (g(i,j)-g(i,j-1))/dy();
          rhs(i,j) = (DfDx + DgDy) / deltaT;
      }
  }
}

//! calculate the final velocities using f, g and p
void Discretization::calculateUV(double deltaT)
{
  // calculate u
  for(int j = 0; j < ujEnd(); j++)
  {
      for(int i = 0; i < uiEnd(); i++)
      {
          u(i,j) = f(i,j) - deltaT*computeDpDx(i,j);
      }
  }

  // calculate v
  for(int j = 0; j < vjEnd(); j++)
  {
      for(int i = 0; i < viEnd(); i++)
      {
          v(i,j) = g(i,j)-deltaT*computeDpDy(i,j);
      }
  }
}

//! compute the 1st derivative ∂p / ∂x
double Discretization::computeDpDx(int i, int j) const
{
  return (p(i+1,j) - p(i,j)) / dx();
}

//! compute the 1st derivative ∂p / ∂y
double Discretization::computeDpDy(int i, int j) const
{
  return (p(i,j+1) - p(i,j)) / dy();
}

//! compute the 2st derivative ∂^2p / ∂x^2
double Discretization::computeD2pDx2(int i, int j) const
{
  return (p(i-1,j) - 2*p(i,j) + p(i+1,j)) / dx2();
}

//! compute the 2st derivative ∂^2p / ∂y^2
double Discretization::computeD2pDy2(int i, int j) const
{
  return (p(i,j-1) - 2*p(i,j) + p(i,j+1)) / dy2();
}

//! compute the 2nd derivative ∂^2 u / ∂x^2
double Discretization::computeD2uDx2(int i, int j) const 
{
  return (u(i+1,j) - 2*u(i,j) + u(i-1,j)) / dx2();
}

//! compute the 2nd derivative ∂^2 u / ∂y^2
double Discretization::computeD2uDy2(int i, int j) const 
{
  return (u(i,j+1) - 2*u(i,j) + u(i,j-1)) / dy2();
}

//! compute the 2nd derivative ∂^2 v / ∂x^2
double Discretization::computeD2vDx2(int i, int j) const 
{
  return (v(i+1,j) - 2*v(i,j) + v(i-1,j)) / dx2();
}

//! compute the 2nd derivative ∂^2 v / ∂y^2
double Discretization::computeD2vDy2(int i, int j) const
{
  return (v(i,j+1) - 2*v(i,j) + v(i,j-1)) / dy2();
}