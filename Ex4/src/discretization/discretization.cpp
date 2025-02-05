#include "discretization/discretization.h"

Discretization::Discretization(PartitionInformation &pi, const Settings &settings) : 
                               StaggeredGrid(pi),
                               settings_(settings),
                               meshWidth2_({pi.meshWidth()[0]*pi.meshWidth()[0],
                                            pi.meshWidth()[1]*pi.meshWidth()[1],
                                            pi.meshWidth()[2]*pi.meshWidth()[2]})
{
  assert(pi.nCellsLocal()[0] > 0 && pi.nCellsLocal()[1] > 0 && pi.nCellsLocal()[2] > 0);
  assert(pi.meshWidth()[0]   > 0 && pi.meshWidth()[1]   > 0 && pi.meshWidth()[2]   > 0);
}

//! calculate the deltaT depending on the velocity and mesh-width, includes boundaries
double Discretization::calculateVelocityDelta() const
{
  // TODO optimize: reset during calculateUV and update during Boundary->setUV?
  double max_u = std::numeric_limits<double>::epsilon();
  double max_v = std::numeric_limits<double>::epsilon();
  double max_w = std::numeric_limits<double>::epsilon();
  // find the minimal delta t due to the fluid velocity
  // velocityDelta may be deactivated in settings by setting "disableAdaptiveDt=true"
  #pragma omp simd collapse(3) reduction(max:max_u)
  for(int k = -1; k < ukN()+1; k++)
  {
    for(int j = -1; j < ujN()+1; j++)
    {
        for(int i = -1; i < uiN()+1; i++)
        {
            max_u = std::max(max_u, u(i,j,k));
        }
    }
  }
  #pragma omp simd collapse(3) reduction(max:max_v)
  for(int k = -1; k < vkN()+1; k++)
  {
    for(int j = -1; j < vjN()+1; j++)
    {
        for(int i = -1; i < viN()+1; i++)
        {
            max_v = std::max(max_v, v(i,j,k));
        }
    }
  }
  #pragma omp simd collapse(3) reduction(max:max_v)
  for(int k = -1; k < wkN()+1; k++)
  {
    for(int j = -1; j < wjN()+1; j++)
    {
        for(int i = -1; i < wiN()+1; i++)
        {
          // max_v instead of max_w
            max_w = std::max(max_w, w(i,j,k));
        }
    }
  }
  std::cout << "max-u " << max_u << "\tmax-v " << max_v << "\tmax-w " << max_w << std::endl;
  double dTx = dx() / max_u;
  double dTy = dy() / max_v;
  double dTz = dz() / max_w;
  return settings_.tau * std::min(dTx, std::min(dTy, dTz));
}

double Discretization::calculateReynoldsDelta() const
{
  // same method to arrive to this result, but 3 coordinates instead of 2
  const double result = settings_.tau * (settings_.re/2) * (dx2()*dy2()*dz2()) / 
    (dy2()*dz2() + dx2()*dz2() + dx2()*dy2());
  return result;
}

//! calculate the preliminary velocities (f, g) and the rhs for pressure solver
void Discretization::calculateFGH(double deltaT)
{
  // calculate f
  #pragma omp simd collapse(3)
  for(int k = 0; k < ukN(); k++)
  {
    for(int j = 0; j < ujN(); j++)
    {
        for(int i = 0; i < uiN(); i++)
        {
            const double D2uDx2 = computeD2uDx2(i,j,k);
            const double D2uDy2 = computeD2uDy2(i,j,k);
            const double D2uDz2 = computeD2uDz2(i,j,k);
            const double Du2Dx  = computeDu2Dx (i,j,k);
            const double DuvDy  = computeDuvDy (i,j,k);
            const double DuwDz  = computeDuwDz (i,j,k);
            f(i,j,k) = u(i,j,k) + deltaT*((D2uDx2+D2uDy2+D2uDz2)/settings_.re 
              - Du2Dx - DuvDy - DuwDz + settings_.g[0]);
        }
    }
  }

  // calculate g
  #pragma omp simd collapse(3)
  for(int k = 0; k < vkN(); k++)
  {
    for(int j = 0; j < vjN(); j++)
    {
        for(int i = 0; i < viN(); i++)
        {
            const double D2vDx2 = computeD2vDx2(i,j,k);
            const double D2vDy2 = computeD2vDy2(i,j,k);
            const double D2vDz2 = computeD2vDz2(i,j,k);
            const double DuvDx  = computeDuvDx (i,j,k);
            const double Dv2Dy  = computeDv2Dy (i,j,k);
            const double DvwDz  = computeDvwDz (i,j,k);
            g(i,j,k) = v(i,j,k) + deltaT*((D2vDx2+D2vDy2+D2vDz2)/settings_.re
              - DuvDx- Dv2Dy - DvwDz + settings_.g[1]);
        }
    }
  }

  // calculate h
  #pragma omp simd collapse(3)
  for(int k = 0; k < wkN(); k++)
  {
    for(int j = 0; j < wjN(); j++)
    {
        for(int i = 0; i < wiN(); i++)
        {
            const double D2wDx2 = computeD2wDx2(i,j,k);
            const double D2wDy2 = computeD2wDy2(i,j,k);
            const double D2wDz2 = computeD2wDz2(i,j,k);
            const double DuwDx  = computeDuwDx (i,j,k);
            const double DvwDy  = computeDvwDy (i,j,k);
            const double Dw2Dz  = computeDw2Dz (i,j,k);
            // h!!! instead of g and w instead of v
            h(i,j,k) = w(i,j,k) + deltaT*((D2wDx2+D2wDy2+D2wDz2)/settings_.re
              - DuwDx - DvwDy - Dw2Dz + settings_.g[2]);
        }
    }
  }

}

void Discretization::calculateRHS(double deltaT)
{
  // using f and g, calculate rhs
  #pragma omp simd collapse(3)
  for(int k = 0; k < pkN(); k++)
  {
    for(int j = 0; j < pjN(); j++)
    {
      for(int i = 0; i < piN(); i++)
      {
        const double DfDx = (f(i,j,k) - f(i-1,j,  k))   / dx();
        const double DgDy = (g(i,j,k) - g(i,  j-1,k))   / dy();
        const double DhDz = (h(i,j,k) - h(i,  j,  k-1)) / dz();
        rhs(i,j,k) = (DfDx + DgDy + DhDz) / deltaT;
      }
    }
  }
}

//! calculate the final velocities using f, g and p
void Discretization::calculateUVW(double deltaT)
{
  // calculate u
  #pragma omp simd collapse(3)
  for(int k = 0; k < ukN(); k++)
  {
    for(int j = 0; j < ujN(); j++)
    {
      for(int i = 0; i < uiN(); i++)
      {
        u(i,j,k) = f(i,j,k) - deltaT*computeDpDx(i,j,k);
      }
    }
  }

  // calculate v
  #pragma omp simd collapse(3)
  for(int k = 0; k < vkN(); k++)
  {
    for(int j = 0; j < vjN(); j++)
    {
      for(int i = 0; i < viN(); i++)
      {
        v(i,j,k) = g(i,j,k) - deltaT*computeDpDy(i,j,k);
      }
    }
  }

  // calculate w
  #pragma omp simd collapse(3)
  for(int k = 0; k < wkN(); k++)
  {
    for(int j = 0; j < wjN(); j++)
    {
      for(int i = 0; i < wiN(); i++)
      {
        w(i,j,k) = h(i,j,k) - deltaT*computeDpDz(i,j,k);
      }
    }
  }
}
