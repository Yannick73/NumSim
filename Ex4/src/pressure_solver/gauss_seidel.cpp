#include "pressure_solver/gauss_seidel.h"

inline void GaussSeidel::step()
{
    const double dx2 = discretization_->dx2();
    const double dy2 = discretization_->dy2();
    const double diff = (dx2 * dy2) / (2. * (dx2 + dy2));   // diffusion factor

    // algorithm updates all cells in place, as described in the lecture
    #pragma omp simd collapse(2)
    for(int j = 0; j < discretization_->pjN(); j++)
    {
        for(int i = 0; i < discretization_->piN(); i++)
        {
            // store all variables with short name for readibilty
            const double p_last = discretization_->p(i,j);
            const double p_xm   = discretization_->p(i-1,j);
            const double p_xp   = discretization_->p(i+1,j);
            const double p_ym   = discretization_->p(i,j-1);
            const double p_yp   = discretization_->p(i,j+1);
            const double rhs    = discretization_->rhs(i,j);
            
            double p_new = diff * ((p_xm+p_xp)/dx2 + (p_ym+p_yp)/dy2 - rhs);
            discretization_->p(i, j) = p_new;
        }
    }
}