#include "pressure_solver/gauss_seidel.h"

inline void GaussSeidel::step()
{
    const double dx2 = discretization_->dx2();
    const double dy2 = discretization_->dy2();
    const double dz2 = discretization_->dz2();
    //TODO: 
    const double diff = (dx2 * dy2 * dz2) / (2. * (dy2*dz2+dx2*dz2+dx2*dy2));   // diffusion factor

    // algorithm updates all cells in place, as described in the lecture
    //#pragma omp simd collapse(2)
    for(int j = 0; j < discretization_->pjN(); j++)
    {
        for(int i = 0; i < discretization_->piN(); i++)
        {   
            for(int k=0; k<discretization_->pkN(); k++) {
                // store all variables with short name for readibilty
                const double p_last = discretization_->p(i,j,k);
                const double p_xm   = discretization_->p(i-1,j,k);
                const double p_xp   = discretization_->p(i+1,j,k);
                const double p_ym   = discretization_->p(i,j-1,k);
                const double p_yp   = discretization_->p(i,j+1,k);
                const double p_zm   = discretization_->p(i,j,k-1);
                const double p_zp   = discretization_->p(i,j,k+1);
                const double rhs    = discretization_->rhs(i,j,k);
                
                //TODO: 
                double p_new = diff * ((p_xm+p_xp)/dx2 + (p_ym+p_yp)/dy2 + (p_zm+p_zp)/dz2 - rhs);
                discretization_->p(i, j, k) = p_new;
            }
        }
    }
}