#include "pressure_solver/sor.h"

SOR::SOR(std::shared_ptr<PartitionShell> partition, 
         double epsilon, int maximumNumberOfIterations, double omega) :
         PressureSolver(partition, epsilon, maximumNumberOfIterations),
         omega_(omega)
{
    // Other omegas are likely to lead to instability
    // This is user facing, so a nice debug message is more approriate,
    // then a simple assert TODO: find other instances related instances
    if(omega <= 0.0 || omega >= 2.0)
    {
        std::stringstream str;
        str << "Omega may only be 0.0 < "  << omega << " < 2.0!\n";
        throw std::out_of_range(str.str());
    }
}

// inlining for virtual methods compiles, the hope is,
// that the following code will be pasted into SOR objects (GaussSeidel respectively)
inline void SOR::step()
{
    //! works very similar to gauss-seidel, but uses a relaxation term with omega instead

    const double dx2 = discretization_->dx2();
    const double dy2 = discretization_->dy2();
    const double dz2 = discretization_->dz2();
    const double factor = (dx2 * dy2 * dz2) / (2. * (dy2*dz2+dx2*dz2+dx2*dy2));

    for(int k = 0; k < discretization_->pkN(); k++)
    {
        for(int j = 0; j < discretization_->pjN(); j++)
        {
            for(int i = 0; i < discretization_->piN(); i++)
            {
                // store all variables with short name for readibilty
                const double p_last = discretization_->p(i,j,k);
                const double p_xm   = discretization_->p(i-1,j,k);
                const double p_xp   = discretization_->p(i+1,j,k);
                const double p_ym   = discretization_->p(i,j-1,k);
                const double p_yp   = discretization_->p(i,j+1,k);
                const double p_zm   = discretization_->p(i,j,k-1);
                const double p_zp   = discretization_->p(i,j,k+1);
                const double rhs    = discretization_->rhs(i,j,k);

                const double p_corretion = factor * ((p_xm+p_xp)/dx2 + (p_ym+p_yp)/dy2 + (p_zm+p_zp)/dz2 - rhs) - p_last;

                const double p_new = p_last + omega_ * p_corretion;
                discretization_->p(i,j,k) = p_new;
            }
        }
    }
}