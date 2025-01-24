#include "pressure_solver/checkerboard.h"

Checkerboard::Checkerboard(std::shared_ptr<PartitionShell> partition, 
         double epsilon, int maximumNumberOfIterations, double omega) :
         PressureSolver(partition, epsilon, maximumNumberOfIterations),
         omega_(omega)
{
    // Other omegas are likely to lead to instability
    if(omega <= 0.0 || omega >= 2.0)
    {
        std::stringstream str;
        str << "Omega may only be 0.0 < "  << omega << " < 2.0!\n";
        throw std::out_of_range(str.str());
    }
}

//! Works similar to SOR, but in two steps
inline void Checkerboard::step()
{
    const double dx2 = discretization_->dx2();
    const double dy2 = discretization_->dy2();
    const double factor = (dx2 * dy2) / (2. * (dx2 + dy2));
    const int nodeOffset = partition_->pi_.nodeOffset()[0] + partition_->pi_.nodeOffset()[1];
    for(int j = 0; j < discretization_->pjN(); j++)
    {
        const int offset = j & 0b1; // same as % 2, but possibly more efficient
        #pragma omp simd
        for(int i = offset; i < discretization_->piN(); i += 2)
        {
            // store all variables with short name for readibilty
            const double p_last = discretization_->p(i,j);
            const double p_xm   = discretization_->p(i-1,j);
            const double p_xp   = discretization_->p(i+1,j);
            const double p_ym   = discretization_->p(i,j-1);
            const double p_yp   = discretization_->p(i,j+1);
            const double rhs    = discretization_->rhs(i,j);

            const double p_corretion = factor * ((p_xm+p_xp)/dx2 + (p_ym+p_yp)/dy2 - rhs) - p_last;

            const double p_new = p_last + omega_ * p_corretion;
            discretization_->p(i, j) = p_new;
        }
    }
    partition_->exchangeP();

    for(int j = 0; j < discretization_->pjN(); j++)
    {
        const int offset = (~j) & 0b1;  // now the LSB is flipped, compared to the first loop
        #pragma omp simd
        for(int i = offset; i < discretization_->piN(); i += 2)
        {
            const double p_last = discretization_->p(i,j);
            const double p_xm   = discretization_->p(i-1,j);
            const double p_xp   = discretization_->p(i+1,j);
            const double p_ym   = discretization_->p(i,j-1);
            const double p_yp   = discretization_->p(i,j+1);
            const double rhs    = discretization_->rhs(i,j);

            const double p_corretion = factor * ((p_xm+p_xp)/dx2 + (p_ym+p_yp)/dy2 - rhs) - p_last;

            const double p_new = p_last + omega_ * p_corretion;
            discretization_->p(i, j) = p_new;
        }
    }
}