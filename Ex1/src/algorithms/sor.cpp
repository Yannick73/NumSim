#include "sor.h"

SOR::SOR(const std::shared_ptr<Discretization> discretization, 
         double epsilon, int maximumNumberOfIterations, double omega) :
         PressureSolver(discretization, epsilon, maximumNumberOfIterations),
         omega_(omega)
{
    // If Omega was 0 the solver would not update, so it is disallowed
    assert(omega > 0);
}

// inlining for virtual methods compiles, the hope is,
// that the following code will be pasted into SOR objects (GaussSeidel respectively)
inline double SOR::step()
{
    //! works very similar to gauss-seidel, but uses a relaxation term with omega instead
    //double res = 0;     // residuum uses Euclidian norm (with p beeing treated as a vector)

    const double dx2 = discretization_->dx2();
    const double dy2 = discretization_->dy2();
    const double diff = (dx2 * dy2) / (2. * (dx2 + dy2));   // diffusion factor

    for(int j = 0; j < discretization_->jEnd(); j++)
    {
        for(int i = 0; i < discretization_->iEnd(); i++)
        {
            // store all variables with short name for readibilty
            const double p_last = discretization_->p(i,j);
            const double p_xm   = discretization_->p(i-1,j);
            const double p_xp   = discretization_->p(i+1,j);
            const double p_ym   = discretization_->p(i,j-1);
            const double p_yp   = discretization_->p(i,j+1);
            const double rhs    = discretization_->rhs(i,j);

            const double p_corretion = diff * ((p_xm+p_xp)/dx2 + (p_ym+p_yp)/dy2 - rhs) - p_last;

            const double p_new = p_last + omega_ * p_corretion;
            discretization_->p(i, j) = p_new;
            
            double e_i = std::fabs(p_last - p_new);
            //res += e_i * e_i;
        }
    }    
    //return std::sqrt(res);
    return calculateResiduum();
}