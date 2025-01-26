#include "storage/field_variable.h"

#define EPS 1e-8

FieldVariable::FieldVariable(std::array<int, 3> size,
                             std::array<double, 3> origin,
                             std::array<double, 3> meshWidth,
                             std::string name)
    : Array3D(size, name), origin_(origin), meshWidth_(meshWidth) { }

//! only works on the assumption, that vertical is already exact or very close
/*double FieldVariable::horizontalInterpolation(double x, double y)
{
    const double i_x = (x + origin_[0]) / meshWidth_[0] - 1.0;
    const double i_y = (y + origin_[1]) / meshWidth_[1] - 1.0;
    const int y_round = std::round(i_y);
    // std::cout << "hI: " << name_ << '(' << x << ',' << y << ')' << ": t("  << i_x << ',' << i_y << ")\n";
    #ifndef NDEBUG
    if(std::abs(y_round - i_y) > EPS)
    {
        std::cerr << "Horizontal interpolation of " << name_ 
                  << " may be imprecise, because the vertical coordinate " << y 
                  << " transformed: " << i_y << " is at least " 
                  << EPS << " away from the next grid point. "
                  "Consider using bilinear interpolation instead.\n";
    }
    #endif
    const int x_l = std::floor(i_x);
    const int x_u = std::ceil(i_x);
    
    // as can be observed, halves the number of accesses, same for math
    const double x1 = (*this)(x_l, y_round);
    const double x2 = (*this)(x_u, y_round);
    return x1 + (i_x - x_l) * (x2 - x1);
}

//! only works on the assumption, that horizontal is already exact or very close
double FieldVariable::verticalInterpolation(double x, double y)
{
    const double i_x = (x + origin_[0]) / meshWidth_[0] - 1.0;
    const int x_round = std::round(i_x);
    const double i_y = (y + origin_[1]) / meshWidth_[1] - 1.0;
    // std::cout << "vI: " << name_ << '(' << x << ',' << y << ')' << ": t("  << i_x << ',' << i_y << ")\n";
    #ifndef NDEBUG
    if(std::abs(x_round - i_x) > EPS)
    {
        std::cerr << "Horizontal interpolation of " << name_ 
                  << " may be imprecise, because the horizontal coordinate " << x
                  << " transformed: " << i_x << " is at least " 
                  << EPS << " away from the next grid point."
                  "Consider using bilinear interpolation instead.\n";
    }
    #endif
    const int y_l = std::floor(i_y);
    const int y_u = std::ceil(i_y);

    const double y1 = (*this)(x_round, y_l);
    const double y2 = (*this)(x_round, y_u);
    return y1 + (i_y - y_l) * (y2 - y1);

}

double FieldVariable::bilinearInterpolation(double x, double y)
{
    // trabsformed grid points
    const double i_x = (x + origin_[0]) / meshWidth_[0] - 1.0;
    const double i_y = (y + origin_[1]) / meshWidth_[1] - 1.0;
    // std::cout << "bI: " << name_ << '(' << x << ',' << y << ')' << ": t("  << i_x << ',' << i_y << ")\n";
    
    // minus one? for C++ zero indexing scheme
    const int x_l = std::floor(i_x);
    // With rounding errors, the ceil of x + eps may be x+1,
    // even though the "correct" result would be x (because x was already an integer)
    const int x_u = std::min(size_[0], (int)std::ceil (i_x));
    const int y_l = std::floor(i_y);
    const int y_u = std::min(size_[1], (int)std::ceil (i_y));
    // - the opposite problem of the upper index beeing smaller then 0 are unlikely due
    //   to the positive offset already beeing 1


    // (*this) looks a little wonky, but enables the usage of the defined i,j indexing
    double q11 = (*this)(x_l, y_l);
    double q21 = (*this)(x_u, y_l);
    double q12 = (*this)(x_l, y_u);
    double q22 = (*this)(x_u, y_u);

    // Bilinear interpolation on grid
    double ty1 = q11 + (q21 - q11)*(i_x - x_l);
    double ty2 = q12 + (q22 - q12)*(i_x - x_l);
   
    return ty1 + (ty2 - ty1)*(i_y - y_l);
}*/