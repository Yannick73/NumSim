#include "storage/field_variable.h"
#include <iostream>

FieldVariable::FieldVariable(std::array<int, 2> size,
                             std::array<double, 2> origin,
                             std::array<double, 2> meshWidth)
    : Array2D(size), origin_(origin), meshWidth_(meshWidth){}

double FieldVariable::interpolateAt(double x, double y) {

    // bilinear transformation seems to have a bug
    const double i_x = (x + origin_[0]) / meshWidth_[0];
    const double i_y = (y + origin_[1]) / meshWidth_[1];
    // 2 for double boundary
    const int x_l = std::floor(i_x);
    const int x_u = std::ceil(i_x);
    const int y_l = std::floor(i_y);
    const int y_u = std::ceil(i_y);

    double q11 = (*this)(x_l, y_l);
    double q21 = (*this)(x_u, y_l);
    double q12 = (*this)(x_l, y_u);
    double q22 = (*this)(x_u, y_u);

    double ty1 = q11 + (q21 - q11)*(i_x - x_l);
    double ty2 = q12 + (q22 - q12)*(i_x - x_l);

    // if it is on the right edge, interpolate between the left points
    if(i_x > size_[0]) {
        ty1 = q11;
        ty2 = q12;
    }

    if(i_y <= size_[1]) {
        return ty1 + (ty2 - ty1)*(i_y - y_l);
    }
    // if the upper edge is outside, use the lower interpolation
    // this also neatly combines the corner case, which was already set on q11
    else {
        return ty1;
    }

    // normal linear transformation (maybe this could be set in the settings?)
    /*const double i_x = (x + origin_[0]) / meshWidth_[0];
    const double i_y = (y + origin_[1]) / meshWidth_[1];
    
    const int x_l = std::floor(i_x);
    const int x_u = std::ceil(i_x);
    const int y_l = std::floor(i_y);
    const int y_u = std::ceil(i_y);

    // call the overloaded indexing of the self object
    const double support_val = (*this)(x_l, y_l);
    const double span_x = (*this)(x_u, y_l) - (*this)(x_l, y_l);
    const double span_y = (*this)(x_l, y_u) - (*this)(x_l, y_l);
    
    // calculate the differences
    const double a_x = (i_x - x_l);
    const double a_y = (i_y - y_l);

    return support_val + a_x*span_x + a_y*span_y;*/
}