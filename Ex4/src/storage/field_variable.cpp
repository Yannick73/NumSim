#include "storage/field_variable.h"

FieldVariable::FieldVariable(std::array<int, 3> size,
                             std::array<double, 3> origin,
                             std::array<double, 3> meshWidth,
                             std::string name)
    : Array3D(size, name), origin_(origin), meshWidth_(meshWidth) { }

double FieldVariable::yzInterpolation(double x, double y, double z)
{
    const double i_x = (x + origin_[0]) / meshWidth_[0] - 1.0;
    const double i_y = (y + origin_[1]) / meshWidth_[1] - 1.0;
    const double i_z = (z + origin_[2]) / meshWidth_[2] - 1.0;
    
    const int x_round = std::round(i_x);
    const int y_upper = std::ceil (i_y);
    const int y_lower = std::floor(i_y);
    const int z_upper = std::ceil (i_z);
    const int z_lower = std::floor(i_z);

    const double q00 = (*this)(x_round, y_lower, z_lower);
    const double q01 = (*this)(x_round, y_lower, z_upper);
    const double q10 = (*this)(x_round, y_upper, z_lower);
    const double q11 = (*this)(x_round, y_upper, z_upper);

    return (q00 + q01 + q10 + q11) / 4.0;
}

double FieldVariable::xzInterpolation(double x, double y, double z)
{
    const double i_x = (x + origin_[0]) / meshWidth_[0] - 1.0;
    const double i_y = (y + origin_[1]) / meshWidth_[1] - 1.0;
    const double i_z = (z + origin_[2]) / meshWidth_[2] - 1.0;
    
    const int y_round = std::round(i_y);
    const int x_upper = std::ceil (i_x);
    const int x_lower = std::floor(i_x);
    const int z_upper = std::ceil (i_z);
    const int z_lower = std::floor(i_z);

    const double q00 = (*this)(x_lower, y_round, z_lower);
    const double q01 = (*this)(x_lower, y_round, z_upper);
    const double q10 = (*this)(x_upper, y_round, z_lower);
    const double q11 = (*this)(x_upper, y_round, z_upper);

    return (q00 + q01 + q10 + q11) / 4.0;
}

double FieldVariable::xyInterpolation(double x, double y, double z)
{
    const double i_x = (x + origin_[0]) / meshWidth_[0] - 1.0;
    const double i_y = (y + origin_[1]) / meshWidth_[1] - 1.0;
    const double i_z = (z + origin_[2]) / meshWidth_[2] - 1.0;
    
    const int z_round = std::round(i_z);
    const int x_upper = std::ceil (i_x);
    const int x_lower = std::floor(i_x);
    const int y_upper = std::ceil (i_y);
    const int y_lower = std::floor(i_y);

    const double q00 = (*this)(x_lower, y_lower, z_round);
    const double q01 = (*this)(x_lower, y_upper, z_round);
    const double q10 = (*this)(x_upper, y_lower, z_round);
    const double q11 = (*this)(x_upper, y_upper, z_round);

    return (q00 + q01 + q10 + q11) / 4.0;
}

// this here could be trilinear interpolation, but it is a simplification,
// because the mesh is regular and the points are always in the same mid positions
double FieldVariable::midInterpolation(double x, double y, double z)
{
    const double i_x = (x + origin_[0]) / meshWidth_[0] - 1.0;
    const double i_y = (y + origin_[1]) / meshWidth_[1] - 1.0;
    const double i_z = (z + origin_[2]) / meshWidth_[2] - 1.0;
    
    const int x_upper = std::ceil (i_x);
    const int x_lower = std::floor(i_x);
    const int y_upper = std::ceil (i_y);
    const int y_lower = std::floor(i_y);
    const int z_upper = std::ceil (i_z);
    const int z_lower = std::floor(i_z);

    const double q000 = (*this)(x_lower, y_lower, z_lower);
    const double q001 = (*this)(x_lower, y_lower, z_upper);
    const double q010 = (*this)(x_lower, y_upper, z_lower);
    const double q011 = (*this)(x_lower, y_upper, z_upper);
    const double q100 = (*this)(x_upper, y_lower, z_lower);
    const double q101 = (*this)(x_upper, y_lower, z_upper);
    const double q110 = (*this)(x_upper, y_upper, z_lower);
    const double q111 = (*this)(x_upper, y_upper, z_upper);

    return (q000 + q001 + q010 + q011 + q100 + q101 + q110 + q111) / 8.0;
}
