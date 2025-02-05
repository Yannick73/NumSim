#include "storage/field_variable.h"

FieldVariable::FieldVariable(std::array<int, 3> size,
                             std::array<double, 3> origin,
                             std::array<double, 3> meshWidth,
                             std::string name)
    : Array3D(size, name), meshWidth_(meshWidth), origin_(origin) { }

double FieldVariable::yzInterpolation(double x, double y, double z)
{
    const double i_x = (x - origin_[0]) / meshWidth_[0];
    const double i_y = (y - origin_[1]) / meshWidth_[1];
    const double i_z = (z - origin_[2]) / meshWidth_[2];
    
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
    const double i_x = (x - origin_[0]) / meshWidth_[0];
    const double i_y = (y - origin_[1]) / meshWidth_[1];
    const double i_z = (z - origin_[2]) / meshWidth_[2];
    
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
    const double i_x = (x - origin_[0]) / meshWidth_[0];
    const double i_y = (y - origin_[1]) / meshWidth_[1];
    const double i_z = (z - origin_[2]) / meshWidth_[2];
    
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
    const double i_x = (x - origin_[0]) / meshWidth_[0];
    const double i_y = (y - origin_[1]) / meshWidth_[1];
    const double i_z = (z - origin_[2]) / meshWidth_[2];
    
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

// debugging code, which is activated by setting -DDTEST=1 during cmake configuration
#ifdef DISCRETIZATION_TEST
double &FieldVariable::operator()(int i, int j, int k)
{
    const std::size_t index = compute_index(i, j, k);
    // make assertion conditional on DEBUG mode to optimize further
    #ifndef NDEBUG
    if(i < 0 || i >= size_[0] || j < 0 || j >= size_[1] || k < 0 || k >= size_[2])
    {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        std::stringstream str;
        str << "Out-of-bound access on " << name_ << "(i,j,k): (" << i << ',' << j << ',' << k
            << "), size: (" << size_[0] << ',' << size_[1] << ',' << size_[2] << ") in R:" << rank << "\n";
        throw std::out_of_range(str.str());
    }
    #endif
    
    const double x = i * meshWidth_[0] + origin_[0];
    const double y = j * meshWidth_[1] + origin_[1];
    const double z = k * meshWidth_[2] + origin_[2];

    std::cout << name_ << '(' << i << ',' << j << ',' << k << ") @pos (" <<
        x << ',' << y << ',' << z << ") = " << data_[index] << "\n";

    return data_[index];
}

double FieldVariable::operator()(int i, int j, int k) const
{
    //std::size_t index = compute_index(i, j, k);
    const std::size_t index = compute_index(i, j, k);
    // assert that indices are in range
    #ifndef NDEBUG
    if(i < 0 || i >= size_[0] || j < 0 || j >= size_[1] || k < 0 || k >= size_[2])
    {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        std::stringstream str;
        str << "Out-of-bound access on " << name_ << "(i,j,k): (" << i << ',' << j << ',' << k
            << "), size: (" << size_[0] << ',' << size_[1] << ',' << size_[2] << ") in R:" << rank << "\n";
        throw std::out_of_range(str.str());
    }
    #endif
    
    const double x = i * meshWidth_[0] + origin_[0];
    const double y = j * meshWidth_[1] + origin_[1];
    const double z = k * meshWidth_[2] + origin_[2];

    std::cout << name_ << '(' << i << ',' << j << ',' << k << ") @pos (" <<
        x << ',' << y << ',' << z << ") = " << data_[index] << "\n";

    return data_[index];
}
#endif
