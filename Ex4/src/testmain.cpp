#include <iostream>
#include <mpi.h>
#include <array>
#include "storage/field_variable.h"

void test_field_variable()
{
    // generate field variable
    int xl = 2;
    int yl = 2;
    int zl = 2;
    std::array<int, 3> dim = {xl, yl, zl};
    std::array<double, 3> origin = {0, 0, 0};
    std::array<double, 3> meshWidth = {0.1, 0.1, 0.1};
    FieldVariable var(dim, origin, meshWidth, "test_var");
    // set some arbitrary values
    for(int z = 0; z < zl; z++) {
        for(int y = 0; y < yl; y++) {
            for(int x = 0; x < xl; x++) {
                var(x, y, z) = x+y+z;
            }
        }
    }
    // then read them back
    for(int z = 0; z < zl; z++) {
        for(int y = 0; y < yl; y++) {
            for(int x = 0; x < xl; x++) {
                std::cout << "(x,y,z) (" << x << ',' << y << ',' << z << "): "
                    << var(x, y, z) << std::endl;
            }
        }
    }
}

int main(int argc, char *argv[])
{
    std::cout << "Hello 3D world\n";
    test_field_variable();
    return 0;
}