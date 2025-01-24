#include <iostream>
#include <mpi.h>
#include <array>
#include "storage/field_variable.h"
#include "discretization/partition_information.h"
#include "discretization/staggered_grid.h"

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

void test_staggered_grid(int rank, int nRanks)
{
    int xl = 2;
    int yl = 2;
    int zl = 2;
    std::array<int, 3> dim = {xl, yl, zl};
    std::array<double, 3> meshWidth = {0.1, 0.1, 0.1};
    PartitionInformation pi(dim, meshWidth, rank, nRanks);
    StaggeredGrid grid(pi);
    // dimensions and origin tested via visual inspection aka using gdb
    std::cout << "Staggered Grid initialized\n";
    grid.u(0,0,0) = 1;
    grid.v(0,0,0) = 2;
    grid.w(0,0,0) = 3;
    grid.f(0,0,0) = 4;
    grid.g(0,0,0) = 5;
    grid.h(0,0,0) = 6;
    grid.p(0,0,0) = 7;
    std::cout <<  "u: " << grid.u(0,0,0) << " v: " << grid.v(0,0,0) << " w: " << grid.w(0,0,0)
              << " f: " << grid.f(0,0,0) << " g: " << grid.g(0,0,0) << " h: " << grid.h(0,0,0)
              << " p: " << grid.p(0,0,0) << std::endl;
}

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    std::cout << "Hello 3D world\n";
    test_field_variable();
    test_staggered_grid(world_rank, world_size);
    return 0;
}