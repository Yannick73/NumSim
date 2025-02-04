#include <iostream>
#include <mpi.h>
#include <array>
#include "storage/field_variable.h"
#include "discretization/partition_information.h"
#include "discretization/staggered_grid.h"
#include "discretization/central_differences.h"

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

void test_discretization(int rank, int nRanks)
{
    Settings settings;
    std::array<double, 3> meshWidth = {settings.physicalSize[0] / settings.nCells[0],
                                             settings.physicalSize[1] / settings.nCells[1],
                                             settings.physicalSize[2] / settings.nCells[2]};
    PartitionInformation pi(settings.nCells, meshWidth, rank, nRanks);
    CentralDifferences cd(pi, settings);
    // in this very first test, we are only really interested in the positions
    cd.computeDuvDx(2,2,2);
    cd.computeDuvDy(2,2,2);
    cd.computeDuwDx(2,2,2);
    cd.computeDuwDz(2,2,2);
    cd.computeDvwDy(2,2,2);
    cd.computeDvwDz(2,2,2);
    // here is tested, whether the indices are all in range
    cd.calculateFGH(0.1);
    cd.calculateRHS(0.1);
    cd.calculateUVW(0.1);

}

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    std::cout << "TEST!\n";
    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    //test_field_variable();
    //test_staggered_grid(world_rank, world_size);
    test_discretization(world_rank, world_size);
    return 0;
}