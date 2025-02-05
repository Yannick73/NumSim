#include <iostream>
#include <array>
#include <memory>
#include <mpi.h>
#include "storage/field_variable.h"
#include "discretization/partition_information.h"
#include "discretization/staggered_grid.h"
#include "discretization/discretization.h"
#include "discretization/central_differences.h"
#include "discretization/async_partition.h"

// setup with "cmake .. -DTEST=1"

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
    // test some of the derivatives
    /*std::cout << "\ncomputDpDx\n";
    cd.computeDpDx(2,2,2);
    std::cout << "\ncomputDpDy\n";
    cd.computeDpDy(2,2,2);
    std::cout << "\ncomputDpDz\n";
    cd.computeDpDz(2,2,2);
    std::cout << "\ncomputeD2pDx2\n";
    cd.computeD2pDx2(2,2,2);
    std::cout << "\ncomputeD2pDy2\n";
    cd.computeD2pDy2(2,2,2);
    std::cout << "\ncomputeD2pDz2\n";
    cd.computeD2pDz2(2,2,2);
    std::cout << "\ncomputeD2uDx2\n";
    cd.computeD2uDx2(2,2,2);
    std::cout << "\ncomputeD2uDy2\n";
    cd.computeD2uDy2(2,2,2);
    std::cout << "\ncomputeD2uDz2\n";
    cd.computeD2uDz2(2,2,2);
    std::cout << "\ncomputeD2vDx2\n";
    cd.computeD2vDx2(2,2,2);
    std::cout << "\ncomputeD2vDy2\n";
    cd.computeD2vDy2(2,2,2);
    std::cout << "\ncomputeD2vDz2\n";
    cd.computeD2vDz2(2,2,2);
    std::cout << "\ncomputeD2wDx2\n";
    cd.computeD2wDx2(2,2,2);
    std::cout << "\ncomputeD2wDy2\n";
    cd.computeD2wDy2(2,2,2);
    std::cout << "\ncomputeD2wDz2\n";
    cd.computeD2wDz2(2,2,2);*/
    // in this very first test, we are only really interested in the positions
    std::cout << cd.computeDu2Dx(2,2,2) << std::endl;
    std::cout << cd.computeDv2Dy(2,2,2) << std::endl;
    std::cout << cd.computeDw2Dz(2,2,2) << std::endl;
    std::cout << cd.computeDuvDx(2,2,2) << std::endl;
    std::cout << cd.computeDuvDy(2,2,2) << std::endl;
    std::cout << cd.computeDuwDx(2,2,2) << std::endl;
    std::cout << cd.computeDuwDz(2,2,2) << std::endl;
    std::cout << cd.computeDvwDy(2,2,2) << std::endl;
    std::cout << cd.computeDvwDz(2,2,2) << std::endl;
    // here is tested, whether the indices are all in range
    /*cd.calculateFGH(0.1);
    cd.calculateRHS(0.1);
    cd.calculateUVW(0.1);*/
}

void test_boundaries(int rank, int nRanks)
{
    Settings settings;
    settings.dirichletBcTop = {0.,0.,0.};
    std::array<double, 3> meshWidth = {settings.physicalSize[0] / settings.nCells[0],
                                             settings.physicalSize[1] / settings.nCells[1],
                                             settings.physicalSize[2] / settings.nCells[2]};
    PartitionInformation pi(settings.nCells, meshWidth, rank, nRanks);
    CentralDifferences cd(pi, settings);
    std::shared_ptr<Discretization> ps(&cd);
    AsyncPartition partition(ps, settings, pi);

    #pragma omp simd collapse(3)
    for(int k = 0; k < cd.ukN(); k++)
    {
        for(int j = 0; j < cd.ujN(); j++)
        {
            for(int i = 0; i < cd.uiN(); i++)
            {
                cd.u(i,j,k) = 1;
                cd.f(i,j,k) = 4;
            }
        }
    }

    // calculate v
    #pragma omp simd collapse(3)
    for(int k = 0; k < cd.vkN(); k++)
    {
        for(int j = 0; j < cd.vjN(); j++)
        {
            for(int i = 0; i < cd.viN(); i++)
            {
                cd.v(i,j,k) = 2;
                cd.g(i,j,k) = 5;
            }
        }
    }

    // calculate w
    #pragma omp simd collapse(3)
    for(int k = 0; k < cd.wkN(); k++)
    {
        for(int j = 0; j < cd.wjN(); j++)
        {
            for(int i = 0; i < cd.wiN(); i++)
            {
                cd.w(i,j,k) = 3;
                cd.h(i,j,k) = 6;
            }
        }
    }

    #pragma omp simd collapse(3)
    for(int k = 0; k < cd.pkN(); k++)
    {
        for(int j = 0; j < cd.pjN(); j++)
        {
            for(int i = 0; i < cd.piN(); i++)
            {
                cd.p(i,j,k) = 7;
            }
        }
    }
    
    partition.setBoundaryUVW();
    partition.setBoundaryFGH();
    partition.setBoundaryP();
    std::cout << "Set boundaries\n";
}

void test_partitioning()
{
    std::array<int, 3> dim = {20, 20, 20};
    std::array<double, 3> meshWidth = {0.1, 0.1, 0.1};
    PartitionInformation pi1(dim, meshWidth, 0, 1);
    PartitionInformation pi2_0(dim, meshWidth, 0, 2);
    PartitionInformation pi2_1(dim, meshWidth, 1, 2);
    PartitionInformation pi4_0(dim, meshWidth, 0, 4);
    PartitionInformation pi4_1(dim, meshWidth, 1, 4);
    PartitionInformation pi4_2(dim, meshWidth, 2, 4);
    PartitionInformation pi4_3(dim, meshWidth, 3, 4);
    PartitionInformation pi6_0(dim, meshWidth, 0, 6);
    PartitionInformation pi6_1(dim, meshWidth, 1, 6);
    PartitionInformation pi6_2(dim, meshWidth, 2, 6);
    PartitionInformation pi6_3(dim, meshWidth, 3, 6);
    PartitionInformation pi6_4(dim, meshWidth, 4, 6);
    PartitionInformation pi6_5(dim, meshWidth, 5, 6);
    PartitionInformation pi8_0(dim, meshWidth, 0, 8);
    PartitionInformation pi8_1(dim, meshWidth, 1, 8);
    PartitionInformation pi8_2(dim, meshWidth, 2, 8);
    PartitionInformation pi8_3(dim, meshWidth, 3, 8);
    PartitionInformation pi8_4(dim, meshWidth, 4, 8);
    PartitionInformation pi8_5(dim, meshWidth, 5, 8);
    PartitionInformation pi8_6(dim, meshWidth, 6, 8);
    PartitionInformation pi8_7(dim, meshWidth, 7, 8);
    for(int i = 0; i < 27; i++)
        PartitionInformation pi27_i(dim, meshWidth, i, 27);
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
    //test_discretization(world_rank, world_size);
    //test_boundaries(world_rank, world_size);
    test_partitioning();
    return 0;
}