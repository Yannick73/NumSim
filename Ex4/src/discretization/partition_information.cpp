#include "discretization/partition_information.h"

#include <vector>
#include <array>

PartitionInformation::PartitionInformation(std::array<int, 3> nCellsGlobal,
                                           std::array<double, 3> meshWidth, 
                                           int rank, int nRanks) :
                                           rank_(rank), nRanks_(nRanks),
                                           nCellsGlobal_(nCellsGlobal),
                                           meshWidth_(meshWidth),
                                           totalNoOfCellsGlobal_(nCellsGlobal[0]*nCellsGlobal[1]*nCellsGlobal[2])
{
    assert(rank < nRanks);
    assert(rank >= 0);

    //calculate possible partitions
    std::vector<std::array<int, 3>> partitionings = {};  //stores partitionings with format {x,y,z} 
    for(int x=1; x<=nRanks; x++) {
        for(int y=1; y<=nRanks; y++) {
            for(int z=1; z<=nRanks; z++) {
                if(x*y*z == nRanks) {
                    partitionings.push_back({x,y,z});
                }
            }
        }
    }

    //determine best partition
    std::array<int, 3> bestPartitioning = {1,1,1};
    int partitioningSurface = nCellsGlobal[0]*nCellsGlobal[1]*nCellsGlobal[2];
    for(int i=0; i<partitionings.size(); i++) {
        int iPartitioningSurface = (nCellsGlobal[0]/partitionings[i][0]) * (nCellsGlobal[1]/partitionings[i][1])  
        *(nCellsGlobal[2]/partitionings[i][2]);
        if(iPartitioningSurface < partitioningSurface) {
            bestPartitioning = partitionings[i];
            partitioningSurface = iPartitioningSurface;
        }
    }
    


    //

    partPosX_ = 0;
    partPosY_ = 0;
    partPosZ_ = 0;
    
    nodeOffset_[0] = 0;
    nodeOffset_[1] = 0;
    nodeOffset_[2] = 0;

    // it may be possible, that on the outer edge, the partition would be bigger, so use the smaller one
    nCellsLocal_[0] = nCellsGlobal[0];
    nCellsLocal_[1] = nCellsGlobal[1];
    nCellsLocal_[2] = nCellsGlobal[2];

    #ifdef GEOMETRY
    int neighbourEdgeCells = 0;
    if(!ownPartitionContainsBottomBoundary())
        neighbourEdgeCells += nCellsLocal_[0];
    if(!ownPartitionContainsTopBoundary())
        neighbourEdgeCells += nCellsLocal_[0];
    if(!ownPartitionContainsLeftBoundary())
        neighbourEdgeCells += nCellsLocal_[1];
    if(!ownPartitionContainsRightBoundary())
        neighbourEdgeCells += nCellsLocal_[1];
    int allNeighbourCells;
    MPI_Allreduce(&neighbourEdgeCells, &allNeighbourCells, 1, MPI_INT32_T, MPI_SUM, MPI_COMM_WORLD);
    if(rank == 0)
    {

	    double area = nCellsLocal_[0] * nCellsLocal_[1];
	    double averageEdges = ((double)allNeighbourCells) / ((double)(nPartX*nPartY));
	    std::cout << "PARTITION Geometry " << nCellsLocal_[0] << 'x' << nCellsLocal_[1] << " A: "
                  << area << " E: " << averageEdges << " E/A: " << averageEdges/area << std::endl;
    }
    #endif
    

    // print debugging information about the partition
    /*std::cout << "R:" << rank_ << "\t@(" << partPosX_ << ',' << partPosY_ << ")"
        << "\tneighbours N:" << northBoundaryRank_ << ", E:" << eastBoundaryRank_
        << ", S:" << southBoundaryRank_ << ", W:" << westBoundaryRank_ << "\tpartition domain: ["
        << nodeOffset_[0] << ", " << nodeOffset_[0]+nCellsLocal_[0]-1
        << "]x[" << nodeOffset_[1] << ", " << nodeOffset_[1]+nCellsLocal_[1]-1 << "]\n";*/
}
