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
    double partitioningSurface = nCellsGlobal[0]*nCellsGlobal[1]*nCellsGlobal[2];
    std::cout << nCellsGlobal[0] << ", " << nCellsGlobal[1] << ", " << nCellsGlobal[2] << std::endl;
    for(int i=0; i<partitionings.size(); i++) {
        double iPartitioningSurface = 2*((double)nCellsGlobal[0]/partitionings[i][0]) * ((double)nCellsGlobal[1]/partitionings[i][1])   //bottom and top surface
        + 2*((double)nCellsGlobal[0]/partitionings[i][0]) *((double)nCellsGlobal[2]/partitionings[i][2])  //front and back surface
        + 2*((double)nCellsGlobal[1]/partitionings[i][1]) * ((double)nCellsGlobal[2]/partitionings[i][2]);    //right and left surface
        std::cout << "partioning " << i << ": " << partitionings[i][0] << ", " << partitionings[i][1] << ", " << partitionings[i][2] << ", " << iPartitioningSurface << std::endl;
        if(iPartitioningSurface < partitioningSurface) {
            bestPartitioning = partitionings[i];
            partitioningSurface = iPartitioningSurface;
        }
    }
    std::cout << "partition: " << bestPartitioning[0] << ", " << bestPartitioning[1] << ", "
    << bestPartitioning[2] << ", "<< bestPartitioning[0]*bestPartitioning[1]*bestPartitioning[2] << std::endl;

    int nPartX = bestPartitioning[0];
    int nPartY = bestPartitioning[1];
    int nPartZ = bestPartitioning[2];

    if(nRanks != nPartX*nPartY*nPartZ)
    {
        std::stringstream str;
        str << "The found partition scheme for " << nRanks << " was " 
            << nPartX << 'x' << nPartY << 'x' << nPartZ << " which does not multiply to nRanks.\n";
        throw std::runtime_error(str.str());
    }

    if(rank_ == 0)
        std::cout << "number of ranks:" << nRanks_ << " partitioning scheme: " << nPartX << 'x' << nPartY << std::endl;

    partPosX_ = rank_ % nPartX;
    partPosY_ = ((int)std::floor(rank_ / (nPartX))) % nPartY;
    partPosZ_ = std::floor(rank_ / (nPartX*nPartY));

    if(partPosX_ > 0)        // not the left-most partition
    {
        // so set the rank information and add an additional ghost layer
        leftRank_  = rank_ - 1;
    }
    if(partPosX_ < nPartX-1) // not the right-most partition
    {
        rightRank_  = rank_ + 1;
        uGhostLayer_ = 1;
    }
    if(partPosY_ > 0)
    {
        bottomRank_ = rank_ - nPartX;   // one line down
    }
    if(partPosY_ < nPartY-1)
    {
        topRank_ = rank_ + nPartX;
        vGhostLayer_ = 1;
    }
    if(partPosZ_ > 0)
    {
        hindRank_ = rank_ - nPartX*nPartY;
    }
    if(partPosZ_ < nPartZ-1)
    {
        frontRank_ = rank_ + nPartX*nPartY;
        wGhostLayer_ = 1;
    }
    
    // slice the whole domain into chunks
    double cellsPerPartitionX = (double)nCellsGlobal_[0] / (double)nPartX;
    double cellsPerPartitionY = (double)nCellsGlobal_[1] / (double)nPartY;
    double cellsPerPartitionZ = (double)nCellsGlobal_[2] / (double)nPartZ;
    nodeOffset_[0]  = std::ceil(cellsPerPartitionX * partPosX_);
    nodeOffset_[1]  = std::ceil(cellsPerPartitionY * partPosY_);
    nodeOffset_[2]  = std::ceil(cellsPerPartitionZ * partPosZ_);
    int xCells = std::ceil(cellsPerPartitionX);
    int yCells = std::ceil(cellsPerPartitionY);
    int zCells = std::ceil(cellsPerPartitionZ);

    // check, that two neighbouring partitions don't overlap (which happens in rare circumstances)
    int xPrevEnd = std::ceil(cellsPerPartitionX * (partPosX_ - 1)) + xCells;
    int yPrevEnd = std::ceil(cellsPerPartitionY * (partPosY_ - 1)) + yCells;
    int zPrevEnd = std::ceil(cellsPerPartitionZ * (partPosZ_ - 1)) + zCells;
    if(partPosX_ > 0 && xPrevEnd > nodeOffset_[0])
    {
        std::cout << "R:" << rank << " overlap on left boundary corrected\n";
        nodeOffset_[0]++;
        xCells--;
    }
    if(partPosY_ > 0 && yPrevEnd > nodeOffset_[1])
    {
        std::cout << "R:" << rank << " overlap on bottom boundary corrected\n";
        nodeOffset_[1]++;
        yCells--;
    }
    if(partPosZ_ > 0 && zPrevEnd > nodeOffset_[2])
    {
        std::cout << "R:" << rank << " overlap on hind boundary corrected\n";
        nodeOffset_[1]++;
        zCells--;
    }

    // it may be possible, that on the outer edge, the partition would be bigger, so use the smaller one
    nCellsLocal_[0] = std::min(nCellsGlobal_[0]-nodeOffset_[0], xCells);
    nCellsLocal_[1] = std::min(nCellsGlobal_[1]-nodeOffset_[1], yCells);
    nCellsLocal_[2] = std::min(nCellsGlobal_[2]-nodeOffset_[2], zCells);

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
	    std::cout << "Pn =  ARTITION Geometry " << nCellsLocal_[0] << 'x' << nCellsLocal_[1] << " A: "
                  << area << " E: " << averageEdges << " E/A: " << averageEdges/area << std::endl;
    }
    #endif
    

    // print debugging information about the partition
    std::cout << "\nR:" << rank_ << "\t@(" << partPosX_ << ',' << partPosY_ << ',' << partPosZ_ << ")"
        << "\tneighbours T:" << topRank_ << ", R:" << rightRank_
        << ", B:" << bottomRank_ << ", L:" << leftRank_ << ", H:" << hindRank_ << ", F:" << frontRank_ << "\tpartition domain: ["
        << nodeOffset_[0] << ", " << nodeOffset_[0]+nCellsLocal_[0]-1
        << "]x[" << nodeOffset_[1] << ", " << nodeOffset_[1]+nCellsLocal_[1]-1 << "]x[" << nodeOffset_[2] << ", " << nodeOffset_[2]+nCellsLocal_[2]-1 << "]\n";
}
