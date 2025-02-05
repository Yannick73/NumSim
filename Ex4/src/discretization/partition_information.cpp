#include "discretization/partition_information.h"

PartitionInformation::PartitionInformation(std::array<int, 3> nCellsGlobal,
                                           std::array<double, 3> meshWidth, 
                                           int rank, int nRanks) :
                                           rank_(rank), nRanks_(nRanks),
                                           nCellsGlobal_(nCellsGlobal),
                                           meshWidth_(meshWidth),
                                           totalNoOfCellsGlobal_(nCellsGlobal[0]*nCellsGlobal[1]*nCellsGlobal[2])
{
    #pragma message("Partinoning scheme still missing")

    assert(rank < nRanks);
    assert(rank >= 0);

    int nPartX = 1;
    int nPartY = 1;
    int nPartZ = 1;

    if(nRanks == 1)
    { }
    else if(nRanks == 2)
        nPartX = 2;
    else if(nRanks == 4)
    {
        nPartX = 2;
        nPartY = 2;
    }
    else if(nRanks == 6)
    {
        nPartX = 3;
        nPartY = 2;
    }
    else if(nRanks == 8)
    {
        nPartX = 2;
        nPartY = 2;
        nPartZ = 2;
    }
    else if(nRanks == 27)
    {
        nPartX = 3;
        nPartY = 3;
        nPartZ = 3;
    }
    else
    {
        std::stringstream message;
        message << "Number of ranks (" << nRanks << ") not yet supported!\n";
        throw std::runtime_error(message.str());
    }

    if(nRanks != nPartX*nPartY*nPartZ)
    {
        std::stringstream str;
        str << "The found partition scheme for " << nRanks << " was " 
            << nPartX << 'x' << nPartY << 'x' << nPartZ << " which does not multiply to nRanks.\n";
        throw std::runtime_error(str.str());
    }

    if(rank_ == 0)
        std::cout << "number of ranks:" << nRanks_ << " partitioning scheme: " << nPartX << 'x' << nPartY << std::endl;

    // ?
    partPosX_ = rank_ % nPartX;
    partPosY_ = ((int)std::floor(rank_ / (nPartX))) % nPartY;
    partPosZ_ = std::floor(rank_ / (nPartX*nPartY));

    // ???
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
        backRank_ = rank_ - nPartX*nPartY;
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
	    std::cout << "PARTITION Geometry " << nCellsLocal_[0] << 'x' << nCellsLocal_[1] << " A: "
                  << area << " E: " << averageEdges << " E/A: " << averageEdges/area << std::endl;
    }
    #endif
    

    // print debugging information about the partition
    std::cout << "\nR:" << rank_ << "\t@(" << partPosX_ << ',' << partPosY_ << ',' << partPosZ_ << ")"
        << "\tneighbours T:" << topRank_ << ", R:" << rightRank_
        << ", B:" << bottomRank_ << ", L:" << leftRank_ << ", H:" << backRank_ << ", F:" << frontRank_ << "\tpartition domain: ["
        << nodeOffset_[0] << ", " << nodeOffset_[0]+nCellsLocal_[0]-1
        << "]x[" << nodeOffset_[1] << ", " << nodeOffset_[1]+nCellsLocal_[1]-1 << "]x[" << nodeOffset_[2] << ", " << nodeOffset_[2]+nCellsLocal_[2]-1 << "]\n";
}
