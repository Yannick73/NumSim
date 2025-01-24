#include "discretization/partition_information.h"

PartitionInformation::PartitionInformation(std::array<int, 2> nCellsGlobal,
                                           std::array<double, 2> meshWidth, 
                                           int rank, int nRanks) :
                                           rank_(rank), nRanks_(nRanks),
                                           nCellsGlobal_(nCellsGlobal),
                                           meshWidth_(meshWidth),
                                           totalNoOfCellsGlobal_(nCellsGlobal[0]*nCellsGlobal[1])
{
    assert(rank < nRanks);
    assert(rank >= 0);

    int nPartX = 1;
    int nPartY = 1;

    // optimize to find the ratio relative to the smaller enumerator and bigger denominator
    double aspectRatio = ((double)std::min(nCellsGlobal[0], nCellsGlobal[1])) 
                         / ((double)std::max(nCellsGlobal[0], nCellsGlobal[1]));
    int enumerator  = 1;
    int denominator = nRanks;
    double bestFit = std::abs(aspectRatio - ((double)enumerator)/((double)denominator));
    // linear optimization problem
    for(int i = 2; i < nRanks; i++)
    {
        // if i is a divider of nRanks, then it may be a candidate for one of the edges
        if(nRanks % i == 0)
        {
            int divisor = nRanks / i;
            double fit = std::abs(aspectRatio - ((double)i)/((double)divisor));
            if(fit < bestFit)
            {
                bestFit = fit;
                enumerator = i;
                denominator = divisor;
            }
        }
    }
    if(nCellsGlobal[0] > nCellsGlobal[1])
    {
        // while the enumerator is likely to be the smaller one, 
        // just make it explicit, to avoid any issues
        nPartX = std::max(enumerator, denominator);
        nPartY = std::min(enumerator, denominator);
    }
    else
    {
        nPartX = std::min(enumerator, denominator);
        nPartY = std::max(enumerator, denominator);
    }
    
    // make this explicit (asserts are ignored in release-mode),
    // if this does not hold, the program does not work
    if(nRanks != nPartX*nPartY)
    {
        std::stringstream str;
        str << "The found partition scheme for " << nRanks << " was " 
            << nPartX << 'x' << nPartY << " which does not multiply to nRanks.\n";
        throw std::runtime_error(str.str());
    }

    if(rank_ == 0)
        std::cout << "number of ranks:" << nRanks_ << " partitioning scheme: " << nPartX << 'x' << nPartY << std::endl;
    
    // counted left to right and then bottom to top (like in the grid)
    partPosX_ = rank_ % nPartX;
    partPosY_ = std::floor(rank_ / nPartX);

    if(partPosX_ > 0)        // not the left-most partition
    {
        // so set the rank information and add an additional ghost layer
        westBoundaryRank_  = rank_ - 1;
    }
    if(partPosX_ < nPartX-1) // not the right-most partition
    {
        eastBoundaryRank_  = rank_ + 1;
        uGhostLayer_ = 1;
    }
    if(partPosY_ > 0)
    {
        southBoundaryRank_ = rank_ - nPartX;   // one line down
    }
    if(partPosY_ < nPartY-1)
    {
        northBoundaryRank_ = rank_ + nPartX;
        vGhostLayer_ = 1;
    }

    // slice the whole domain into chunks
    double cellsPerPartitionX = (double)nCellsGlobal_[0] / (double)nPartX;
    double cellsPerPartitionY = (double)nCellsGlobal_[1] / (double)nPartY;
    nodeOffset_[0]  = std::ceil(cellsPerPartitionX * partPosX_);
    nodeOffset_[1]  = std::ceil(cellsPerPartitionY * partPosY_);
    int xCells = std::ceil(cellsPerPartitionX);
    int yCells = std::ceil(cellsPerPartitionY);
    
    // check, that two neighbouring partitions don't overlap (which happens in rare circumstances)
    int xPrevEnd = std::ceil(cellsPerPartitionX * (partPosX_ - 1)) + xCells;
    int yPrevEnd = std::ceil(cellsPerPartitionY * (partPosY_ - 1)) + yCells;
    if(partPosX_ > 0 && xPrevEnd > nodeOffset_[0])
    {
        std::cout << "R:" << rank << " overlap on west boundary corrected\n";
        nodeOffset_[0]++;
        xCells--;
    }
    if(partPosY_ > 0 && yPrevEnd > nodeOffset_[1])
    {
        std::cout << "R:" << rank << " overlap on south boundary corrected\n";
        nodeOffset_[1]++;
        yCells--;
    }

    // it may be possible, that on the outer edge, the partition would be bigger, so use the smaller one
    nCellsLocal_[0] = std::min(nCellsGlobal_[0]-nodeOffset_[0], xCells);
    nCellsLocal_[1] = std::min(nCellsGlobal_[1]-nodeOffset_[1], yCells);

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
    std::cout << "R:" << rank_ << "\t@(" << partPosX_ << ',' << partPosY_ << ")"
        << "\tneighbours N:" << northBoundaryRank_ << ", E:" << eastBoundaryRank_
        << ", S:" << southBoundaryRank_ << ", W:" << westBoundaryRank_ << "\tpartition domain: ["
        << nodeOffset_[0] << ", " << nodeOffset_[0]+nCellsLocal_[0]-1
        << "]x[" << nodeOffset_[1] << ", " << nodeOffset_[1]+nCellsLocal_[1]-1 << "]\n";
}
