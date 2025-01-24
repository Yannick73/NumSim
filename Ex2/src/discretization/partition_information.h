#pragma once

#include <array>
#include <cassert>
#include <exception>
#include <sstream>
#include <cmath>
#include <iostream>
#include <mpi.h>

//! this primarily serves the purpose of a data object calculating all relevant data
class PartitionInformation
{
public:
    PartitionInformation(std::array<int, 2> nCellsGlobal, 
                         std::array<double, 2> meshWidth, int rank, int nRanks);
    //! get the local number of cells in the own subdomain
    inline std::array<int,2> nCellsLocal()  const { return nCellsLocal_; }

    //! all cells within the partition (may not be used, so not optimized)
    inline int totalNoOfCellsLocal() const { return nCellsLocal_[0] * nCellsLocal_[1]; }

    //! get the global number of cells in the whole computational domain
    //! used in OutputWriterParaviewParallel
    inline std::array<int,2> nCellsGlobal() const { return nCellsGlobal_; }

    //! all cells within the overall system
    inline int totalNoOfCellsGlobal()  const { return totalNoOfCellsGlobal_; }

    //! get the own MPI rank no
    //! used in OutputWriterParaviewParallel and OutputWriterTextParallel
    inline int ownRankNo() const { return rank_; }

    //! number of MPI ranks
    inline int nRanks() const { return nRanks_; }

    //! if the own partition has part of the bottom boundary of the whole domain
    inline bool ownPartitionContainsBottomBoundary() const { return southBoundaryRank_ == -1; }

    //! if the own partition has part of the top boundary of the whole domain
    //! used in OutputWriterParaviewParallel
    inline bool ownPartitionContainsTopBoundary()    const { return northBoundaryRank_ == -1; }

    //! if the own partition has part of the left boundary of the whole domain
    inline bool ownPartitionContainsLeftBoundary()   const { return westBoundaryRank_ == -1; }

    //! if the own partition has part of the right boundary of the whole domain
    //! used in OutputWriterParaviewParallel
    inline bool ownPartitionContainsRightBoundary()  const { return eastBoundaryRank_ == -1; }

    //! get the rank no of the left neighbouring rank
    inline int leftNeighbourRankNo()   const { return westBoundaryRank_; }

    //! get the rank no of the right neighbouring rank
    inline int rightNeighbourRankNo()  const { return eastBoundaryRank_; }

    //! get the rank no of the top neighbouring rank
    inline int topNeighbourRankNo()    const { return northBoundaryRank_; }

    //! get the rank no of the bottom neighbouring rank
    inline int bottomNeighbourRankNo() const { return southBoundaryRank_; }

    //! get the offset values for counting local nodes in x and y direction. 
    //! (i_local,j_local) + nodeOffset = (i_global,j_global)
    //! used in OutputWriterParaviewParallel
    inline std::array<int,2> nodeOffset() const { return nodeOffset_; }

    //! Get the meshWidth used
    inline std::array<double, 2> meshWidth() const { return meshWidth_; }

    //! Get the number of ghost layers
    inline int uGhostLayer() const { return uGhostLayer_; }
    inline int vGhostLayer() const { return vGhostLayer_; }

    //! Get the partition coordinates
    inline int getPartPosX() const { return partPosX_; }
    inline int getPartPosY() const { return partPosY_; }

private:
    //! rank information
    const int rank_;
    const int nRanks_;

    //! dimension information
    const std::array<int, 2> nCellsGlobal_;
    const std::array<double, 2> meshWidth_;
    std::array<int, 2> nCellsLocal_;
    std::array<int, 2> nodeOffset_;

    //! used in every residuum calculation, so optimize by pre-calculation
    const int totalNoOfCellsGlobal_;

    //! boundary information, -1 meaning, no neighbouring partition
    int northBoundaryRank_ = -1;
    int eastBoundaryRank_  = -1;
    int southBoundaryRank_ = -1;
    int westBoundaryRank_  = -1;

    int uGhostLayer_ = 0;
    int vGhostLayer_ = 0;

    int partPosX_;
    int partPosY_;
};
