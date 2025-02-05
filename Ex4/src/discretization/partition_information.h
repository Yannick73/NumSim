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
    PartitionInformation(std::array<int, 3> nCellsGlobal, 
                         std::array<double, 3> meshWidth, int rank, int nRanks);
    //! get the local number of cells in the own subdomain
    inline std::array<int, 3> nCellsLocal()  const { return nCellsLocal_; }

    //! all cells within the partition (may not be used, so not optimized)
    inline std::size_t totalNoOfCellsLocal() const 
        { return nCellsLocal_[0] * nCellsLocal_[1] * nCellsLocal_[2]; }

    //! get the global number of cells in the whole computational domain
    //! used in OutputWriterParaviewParallel
    inline std::array<int,3> nCellsGlobal() const { return nCellsGlobal_; }

    //! all cells within the overall system
    inline int totalNoOfCellsGlobal()  const { return totalNoOfCellsGlobal_; }

    //! get the own MPI rank no
    //! used in OutputWriterParaviewParallel and OutputWriterTextParallel
    inline int ownRankNo() const { return rank_; }

    //! number of MPI ranks
    inline int nRanks() const { return nRanks_; }

    //! if the own partition has part of the bottom boundary of the whole domain
    inline bool ownBottomBoundary() const { return bottomRank_  == -1; }
    inline bool ownTopBoundary()    const { return topRank_     == -1; }
    inline bool ownLeftBoundary()   const { return leftRank_    == -1; }
    inline bool ownRightBoundary()  const { return rightRank_   == -1; }
    inline bool ownFrontBoundary()  const { return frontRank_   == -1; }
    inline bool ownHindBoundary()   const { return hindRank_    == -1; }

    //! get the rank no of the left neighbouring rank
    inline int bottomRank() const { return bottomRank_; }
    inline int topRank()    const { return topRank_; }
    inline int leftRank()   const { return leftRank_; }
    inline int rightRank()  const { return rightRank_; }
    inline int frontRank()  const { return frontRank_; }
    inline int hindRank()   const { return hindRank_; }

    //! get the offset values for counting local nodes in x and y direction. 
    //! (i_local,j_local) + nodeOffset = (i_global,j_global)
    //! used in OutputWriterParaviewParallel
    inline std::array<int, 3> nodeOffset() const { return nodeOffset_; }

    //! Get the meshWidth used
    inline std::array<double, 3> meshWidth() const { return meshWidth_; }

    //! Get the number of ghost layers
    inline int uGhostLayer() const { return uGhostLayer_; }
    inline int vGhostLayer() const { return vGhostLayer_; }
    inline int wGhostLayer() const { return wGhostLayer_; }
    //! Get the partition coordinates
    inline int getPartPosX() const { return partPosX_; }
    inline int getPartPosY() const { return partPosY_; }
    inline int getPartPosZ() const { return partPosZ_; }

private:
    //! rank information
    const int rank_;
    const int nRanks_;

    //! dimension information
    const std::array<int, 3> nCellsGlobal_;
    const std::array<double, 3> meshWidth_;
    std::array<int, 3> nCellsLocal_;
    std::array<int, 3> nodeOffset_;

    //! used in every residuum calculation, so optimize by pre-calculation
    const std::size_t totalNoOfCellsGlobal_;

    //! boundary information, -1 meaning, no neighbouring partition
    int bottomRank_ = -1;
    int topRank_    = -1;
    int leftRank_   = -1;
    int rightRank_  = -1;
    int frontRank_  = -1;
    int hindRank_   = -1;

    int uGhostLayer_ = 0;
    int vGhostLayer_ = 0;
    int wGhostLayer_ = 0;

    int partPosX_;
    int partPosY_;
    int partPosZ_;
};
