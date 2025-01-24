#include "boundary/neighbour_boundary.h"

/*
    Scenario 3_5 (8x4 with 160x120 cells) 20 cores on simcl1
    Edge prio means order of edge depending on the partitionX position
    Exchange prio means ordering receive->send and send->receive for corresponding corners

    Edge prio & exchange prio: 18.8s, 18.8s, 18.72
    no Edge prio & exchange prio: 18.28s, 18.24s, 18.12s
    Edge prio & no exchange prio: 16.18s, 16.08s, 16.17s
    no Edge prio & no exchange prio: 17.52s, 17.45s, 17.5

    Previously foundout in an informal test, that using immediate receive->send with setting the first 
    incoming data first took longer. But should be investigated proberly,
    as this model would completely deadlock immune (as opposed to the current one).
*/

void NeighbourNorth::setP()
{
    // generate buffers for MPI exchange, scoped output buffer is adequate,
    // as the send is finished at the time, when the receive is done, given the corresponding rank does the same
    std::vector<double> oBuf(pBufLen_, 0.0);
    std::vector<double> iBuf(pBufLen_, 0.0);
    #pragma omp simd
    // set the content of the output with our own cells
    for(int i = westOffset_; i < piLen_ - eastOffset_; i++)
    {
        oBuf[i] = p_(i,pjLen_-2);
    }
    // send request of non-blocking send is ignored, but could be used with wait, if in doubt
    MPI_Request sendReq;
    MPI_Isend(oBuf.data(), pBufLen_, MPI_DOUBLE, neighbourRank_, data::P, MPI_COMM_WORLD, &sendReq);
    MPI_Recv (iBuf.data(), pBufLen_, MPI_DOUBLE, neighbourRank_, data::P, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    // as soon as the data from the neighbour is received, it is written on our own boundary
    #pragma omp simd
    for(int i = westOffset_; i < piLen_ - eastOffset_; i++)
    {
        p_(i,pjLen_-1) = iBuf[i];
    }
}

// these functions were developed for CG, which is still under developement
void NeighbourNorth::setRA()
{
    std::vector<double> oBuf(2*pBufLen_, 0.0);
    std::vector<double> iBuf(2*pBufLen_, 0.0);
    #pragma omp simd
    for(int i = westOffset_; i < piLen_ - eastOffset_; i++)
    {
        oBuf[2*i]   = (*r_)(i,pjLen_-2);
        oBuf[2*i+1] = (*a_)(i,pjLen_-2);
    }

    MPI_Request sendReq;
    MPI_Isend(oBuf.data(), 2*pBufLen_, MPI_DOUBLE, neighbourRank_, data::RA, MPI_COMM_WORLD, &sendReq);
    MPI_Recv (iBuf.data(), 2*pBufLen_, MPI_DOUBLE, neighbourRank_, data::RA, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    #pragma omp simd
    for(int i = westOffset_; i < piLen_ - eastOffset_; i++)
    {
        (*r_)(i,pjLen_-1) = iBuf[2*i];
        (*a_)(i,pjLen_-1) = iBuf[2*i+1];
    }
}

void NeighbourNorth::setUV()
{
    std::vector<double> oBuf(velBufLen_, 0.0);
    std::vector<double> iBuf(velBufLen_, 0.0);
    // first fill our u data in buffer
    #pragma omp simd
    for(int i = westOffset_; i < uiLen_ - eastOffset_; i++)
    {
        oBuf[i]   = u_(i,ujLen_-2);
    }
    // and then with an offset fill in v data
    #pragma omp simd
    for(int i = westOffset_; i < viLen_ - eastOffset_; i++)
    {
        oBuf[i + uiLen_]   = v_(i,vjLen_-2);
    }
    // before sending the combined buffer (to save on communication)
    MPI_Request sendReq;
    MPI_Isend(oBuf.data(), velBufLen_, MPI_DOUBLE, neighbourRank_, data::UV, MPI_COMM_WORLD, &sendReq);
    MPI_Recv (iBuf.data(), velBufLen_, MPI_DOUBLE, neighbourRank_, data::UV, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    // same goes for receiveing
    #pragma omp simd
    for(int i = westOffset_; i < uiLen_ - eastOffset_; i++)
    {
        u_(i,ujLen_-1) = iBuf[i];
    }
    #pragma omp simd
    for(int i = westOffset_; i < viLen_ - eastOffset_; i++)
    {
        v_(i,vjLen_-1) = iBuf[i + uiLen_];
    }
}

// same rules as with UV apply
void NeighbourNorth::setFG()
{
    std::vector<double> oBuf(velBufLen_, 0.0);
    std::vector<double> iBuf(velBufLen_, 0.0);
    #pragma omp simd
    for(int i = westOffset_; i < uiLen_ - eastOffset_; i++)
    {
        oBuf[i] = f_(i,ujLen_-2);
    }
    #pragma omp simd
    for(int i = westOffset_; i < viLen_ - eastOffset_; i++)
    {
        oBuf[i + uiLen_]   = g_(i,vjLen_-2);
    }
    MPI_Request sendReq;
    MPI_Isend(oBuf.data(), velBufLen_, MPI_DOUBLE, neighbourRank_, data::FG, MPI_COMM_WORLD, &sendReq);
    MPI_Recv (iBuf.data(), velBufLen_, MPI_DOUBLE, neighbourRank_, data::FG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    #pragma omp simd
    for(int i = westOffset_; i < uiLen_ - eastOffset_; i++)
    {
        f_(i,ujLen_-1) = iBuf[i];
    }
    #pragma omp simd
    for(int i = westOffset_; i < viLen_ - eastOffset_; i++)
    {
        g_(i,vjLen_-1) = iBuf[i + uiLen_];
    }
}

// works the same way as NeighbourNorth, but using other indices and always setting the corners
void NeighbourEast::setP()
{
    std::vector<double> oBuf(pBufLen_, 0.0);
    std::vector<double> iBuf(pBufLen_, 0.0);
    #pragma omp simd
    for(int j = 0; j < pjLen_; j++)
    {
        oBuf[j] = p_(piLen_-2,j);
    }
    MPI_Request sendReq;
    MPI_Isend(oBuf.data(), pBufLen_, MPI_DOUBLE, neighbourRank_, data::P, MPI_COMM_WORLD, &sendReq);
    MPI_Recv (iBuf.data(), pBufLen_, MPI_DOUBLE, neighbourRank_, data::P, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    #pragma omp simd
    for(int j = 0; j < pjLen_; j++)
    {
        p_(piLen_-1,j) = iBuf[j];
    }
}

void NeighbourEast::setRA()
{
    std::vector<double> oBuf(2*pBufLen_, 0.0);
    std::vector<double> iBuf(2*pBufLen_, 0.0);
    #pragma omp simd
    for(int j = 0; j < pjLen_; j++)
    {
        oBuf[2*j]   = (*r_)(piLen_-2,j);
        oBuf[2*j+1] = (*a_)(piLen_-2,j);
    }

    MPI_Request sendReq;
    MPI_Isend(oBuf.data(), 2*pBufLen_, MPI_DOUBLE, neighbourRank_, data::RA, MPI_COMM_WORLD, &sendReq);
    MPI_Recv (iBuf.data(), 2*pBufLen_, MPI_DOUBLE, neighbourRank_, data::RA, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    #pragma omp simd
    for(int j = 0; j < pjLen_; j++)
    {
        (*r_)(piLen_-1,j) = iBuf[2*j];
        (*a_)(piLen_-1,j) = iBuf[2*j+1];
    }
}

void NeighbourEast::setUV()
{
    std::vector<double> oBuf(velBufLen_, 0.0);
    std::vector<double> iBuf(velBufLen_, 0.0);
    #pragma omp simd
    for(int j = 0; j < ujLen_; j++)
    {
        oBuf[j]   = u_(uiLen_-2,j);
    }
    #pragma omp simd
    for(int j = 0; j < vjLen_; j++)
    {
        oBuf[j + ujLen_]   = v_(viLen_-2,j);
    }
    MPI_Request sendReq;
    MPI_Isend(oBuf.data(), velBufLen_, MPI_DOUBLE, neighbourRank_, data::UV, MPI_COMM_WORLD, &sendReq);
    MPI_Recv (iBuf.data(), velBufLen_, MPI_DOUBLE, neighbourRank_, data::UV, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    
    #pragma omp simd
    for(int j = 0; j < ujLen_; j++)
    {
        u_(uiLen_-1,j) = iBuf[j];
    }
    #pragma omp simd
    for(int j = 0; j < vjLen_; j++)
    {
        v_(viLen_-1,j) = iBuf[j + ujLen_];
    }
}

void NeighbourEast::setFG()
{
    std::vector<double> oBuf(velBufLen_, 0.0);
    std::vector<double> iBuf(velBufLen_, 0.0);
    #pragma omp simd
    for(int j = 0; j < ujLen_; j++)
    {
        oBuf[j]   = f_(uiLen_-2,j);
    }
    #pragma omp simd
    for(int j = 0; j < vjLen_; j++)
    {
        oBuf[j + ujLen_]   = g_(viLen_-2,j);
    }
    MPI_Request sendReq;
    MPI_Isend(oBuf.data(), velBufLen_, MPI_DOUBLE, neighbourRank_, data::FG, MPI_COMM_WORLD, &sendReq);
    MPI_Recv (iBuf.data(), velBufLen_, MPI_DOUBLE, neighbourRank_, data::FG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    
    #pragma omp simd
    for(int j = 0; j < ujLen_; j++)
    {
        f_(uiLen_-1,j) = iBuf[j];
    }
    #pragma omp simd
    for(int j = 0; j < vjLen_; j++)
    {
        g_(viLen_-1,j) = iBuf[j + ujLen_];
    }
}

// works the same way as NeighbourNorth, but using other indices
void NeighbourSouth::setP()
{
    std::vector<double> oBuf(pBufLen_, 0.0);
    std::vector<double> iBuf(pBufLen_, 0.0);
    #pragma omp simd
    for(int i = westOffset_; i < piLen_ - eastOffset_; i++)
    {
        oBuf[i] = p_(i,1);
    }
    MPI_Request sendReq;
    MPI_Isend(oBuf.data(), pBufLen_, MPI_DOUBLE, neighbourRank_, data::P, MPI_COMM_WORLD, &sendReq);
    MPI_Recv (iBuf.data(), pBufLen_, MPI_DOUBLE, neighbourRank_, data::P, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    #pragma omp simd
    for(int i = westOffset_; i < piLen_ - eastOffset_; i++)
    {
        p_(i,0) = iBuf[i];
    }
}



void NeighbourSouth::setRA()
{
    std::vector<double> oBuf(2*pBufLen_, 0.0);
    std::vector<double> iBuf(2*pBufLen_, 0.0);
    #pragma omp simd
    for(int i = westOffset_; i < piLen_ - eastOffset_; i++)
    {
        oBuf[2*i]   = (*r_)(i,1);
        oBuf[2*i+1] = (*a_)(i,1);
    }
    
    MPI_Request sendReq;
    MPI_Isend(oBuf.data(), 2*pBufLen_, MPI_DOUBLE, neighbourRank_, data::RA, MPI_COMM_WORLD, &sendReq);
    MPI_Recv (iBuf.data(), 2*pBufLen_, MPI_DOUBLE, neighbourRank_, data::RA, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    #pragma omp simd
    for(int i = westOffset_; i < piLen_ - eastOffset_; i++)
    {
        (*r_)(i,0) = iBuf[2*i];
        (*a_)(i,0) = iBuf[2*i+1];
    }
}

void NeighbourSouth::setUV()
{
    std::vector<double> oBuf(velBufLen_, 0.0);
    std::vector<double> iBuf(velBufLen_, 0.0);
    #pragma omp simd
    for(int i = westOffset_; i < uiLen_ - eastOffset_; i++)
    {
        oBuf[i] = u_(i,1);
    }
    #pragma omp simd
    for(int i = westOffset_; i < viLen_ - eastOffset_; i++)
    {
        oBuf[i + uiLen_]   = v_(i,1);
    }
    MPI_Request sendReq;
    MPI_Isend(oBuf.data(), velBufLen_, MPI_DOUBLE, neighbourRank_, data::UV, MPI_COMM_WORLD, &sendReq);
    MPI_Recv (iBuf.data(), velBufLen_, MPI_DOUBLE, neighbourRank_, data::UV, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    #pragma omp simd
    for(int i = westOffset_; i < uiLen_ - eastOffset_; i++)
    {
        u_(i,0) = iBuf[i];
    }
    #pragma omp simd
    for(int i = westOffset_; i < viLen_ - eastOffset_; i++)
    {
        v_(i,0) = iBuf[i + uiLen_];
    }
}

void NeighbourSouth::setFG()
{
    std::vector<double> oBuf(velBufLen_, 0.0);
    std::vector<double> iBuf(velBufLen_, 0.0);
    #pragma omp simd
    for(int i = westOffset_; i < uiLen_ - eastOffset_; i++)
    {
        oBuf[i]   = f_(i,1);
    }
    #pragma omp simd
    for(int i = westOffset_; i < viLen_ - eastOffset_; i++)
    {
        oBuf[i + uiLen_]   = g_(i,1);
    }
    MPI_Request sendReq;
    MPI_Isend(oBuf.data(), velBufLen_, MPI_DOUBLE, neighbourRank_, data::FG, MPI_COMM_WORLD, &sendReq);
    MPI_Recv (iBuf.data(), velBufLen_, MPI_DOUBLE, neighbourRank_, data::FG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    #pragma omp simd
    for(int i = westOffset_; i < uiLen_ - eastOffset_; i++)
    {
        f_(i,0) = iBuf[i];
    }
    #pragma omp simd
    for(int i = westOffset_; i < viLen_ - eastOffset_; i++)
    {
        g_(i,0) = iBuf[i + uiLen_];
    }
}

// works the same way as NeighbourNorth, but using other indices and always setting the corners
void NeighbourWest::setP()
{
    std::vector<double> oBuf(pBufLen_, 0.0);
    std::vector<double> iBuf(pBufLen_, 0.0);
    #pragma omp simd
    for(int j = 0; j < pjLen_; j++)
    {
        oBuf[j] = p_(1,j);
    }
    MPI_Request sendReq;
    MPI_Isend(oBuf.data(), pBufLen_, MPI_DOUBLE, neighbourRank_, data::P, MPI_COMM_WORLD, &sendReq);
    MPI_Recv (iBuf.data(), pBufLen_, MPI_DOUBLE, neighbourRank_, data::P, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    #pragma omp simd
    for(int j = 0; j < pjLen_; j++)
    {
        p_(0,j) = iBuf[j];
    }
}

void NeighbourWest::setRA()
{
    std::vector<double> oBuf(2*pBufLen_, 0.0);
    std::vector<double> iBuf(2*pBufLen_, 0.0);
    #pragma omp simd
    for(int j = 0; j < pjLen_; j++)
    {
        oBuf[2*j]   = (*r_)(1,j);
        oBuf[2*j+1] = (*a_)(1,j);
    }
    MPI_Request sendReq;
    MPI_Isend(oBuf.data(), 2*pBufLen_, MPI_DOUBLE, neighbourRank_, data::RA, MPI_COMM_WORLD, &sendReq);
    MPI_Recv (iBuf.data(), 2*pBufLen_, MPI_DOUBLE, neighbourRank_, data::RA, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    #pragma omp simd
    for(int j = 0; j < pjLen_; j++)
    {
        (*r_)(0,j) = iBuf[2*j];
        (*a_)(0,j) = iBuf[2*j+1];
    }
}

void NeighbourWest::setUV()
{
    std::vector<double> oBuf(velBufLen_, 0.0);
    std::vector<double> iBuf(velBufLen_, 0.0);
    #pragma omp simd
    for(int j = 0; j < ujLen_; j++)
    {
        oBuf[j]   = u_(1,j);
    }
    #pragma omp simd
    for(int j = 0; j < vjLen_; j++)
    {
        oBuf[j + ujLen_]   = v_(1,j);
    }
    MPI_Request sendReq;
    MPI_Isend(oBuf.data(), velBufLen_, MPI_DOUBLE, neighbourRank_, data::UV, MPI_COMM_WORLD, &sendReq);
    MPI_Recv (iBuf.data(), velBufLen_, MPI_DOUBLE, neighbourRank_, data::UV, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    #pragma omp simd
    for(int j = 0; j < ujLen_; j++)
    {
        u_(0,j) = iBuf[j];
    }
    #pragma omp simd
    for(int j = 0; j < vjLen_; j++)
    {
        v_(0,j) = iBuf[j + ujLen_];
    }
}

void NeighbourWest::setFG()
{
    std::vector<double> oBuf(velBufLen_, 0.0);
    std::vector<double> iBuf(velBufLen_, 0.0);
    #pragma omp simd
    for(int j = 0; j < ujLen_; j++)
    {
        oBuf[j]   = f_(1,j);
    }
    #pragma omp simd
    for(int j = 0; j < vjLen_; j++)
    {
        oBuf[j + ujLen_]   = g_(1,j);
    }
    MPI_Request sendReq;
    MPI_Isend(oBuf.data(), velBufLen_, MPI_DOUBLE, neighbourRank_, data::FG, MPI_COMM_WORLD, &sendReq);
    MPI_Recv (iBuf.data(), velBufLen_, MPI_DOUBLE, neighbourRank_, data::FG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    #pragma omp simd
    for(int j = 0; j < ujLen_; j++)
    {
        f_(0,j) = iBuf[j];
    }
    #pragma omp simd
    for(int j = 0; j < vjLen_; j++)
    {
        g_(0,j) = iBuf[j + ujLen_];
    }
}
