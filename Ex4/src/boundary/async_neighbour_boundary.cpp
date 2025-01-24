#include "boundary/async_neighbour_boundary.h"

void AsyncNeighbourNorth::exchangeP()
{
    mpiHandler_.startReceive(recvBuf_, pBufLen_);
    #pragma omp simd
    for(int i = westOffset_; i < piLen_ - eastOffset_; i++)
    {
        sendBuf_[i] = p_(i,pjLen_-2);
    }
    mpiHandler_.send(sendBuf_, pBufLen_);
}

void AsyncNeighbourNorth::setRecvP()
{
    #pragma omp simd
    for(int i = westOffset_; i < piLen_ - eastOffset_; i++)
    {
        p_(i,pjLen_-1) = recvBuf_[i];
    }
}

void AsyncNeighbourNorth::exchangeRA()
{
    mpiHandler_.startReceive(recvBuf_, 2*pBufLen_);
    #pragma omp simd
    for(int i = westOffset_; i < piLen_ - eastOffset_; i++)
    {
        sendBuf_[2*i]   = (*r_)(i,pjLen_-2);
        sendBuf_[2*i+1] = (*a_)(i,pjLen_-2);
    }
    mpiHandler_.send(sendBuf_, 2*pBufLen_);
}

void AsyncNeighbourNorth::setRecvRA()
{
    #pragma omp simd
    for(int i = westOffset_; i < piLen_ - eastOffset_; i++)
    {
        (*r_)(i,pjLen_-1) = recvBuf_[2*i];
        (*a_)(i,pjLen_-1) = recvBuf_[2*i+1];
    }
}

void AsyncNeighbourNorth::exchangeUV()
{
    mpiHandler_.startReceive(recvBuf_, velBufLen_);
    #pragma omp simd
    for(int i = westOffset_; i < uiLen_ - eastOffset_; i++)
    {
        sendBuf_[i]   = u_(i,ujLen_-2);
    }
    #pragma omp simd
    for(int i = westOffset_; i < viLen_ - eastOffset_; i++)
    {
        sendBuf_[i + uiLen_]   = v_(i,vjLen_-2);
    }
    mpiHandler_.send(sendBuf_, velBufLen_);
}

void AsyncNeighbourNorth::setRecvUV()
{
    #pragma omp simd
    for(int i = westOffset_; i < uiLen_ - eastOffset_; i++)
    {
        u_(i,ujLen_-1) = recvBuf_[i];
    }
    #pragma omp simd
    for(int i = westOffset_; i < viLen_ - eastOffset_; i++)
    {
        v_(i,vjLen_-1) = recvBuf_[i + uiLen_];
    }
}

void AsyncNeighbourNorth::exchangeFG()
{
    mpiHandler_.startReceive(recvBuf_, velBufLen_);
    #pragma omp simd
    for(int i = westOffset_; i < uiLen_ - eastOffset_; i++)
    {
        sendBuf_[i] = f_(i,ujLen_-2);
    }
    #pragma omp simd
    for(int i = westOffset_; i < viLen_ - eastOffset_; i++)
    {
        sendBuf_[i + uiLen_]   = g_(i,vjLen_-2);
    }
    mpiHandler_.send(sendBuf_, velBufLen_);
}

void AsyncNeighbourNorth::setRecvFG()
{
    #pragma omp simd
    for(int i = westOffset_; i < uiLen_ - eastOffset_; i++)
    {
        f_(i,ujLen_-1) = recvBuf_[i];
    }
    #pragma omp simd
    for(int i = westOffset_; i < viLen_ - eastOffset_; i++)
    {
        g_(i,vjLen_-1) = recvBuf_[i + uiLen_];
    }
}

void AsyncNeighbourEast::exchangeP()
{
    mpiHandler_.startReceive(recvBuf_, pBufLen_);
    #pragma omp simd
    for(int j = 0; j < pjLen_; j++)
    {
        sendBuf_[j] = p_(piLen_-2,j);
    }
    MPI_Request sendReq;
    //! buffer, length, dataType, destRank, tag (variableType), MPI-context and callback (unused)
    mpiHandler_.send(sendBuf_, pBufLen_);
}

void AsyncNeighbourEast::setRecvP()
{
    for(int j = 0; j < pjLen_; j++)
    {
        p_(piLen_-1,j) = recvBuf_[j];
    }
}

void AsyncNeighbourEast::exchangeRA()
{
    mpiHandler_.startReceive(recvBuf_, 2*pBufLen_);
    #pragma omp simd
    for(int j = 0; j < pjLen_; j++)
    {
        sendBuf_[2*j]   = (*r_)(piLen_-2,j);
        sendBuf_[2*j+1] = (*a_)(piLen_-2,j);
    }
    mpiHandler_.send(sendBuf_, 2*pBufLen_);
}

void AsyncNeighbourEast::setRecvRA()
{
    #pragma omp simd
    for(int j = 0; j < pjLen_; j++)
    {
        (*r_)(piLen_-1,j) = recvBuf_[2*j];
        (*a_)(piLen_-1,j) = recvBuf_[2*j+1];
    }
}

void AsyncNeighbourEast::exchangeUV()
{
    mpiHandler_.startReceive(recvBuf_, velBufLen_);
    #pragma omp simd
    for(int j = 0; j < ujLen_; j++)
    {
        sendBuf_[j]   = u_(uiLen_-2,j);
    }
    #pragma omp simd
    for(int j = 0; j < vjLen_; j++)
    {
        sendBuf_[j + ujLen_]   = v_(viLen_-2,j);
    }
    mpiHandler_.send(sendBuf_, velBufLen_);
}

void AsyncNeighbourEast::setRecvUV()
{
    #pragma omp simd
    for(int j = 0; j < ujLen_; j++)
    {
        u_(uiLen_-1,j) = recvBuf_[j];
    }
    #pragma omp simd
    for(int j = 0; j < vjLen_; j++)
    {
        v_(viLen_-1,j) = recvBuf_[j + ujLen_];
    }
}

void AsyncNeighbourEast::exchangeFG()
{
    mpiHandler_.startReceive(recvBuf_, velBufLen_);
    #pragma omp simd
    for(int j = 0; j < ujLen_; j++)
    {
        sendBuf_[j]   = f_(uiLen_-2,j);
    }
    #pragma omp simd
    for(int j = 0; j < vjLen_; j++)
    {
        sendBuf_[j + ujLen_]   = g_(viLen_-2,j);
    }
    mpiHandler_.send(sendBuf_, velBufLen_);
}

void AsyncNeighbourEast::setRecvFG()
{
    #pragma omp simd
    for(int j = 0; j < ujLen_; j++)
    {
        f_(uiLen_-1,j) = recvBuf_[j];
    }
    #pragma omp simd
    for(int j = 0; j < vjLen_; j++)
    {
        g_(viLen_-1,j) = recvBuf_[j + ujLen_];
    }
}

void AsyncNeighbourSouth::exchangeP()
{
    mpiHandler_.startReceive(recvBuf_, pBufLen_);
    #pragma omp simd
    for(int i = westOffset_; i < piLen_ - eastOffset_; i++)
    {
        sendBuf_[i] = p_(i,1);
    }
    mpiHandler_.send(sendBuf_, pBufLen_);
}

void AsyncNeighbourSouth::setRecvP()
{
    #pragma omp simd
    for(int i = westOffset_; i < piLen_ - eastOffset_; i++)
    {
        p_(i,0) = recvBuf_[i];
    }
}

void AsyncNeighbourSouth::exchangeRA()
{
    mpiHandler_.startReceive(recvBuf_, 2*pBufLen_);
    #pragma omp simd
    for(int i = westOffset_; i < piLen_ - eastOffset_; i++)
    {
        sendBuf_[2*i]   = (*r_)(i,1);
        sendBuf_[2*i+1] = (*a_)(i,1);
    }
    mpiHandler_.send(sendBuf_, 2*pBufLen_);
}

void AsyncNeighbourSouth::setRecvRA()
{
    #pragma omp simd
    for(int i = westOffset_; i < piLen_ - eastOffset_; i++)
    {
        (*r_)(i,0) = recvBuf_[2*i];
        (*a_)(i,0) = recvBuf_[2*i+1];
    }
}

void AsyncNeighbourSouth::exchangeUV()
{
    mpiHandler_.startReceive(recvBuf_, velBufLen_);
    #pragma omp simd
    for(int i = westOffset_; i < uiLen_ - eastOffset_; i++)
    {
        sendBuf_[i] = u_(i,1);
    }
    #pragma omp simd
    for(int i = westOffset_; i < viLen_ - eastOffset_; i++)
    {
        sendBuf_[i + uiLen_]   = v_(i,1);
    }
    mpiHandler_.send(sendBuf_, velBufLen_);
}

void AsyncNeighbourSouth::setRecvUV()
{
    #pragma omp simd
    for(int i = westOffset_; i < uiLen_ - eastOffset_; i++)
    {
        u_(i,0) = recvBuf_[i];
    }
    #pragma omp simd
    for(int i = westOffset_; i < viLen_ - eastOffset_; i++)
    {
        v_(i,0) = recvBuf_[i + uiLen_];
    }
}

void AsyncNeighbourSouth::exchangeFG()
{
    mpiHandler_.startReceive(recvBuf_, velBufLen_);
    #pragma omp simd
    for(int i = westOffset_; i < uiLen_ - eastOffset_; i++)
    {
        sendBuf_[i]   = f_(i,1);
    }
    #pragma omp simd
    for(int i = westOffset_; i < viLen_ - eastOffset_; i++)
    {
        sendBuf_[i + uiLen_]   = g_(i,1);
    }
    mpiHandler_.send(sendBuf_, velBufLen_);
}

void AsyncNeighbourSouth::setRecvFG()
{
    #pragma omp simd
    for(int i = westOffset_; i < uiLen_ - eastOffset_; i++)
    {
        f_(i,0) = recvBuf_[i];
    }
    #pragma omp simd
    for(int i = westOffset_; i < viLen_ - eastOffset_; i++)
    {
        g_(i,0) = recvBuf_[i + uiLen_];
    }
}

void AsyncNeighbourWest::exchangeP()
{
    mpiHandler_.startReceive(recvBuf_, pBufLen_);
    #pragma omp simd
    for(int j = 0; j < pjLen_; j++)
    {
        sendBuf_[j] = p_(1,j);
    }
    mpiHandler_.send(sendBuf_, pBufLen_);
}

void AsyncNeighbourWest::setRecvP()
{
    #pragma omp simd
    for(int j = 0; j < pjLen_; j++)
    {
        p_(0,j) = recvBuf_[j];
    }
}

void AsyncNeighbourWest::exchangeRA()
{
    mpiHandler_.startReceive(recvBuf_, 2*pBufLen_);
    #pragma omp simd
    for(int j = 0; j < pjLen_; j++)
    {
        sendBuf_[2*j]   = (*r_)(1,j);
        sendBuf_[2*j+1] = (*a_)(1,j);
    }
    mpiHandler_.send(sendBuf_, 2*pBufLen_);
}

void AsyncNeighbourWest::setRecvRA()
{
    #pragma omp simd
    for(int j = 0; j < pjLen_; j++)
    {
        (*r_)(0,j) = recvBuf_[2*j];
        (*a_)(0,j) = recvBuf_[2*j+1];
    }
}

void AsyncNeighbourWest::exchangeUV()
{
    mpiHandler_.startReceive(recvBuf_, velBufLen_);
    #pragma omp simd
    for(int j = 0; j < ujLen_; j++)
    {
        sendBuf_[j]   = u_(1,j);
    }
    #pragma omp simd
    for(int j = 0; j < vjLen_; j++)
    {
        sendBuf_[j + ujLen_]   = v_(1,j);
    }
    mpiHandler_.send(sendBuf_, velBufLen_);
}

void AsyncNeighbourWest::setRecvUV()
{
    #pragma omp simd
    for(int j = 0; j < ujLen_; j++)
    {
        u_(0,j) = recvBuf_[j];
    }
    #pragma omp simd
    for(int j = 0; j < vjLen_; j++)
    {
        v_(0,j) = recvBuf_[j + ujLen_];
    }
}

void AsyncNeighbourWest::exchangeFG()
{
    mpiHandler_.startReceive(recvBuf_, velBufLen_);
    #pragma omp simd
    for(int j = 0; j < ujLen_; j++)
    {
        sendBuf_[j]   = f_(1,j);
    }
    #pragma omp simd
    for(int j = 0; j < vjLen_; j++)
    {
        sendBuf_[j + ujLen_]   = g_(1,j);
    }
    mpiHandler_.send(sendBuf_, velBufLen_);
}

void AsyncNeighbourWest::setRecvFG()
{
    #pragma omp simd
    for(int j = 0; j < ujLen_; j++)
    {
        f_(0,j) = recvBuf_[j];
    }
    #pragma omp simd
    for(int j = 0; j < vjLen_; j++)
    {
        g_(0,j) = recvBuf_[j + ujLen_];
    }
}