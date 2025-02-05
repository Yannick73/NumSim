#include "boundary/async_neighbour_boundary.h"

void AsyncNeighbourTop::exchangeP()
{
    mpiHandler_.startReceive(pRecvBuf_.data(), pRecvBuf_.length());
    #pragma omp simd collapse(2)
    for(int k = hindOffset_; k < pkLen_ - frontOffset_; k++)
    {
        for(int i = leftOffset_; i < piLen_ - rightOffset_; i++)
        {
            pSendBuf_(i,k) = p_(i,pjLen_-2,k);
        }
    }
    mpiHandler_.send(pSendBuf_.data(), pSendBuf_.length());
}

void AsyncNeighbourTop::setRecvP()
{
    #pragma omp simd collapse(2)
    for(int k = hindOffset_; k < pkLen_ - frontOffset_; k++)
    {
        for(int i = leftOffset_; i < piLen_ - rightOffset_; i++)
        {
            p_(i,pjLen_-1,k) = pRecvBuf_(i,k);
        }
    }
}

void AsyncNeighbourTop::exchangeUVW()
{
    mpiHandler_.startReceive(velRecvBuf_.data(), velRecvBuf_.length());
    #pragma omp simd collapse(2)
    for(int k = hindOffset_; k < ukLen_ - frontOffset_; k++)
    {
        for(int i = leftOffset_; i < uiLen_ - rightOffset_; i++)
        {
            velSendBuf_(i,k,ind::X) = u_(i,ujLen_-2,k);
        }
    }
    #pragma omp simd collapse(2)
    for(int k = hindOffset_; k < vkLen_ - frontOffset_; k++)
    {
        for(int i = leftOffset_; i < viLen_ - rightOffset_; i++)
        {
            velSendBuf_(i,k,ind::Y) = v_(i,vjLen_-2,k);
        }
    }
    #pragma omp simd collapse(2)
    for(int k = hindOffset_; k < wkLen_ - frontOffset_; k++)
    {
        for(int i = leftOffset_; i < wiLen_ - rightOffset_; i++)
        {
            velSendBuf_(i,k,ind::Z) = w_(i,wjLen_-2,k);
        }
    }
    mpiHandler_.send(velSendBuf_.data(), velSendBuf_.length());
}

void AsyncNeighbourTop::setRecvUVW()
{
    #pragma omp simd collapse(2)
    for(int k = hindOffset_; k < ukLen_ - frontOffset_; k++)
    {
        for(int i = leftOffset_; i < uiLen_ - rightOffset_; i++)
        {
            u_(i,ujLen_-1,k) = velRecvBuf_(i,k,ind::X);
        }
    }
    #pragma omp simd collapse(2)
    for(int k = hindOffset_; k < vkLen_ - frontOffset_; k++)
    {
        for(int i = leftOffset_; i < viLen_ - rightOffset_; i++)
        {
            v_(i,vjLen_-1,k) = velRecvBuf_(i,k,ind::Y);
        }
    }
    #pragma omp simd collapse(2)
    for(int k = hindOffset_; k < wkLen_ - frontOffset_; k++)
    {
        for(int i = leftOffset_; i < wiLen_ - rightOffset_; i++)
        {
            w_(i,wjLen_-1,k) = velRecvBuf_(i,k,ind::Z);
        }
    }
}

void AsyncNeighbourTop::exchangeFGH()
{
    mpiHandler_.startReceive(velRecvBuf_.data(), velRecvBuf_.length());
    #pragma omp simd collapse(2)
    for(int k = hindOffset_; k < ukLen_ - frontOffset_; k++)
    {
        for(int i = leftOffset_; i < uiLen_ - rightOffset_; i++)
        {
            velSendBuf_(i,k,ind::X) = f_(i,ujLen_-2,k);
        }
    }
    #pragma omp simd collapse(2)
    for(int k = hindOffset_; k < vkLen_ - frontOffset_; k++)
    {
        for(int i = leftOffset_; i < viLen_ - rightOffset_; i++)
        {
            velSendBuf_(i,k,ind::Y) = g_(i,vjLen_-2,k);
        }
    }
    #pragma omp simd collapse(2)
    for(int k = hindOffset_; k < wkLen_ - frontOffset_; k++)
    {
        for(int i = leftOffset_; i < wiLen_ - rightOffset_; i++)
        {
            velSendBuf_(i,k,ind::Z) = h_(i,wjLen_-2,k);
        }
    }
    mpiHandler_.send(velSendBuf_.data(), velSendBuf_.length());
}

void AsyncNeighbourTop::setRecvFGH()
{
    #pragma omp simd collapse(2)
    for(int k = hindOffset_; k < ukLen_ - frontOffset_; k++)
    {
        for(int i = leftOffset_; i < uiLen_ - rightOffset_; i++)
        {
            f_(i,ujLen_-1,k) = velRecvBuf_(i,k,ind::X);
        }
    }
    #pragma omp simd collapse(2)
    for(int k = hindOffset_; k < vkLen_ - frontOffset_; k++)
    {
        for(int i = leftOffset_; i < viLen_ - rightOffset_; i++)
        {
            g_(i,vjLen_-1,k) = velRecvBuf_(i,k,ind::Y);
        }
    }
    #pragma omp simd collapse(2)
    for(int k = hindOffset_; k < wkLen_ - frontOffset_; k++)
    {
        for(int i = leftOffset_; i < wiLen_ - rightOffset_; i++)
        {
            h_(i,wjLen_-1,k) = velRecvBuf_(i,k,ind::Z);
        }
    }
}

void AsyncNeighbourRight::exchangeP()
{
    mpiHandler_.startReceive(pRecvBuf_.data(), pRecvBuf_.length());
    #pragma omp simd collapse(2)
    for(int k = 0; k < pkLen_; k++)
    {
        for(int j = 0; j < pjLen_; j++)
        {
            pSendBuf_(j,k) = p_(piLen_-2,j,k);
        }
    }
    mpiHandler_.send(pSendBuf_.data(), pSendBuf_.length());
}

void AsyncNeighbourRight::setRecvP()
{
    #pragma omp simd collapse(2)
    for(int k = 0; k < pkLen_; k++)
    {
        for(int j = 0; j < pjLen_; j++)
        {
            p_(piLen_-1,j,k) = pRecvBuf_(j,k);
        }
    }
}

void AsyncNeighbourRight::exchangeUVW()
{
    mpiHandler_.startReceive(velRecvBuf_.data(), velRecvBuf_.length());
    #pragma omp simd collapse(2)
    for(int k = 0; k < ukLen_; k++)
    {
        for(int j = 0; j < ujLen_; j++)
        {
            velSendBuf_(j,k,ind::X) = u_(viLen_-2,j,k);
        }
    }
    #pragma omp simd collapse(2)
    for(int k = 0; k < vkLen_; k++)
    {
        for(int j = 0; j < vjLen_; j++)
        {
            velSendBuf_(j,k,ind::Y) = v_(viLen_-2,j,k);
        }
    }
    #pragma omp simd collapse(2)
    for(int k = 0; k < wkLen_; k++)
    {
        for(int j = 0; j < wjLen_; j++)
        {
            velSendBuf_(j,k,ind::Z) = w_(viLen_-2,j,k);
        }
    }
    mpiHandler_.send(velSendBuf_.data(), velSendBuf_.length());
}

void AsyncNeighbourRight::setRecvUVW()
{
    #pragma omp simd collapse(2)
    for(int k = 0; k < ukLen_; k++)
    {
        for(int j = 0; j < ujLen_; j++)
        {
            u_(viLen_-1,j,k) = velRecvBuf_(j,k,ind::X);
        }
    }
    #pragma omp simd collapse(2)
    for(int k = 0; k < vkLen_; k++)
    {
        for(int j = 0; j < vjLen_; j++)
        {
            v_(viLen_-1,j,k) = velRecvBuf_(j,k,ind::Y);
        }
    }
    #pragma omp simd collapse(2)
    for(int k = 0; k < wkLen_; k++)
    {
        for(int j = 0; j < wjLen_; j++)
        {
            w_(viLen_-1,j,k) = velRecvBuf_(j,k,ind::Z);
        }
    }
}

void AsyncNeighbourRight::exchangeFGH()
{
    mpiHandler_.startReceive(velRecvBuf_.data(), velRecvBuf_.length());
    #pragma omp simd collapse(2)
    for(int k = 0; k < ukLen_; k++)
    {
        for(int j = 0; j < ujLen_; j++)
        {
            velSendBuf_(j,k,ind::X) = f_(viLen_-2,j,k);
        }
    }
    #pragma omp simd collapse(2)
    for(int k = 0; k < vkLen_; k++)
    {
        for(int j = 0; j < vjLen_; j++)
        {
            velSendBuf_(j,k,ind::Y) = g_(viLen_-2,j,k);
        }
    }
    #pragma omp simd collapse(2)
    for(int k = 0; k < wkLen_; k++)
    {
        for(int j = 0; j < wjLen_; j++)
        {
            velSendBuf_(j,k,ind::Z) = h_(viLen_-2,j,k);
        }
    }
    mpiHandler_.send(velSendBuf_.data(), velSendBuf_.length());
}

void AsyncNeighbourRight::setRecvFGH()
{
    #pragma omp simd collapse(2)
    for(int k = 0; k < ukLen_; k++)
    {
        for(int j = 0; j < ujLen_; j++)
        {
            f_(viLen_-1,j,k) = velRecvBuf_(j,k,ind::X);
        }
    }
    #pragma omp simd collapse(2)
    for(int k = 0; k < vkLen_; k++)
    {
        for(int j = 0; j < vjLen_; j++)
        {
            g_(viLen_-1,j,k) = velRecvBuf_(j,k,ind::Y);
        }
    }
    #pragma omp simd collapse(2)
    for(int k = 0; k < wkLen_; k++)
    {
        for(int j = 0; j < wjLen_; j++)
        {
            h_(viLen_-1,j,k) = velRecvBuf_(j,k,ind::Z);
        }
    }
}

void AsyncNeighbourBottom::exchangeP()
{
    mpiHandler_.startReceive(pRecvBuf_.data(), pRecvBuf_.length());
    #pragma omp simd collapse(2)
    for(int k = hindOffset_; k < pkLen_ - frontOffset_; k++)
    {
        for(int i = leftOffset_; i < piLen_ - rightOffset_; i++)
        {
            pSendBuf_(i,k) = p_(i,1,k);
        }
    }
    mpiHandler_.send(pSendBuf_.data(), pSendBuf_.length());
}

void AsyncNeighbourBottom::setRecvP()
{
    #pragma omp simd collapse(2)
    for(int k = hindOffset_; k < pkLen_ - frontOffset_; k++)
    {
        for(int i = leftOffset_; i < piLen_ - rightOffset_; i++)
        {
            p_(i,0,k) = pRecvBuf_(i,k);
        }
    }
}

void AsyncNeighbourBottom::exchangeUVW()
{
    mpiHandler_.startReceive(velRecvBuf_.data(), velRecvBuf_.length());
    #pragma omp simd collapse(2)
    for(int k = hindOffset_; k < ukLen_ - frontOffset_; k++)
    {
        for(int i = leftOffset_; i < uiLen_ - rightOffset_; i++)
        {
            velSendBuf_(i,k,ind::X) = u_(i,1,k);
        }
    }
    #pragma omp simd collapse(2)
    for(int k = hindOffset_; k < vkLen_ - frontOffset_; k++)
    {
        for(int i = leftOffset_; i < viLen_ - rightOffset_; i++)
        {
            velSendBuf_(i,k,ind::Y) = v_(i,1,k);
        }
    }
    #pragma omp simd collapse(2)
    for(int k = hindOffset_; k < wkLen_ - frontOffset_; k++)
    {
        for(int i = leftOffset_; i < wiLen_ - rightOffset_; i++)
        {
            velSendBuf_(i,k,ind::Z) = w_(i,1,k);
        }
    }
    mpiHandler_.send(velSendBuf_.data(), velSendBuf_.length());
}

void AsyncNeighbourBottom::setRecvUVW()
{
    #pragma omp simd collapse(2)
    for(int k = hindOffset_; k < ukLen_ - frontOffset_; k++)
    {
        for(int i = leftOffset_; i < uiLen_ - rightOffset_; i++)
        {
            u_(i,0,k) = velRecvBuf_(i,k,ind::X);
        }
    }
    #pragma omp simd collapse(2)
    for(int k = hindOffset_; k < vkLen_ - frontOffset_; k++)
    {
        for(int i = leftOffset_; i < viLen_ - rightOffset_; i++)
        {
            v_(i,0,k) = velRecvBuf_(i,k,ind::Y);
        }
    }
    #pragma omp simd collapse(2)
    for(int k = hindOffset_; k < wkLen_ - frontOffset_; k++)
    {
        for(int i = leftOffset_; i < wiLen_ - rightOffset_; i++)
        {
            w_(i,0,k) = velRecvBuf_(i,k,ind::Z);
        }
    }
}

void AsyncNeighbourBottom::exchangeFGH()
{
    mpiHandler_.startReceive(velRecvBuf_.data(), velRecvBuf_.length());
    #pragma omp simd collapse(2)
    for(int k = hindOffset_; k < ukLen_ - frontOffset_; k++)
    {
        for(int i = leftOffset_; i < uiLen_ - rightOffset_; i++)
        {
            velSendBuf_(i,k,ind::X) = f_(i,1,k);
        }
    }
    #pragma omp simd collapse(2)
    for(int k = hindOffset_; k < vkLen_ - frontOffset_; k++)
    {
        for(int i = leftOffset_; i < viLen_ - rightOffset_; i++)
        {
            velSendBuf_(i,k,ind::Y) = g_(i,1,k);
        }
    }
    #pragma omp simd collapse(2)
    for(int k = hindOffset_; k < wkLen_ - frontOffset_; k++)
    {
        for(int i = leftOffset_; i < wiLen_ - rightOffset_; i++)
        {
            velSendBuf_(i,k,ind::Z) = h_(i,1,k);
        }
    }
    mpiHandler_.send(velSendBuf_.data(), velSendBuf_.length());
}

void AsyncNeighbourBottom::setRecvFGH()
{
    #pragma omp simd collapse(2)
    for(int k = hindOffset_; k < ukLen_ - frontOffset_; k++)
    {
        for(int i = leftOffset_; i < uiLen_ - rightOffset_; i++)
        {
            f_(i,0,k) = velRecvBuf_(i,k,ind::X);
        }
    }
    #pragma omp simd collapse(2)
    for(int k = hindOffset_; k < vkLen_ - frontOffset_; k++)
    {
        for(int i = leftOffset_; i < viLen_ - rightOffset_; i++)
        {
            g_(i,0,k) = velRecvBuf_(i,k,ind::Y);
        }
    }
    #pragma omp simd collapse(2)
    for(int k = hindOffset_; k < wkLen_ - frontOffset_; k++)
    {
        for(int i = leftOffset_; i < wiLen_ - rightOffset_; i++)
        {
            h_(i,0,k) = velRecvBuf_(i,k,ind::Z);
        }
    }
}

void AsyncNeighbourLeft::exchangeP()
{
    mpiHandler_.startReceive(pRecvBuf_.data(), pRecvBuf_.length());
    #pragma omp simd collapse(2)
    for(int k = 0; k < pkLen_; k++)
    {
        for(int j = 0; j < pjLen_; j++)
        {
            pSendBuf_(j,k) = p_(1,j,k);
        }
    }
    mpiHandler_.send(pSendBuf_.data(), pSendBuf_.length());
}

void AsyncNeighbourLeft::setRecvP()
{
    #pragma omp simd collapse(2)
    for(int k = 0; k < pkLen_; k++)
    {
        for(int j = 0; j < pjLen_; j++)
        {
            p_(0,j,k) = pRecvBuf_(j,k);
        }
    }
}

void AsyncNeighbourLeft::exchangeUVW()
{
    mpiHandler_.startReceive(velRecvBuf_.data(), velRecvBuf_.length());
    #pragma omp simd collapse(2)
    for(int k = 0; k < ukLen_; k++)
    {
        for(int j = 0; j < ujLen_; j++)
        {
            velSendBuf_(j,k,ind::X) = u_(1,j,k);
        }
    }
    #pragma omp simd collapse(2)
    for(int k = 0; k < vkLen_; k++)
    {
        for(int j = 0; j < vjLen_; j++)
        {
            velSendBuf_(j,k,ind::Y) = v_(1,j,k);
        }
    }
    #pragma omp simd collapse(2)
    for(int k = 0; k < wkLen_; k++)
    {
        for(int j = 0; j < wjLen_; j++)
        {
            velSendBuf_(j,k,ind::Z) = w_(1,j,k);
        }
    }
    mpiHandler_.send(velSendBuf_.data(), velSendBuf_.length());
}

void AsyncNeighbourLeft::setRecvUVW()
{
    #pragma omp simd collapse(2)
    for(int k = 0; k < ukLen_; k++)
    {
        for(int j = 0; j < ujLen_; j++)
        {
            u_(0,j,k) = velRecvBuf_(j,k,ind::X);
        }
    }
    #pragma omp simd collapse(2)
    for(int k = 0; k < vkLen_; k++)
    {
        for(int j = 0; j < vjLen_; j++)
        {
            v_(0,j,k) = velRecvBuf_(j,k,ind::Y);
        }
    }
    #pragma omp simd collapse(2)
    for(int k = 0; k < wkLen_; k++)
    {
        for(int j = 0; j < wjLen_; j++)
        {
            w_(0,j,k) = velRecvBuf_(j,k,ind::Z);
        }
    }
}

void AsyncNeighbourLeft::exchangeFGH()
{
    mpiHandler_.startReceive(velRecvBuf_.data(), velRecvBuf_.length());
    #pragma omp simd collapse(2)
    for(int k = 0; k < ukLen_; k++)
    {
        for(int j = 0; j < ujLen_; j++)
        {
            velSendBuf_(j,k,ind::X) = f_(1,j,k);
        }
    }
    #pragma omp simd collapse(2)
    for(int k = 0; k < vkLen_; k++)
    {
        for(int j = 0; j < vjLen_; j++)
        {
            velSendBuf_(j,k,ind::Y) = g_(1,j,k);
        }
    }
    #pragma omp simd collapse(2)
    for(int k = 0; k < wkLen_; k++)
    {
        for(int j = 0; j < wjLen_; j++)
        {
            velSendBuf_(j,k,ind::Z) = h_(1,j,k);
        }
    }
    mpiHandler_.send(velSendBuf_.data(), velSendBuf_.length());
}

void AsyncNeighbourLeft::setRecvFGH()
{
    #pragma omp simd collapse(2)
    for(int k = 0; k < ukLen_; k++)
    {
        for(int j = 0; j < ujLen_; j++)
        {
            f_(0,j,k) = velRecvBuf_(j,k,ind::X);
        }
    }
    #pragma omp simd collapse(2)
    for(int k = 0; k < vkLen_; k++)
    {
        for(int j = 0; j < vjLen_; j++)
        {
            g_(0,j,k) = velRecvBuf_(j,k,ind::Y);
        }
    }
    #pragma omp simd collapse(2)
    for(int k = 0; k < wkLen_; k++)
    {
        for(int j = 0; j < wjLen_; j++)
        {
            h_(0,j,k) = velRecvBuf_(j,k,ind::Z);
        }
    }
}

void AsyncNeighbourHind::exchangeP()
{
    mpiHandler_.startReceive(pRecvBuf_.data(), pRecvBuf_.length());
    #pragma omp simd collapse(2)
    for(int j = 0; j < pjLen_; j++)
    {
        for(int i = leftOffset_; i < piLen_ - rightOffset_; i++)
        {
            pSendBuf_(i,j) = p_(i,j,1);
        }
    }
    mpiHandler_.send(pSendBuf_.data(), pSendBuf_.length());
}

void AsyncNeighbourHind::setRecvP()
{
    #pragma omp simd collapse(2)
    for(int j = 0; j < pjLen_; j++)
    {
        for(int i = leftOffset_; i < piLen_ - rightOffset_; i++)
        {
            p_(i,j,0) = pRecvBuf_(i,j);
        }
    }
}

void AsyncNeighbourHind::exchangeUVW()
{
    mpiHandler_.startReceive(velRecvBuf_.data(), velRecvBuf_.length());
    #pragma omp simd collapse(2)
    for(int j = 0; j < ujLen_; j++)
    {
        for(int i = leftOffset_; i < uiLen_ - rightOffset_; i++)
        {
            velSendBuf_(i,j,ind::X) = u_(i,j,1);
        }
    }
    #pragma omp simd collapse(2)
    for(int j = 0; j < vjLen_; j++)
    {
        for(int i = leftOffset_; i < viLen_ - rightOffset_; i++)
        {
            velSendBuf_(i,j,ind::Y) = v_(i,j,1);
        }
    }
    #pragma omp simd collapse(2)
    for(int j = 0; j < wjLen_; j++)
    {
        for(int i = leftOffset_; i < wiLen_ - rightOffset_; i++)
        {
            velSendBuf_(i,j,ind::Z) = w_(i,j,1);
        }
    }
    mpiHandler_.send(velSendBuf_.data(), velSendBuf_.length());
}

void AsyncNeighbourHind::setRecvUVW()
{
    #pragma omp simd collapse(2)
    for(int j = 0; j < ujLen_; j++)
    {
        for(int i = leftOffset_; i < uiLen_ - rightOffset_; i++)
        {
            u_(i,j,0) = velRecvBuf_(i,j,ind::X);
        }
    }
    #pragma omp simd collapse(2)
    for(int j = 0; j < vjLen_; j++)
    {
        for(int i = leftOffset_; i < viLen_ - rightOffset_; i++)
        {
            v_(i,j,0) = velRecvBuf_(i,j,ind::Y);
        }
    }
    #pragma omp simd collapse(2)
    for(int j = 0; j < wjLen_; j++)
    {
        for(int i = leftOffset_; i < wiLen_ - rightOffset_; i++)
        {
            w_(i,j,0) = velRecvBuf_(i,j,ind::Z);
        }
    }
}

void AsyncNeighbourHind::exchangeFGH()
{
    mpiHandler_.startReceive(velRecvBuf_.data(), velRecvBuf_.length());
    #pragma omp simd collapse(2)
    for(int j = 0; j < ujLen_; j++)
    {
        for(int i = leftOffset_; i < uiLen_ - rightOffset_; i++)
        {
            velSendBuf_(i,j,ind::X) = f_(i,j,1);
        }
    }
    #pragma omp simd collapse(2)
    for(int j = 0; j < vjLen_; j++)
    {
        for(int i = leftOffset_; i < viLen_ - rightOffset_; i++)
        {
            velSendBuf_(i,j,ind::Y) = g_(i,j,1);
        }
    }
    #pragma omp simd collapse(2)
    for(int j = 0; j < wjLen_; j++)
    {
        for(int i = leftOffset_; i < wiLen_ - rightOffset_; i++)
        {
            velSendBuf_(i,j,ind::Z) = h_(i,j,1);
        }
    }
    mpiHandler_.send(velSendBuf_.data(), velSendBuf_.length());
}

void AsyncNeighbourHind::setRecvFGH()
{
    #pragma omp simd collapse(2)
    for(int j = 0; j < ujLen_; j++)
    {
        for(int i = leftOffset_; i < uiLen_ - rightOffset_; i++)
        {
            f_(i,j,0) = velRecvBuf_(i,j,ind::X);
        }
    }
    #pragma omp simd collapse(2)
    for(int j = 0; j < vjLen_; j++)
    {
        for(int i = leftOffset_; i < viLen_ - rightOffset_; i++)
        {
            g_(i,j,0) = velRecvBuf_(i,j,ind::Y);
        }
    }
    #pragma omp simd collapse(2)
    for(int j = 0; j < wjLen_; j++)
    {
        for(int i = leftOffset_; i < wiLen_ - rightOffset_; i++)
        {
            h_(i,j,0) = velRecvBuf_(i,j,ind::Z);
        }
    }
}

void AsyncNeighbourFront::exchangeP()
{
    mpiHandler_.startReceive(pRecvBuf_.data(), pRecvBuf_.length());
    #pragma omp simd collapse(2)
    for(int j = 0; j < pjLen_; j++)
    {
        for(int i = leftOffset_; i < piLen_ - rightOffset_; i++)
        {
            pSendBuf_(i,j) = p_(i,j,pkLen_-2);
        }
    }
    mpiHandler_.send(pSendBuf_.data(), pSendBuf_.length());
}

void AsyncNeighbourFront::setRecvP()
{
    #pragma omp simd collapse(2)
    for(int j = 0; j < pjLen_; j++)
    {
        for(int i = leftOffset_; i < piLen_ - rightOffset_; i++)
        {
            p_(i,j,pkLen_-1) = pRecvBuf_(i,j);
        }
    }
}

void AsyncNeighbourFront::exchangeUVW()
{
    mpiHandler_.startReceive(velRecvBuf_.data(), velRecvBuf_.length());
    #pragma omp simd collapse(2)
    for(int j = 0; j < ujLen_; j++)
    {
        for(int i = leftOffset_; i < uiLen_ - rightOffset_; i++)
        {
            velSendBuf_(i,j,ind::X) = u_(i,j,ukLen_-2);
        }
    }
    #pragma omp simd collapse(2)
    for(int j = 0; j < vjLen_; j++)
    {
        for(int i = leftOffset_; i < viLen_ - rightOffset_; i++)
        {
            velSendBuf_(i,j,ind::Y) = v_(i,j,vkLen_-2);
        }
    }
    #pragma omp simd collapse(2)
    for(int j = 0; j < wjLen_; j++)
    {
        for(int i = leftOffset_; i < wiLen_ - rightOffset_; i++)
        {
            velSendBuf_(i,j,ind::Z) = w_(i,j,wkLen_-2);
        }
    }
    mpiHandler_.send(velSendBuf_.data(), velSendBuf_.length());
}

void AsyncNeighbourFront::setRecvUVW()
{
    #pragma omp simd collapse(2)
    for(int j = 0; j < ujLen_; j++)
    {
        for(int i = leftOffset_; i < uiLen_ - rightOffset_; i++)
        {
            u_(i,j,ukLen_-1) = velRecvBuf_(i,j,ind::X);
        }
    }
    #pragma omp simd collapse(2)
    for(int j = 0; j < vjLen_; j++)
    {
        for(int i = leftOffset_; i < viLen_ - rightOffset_; i++)
        {
            v_(i,j,vkLen_-1) = velRecvBuf_(i,j,ind::Y);
        }
    }
    #pragma omp simd collapse(2)
    for(int j = 0; j < wjLen_; j++)
    {
        for(int i = leftOffset_; i < wiLen_ - rightOffset_; i++)
        {
            w_(i,j,wkLen_-1) = velRecvBuf_(i,j,ind::Z);
        }
    }
}

void AsyncNeighbourFront::exchangeFGH()
{
    mpiHandler_.startReceive(velRecvBuf_.data(), velRecvBuf_.length());
    #pragma omp simd collapse(2)
    for(int j = 0; j < ujLen_; j++)
    {
        for(int i = leftOffset_; i < uiLen_ - rightOffset_; i++)
        {
            velSendBuf_(i,j,ind::X) = f_(i,j,ukLen_-2);
        }
    }
    #pragma omp simd collapse(2)
    for(int j = 0; j < vjLen_; j++)
    {
        for(int i = leftOffset_; i < viLen_ - rightOffset_; i++)
        {
            velSendBuf_(i,j,ind::Y) = g_(i,j,vkLen_-2);
        }
    }
    #pragma omp simd collapse(2)
    for(int j = 0; j < wjLen_; j++)
    {
        for(int i = leftOffset_; i < wiLen_ - rightOffset_; i++)
        {
            velSendBuf_(i,j,ind::Z) = h_(i,j,wkLen_-2);
        }
    }
    mpiHandler_.send(velSendBuf_.data(), velSendBuf_.length());
}

void AsyncNeighbourFront::setRecvFGH()
{
    #pragma omp simd collapse(2)
    for(int j = 0; j < ujLen_; j++)
    {
        for(int i = leftOffset_; i < uiLen_ - rightOffset_; i++)
        {
            f_(i,j,ukLen_-1) = velRecvBuf_(i,j,ind::X);
        }
    }
    #pragma omp simd collapse(2)
    for(int j = 0; j < vjLen_; j++)
    {
        for(int i = leftOffset_; i < viLen_ - rightOffset_; i++)
        {
            g_(i,j,vkLen_-1) = velRecvBuf_(i,j,ind::Y);
        }
    }
    #pragma omp simd collapse(2)
    for(int j = 0; j < wjLen_; j++)
    {
        for(int i = leftOffset_; i < wiLen_ - rightOffset_; i++)
        {
            h_(i,j,wkLen_-1) = velRecvBuf_(i,j,ind::Z);
        }
    }
}