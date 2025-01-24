#include "boundary/dirichlet.h"

void DirichletNorth::setUV()
{
    #pragma omp simd
    for(int i = 1; i < uiLen_-1; i++) 
    {
        u_(i,ujLen_-1) = 2*velX_ - u_(i,ujLen_-2);
    }
    #pragma omp simd
    for(int i = 1; i < viLen_-1; i++) 
    {
        v_(i,vjLen_-1) = velY_;
    }
}

// In case of Dirichlet, setting FG and UV may be done in one run,
// but for neighbour this has to be split
void DirichletNorth::setFG()
{
    #pragma omp simd
    for(int i = 1; i < uiLen_-1; i++) 
    {
        f_(i,ujLen_-1) = u_(i,ujLen_-1);
    }
    #pragma omp simd
    for(int i = 1; i < viLen_-1; i++) 
    {
        g_(i,vjLen_-1) = v_(i,vjLen_-1);
    }
}

void DirichletNorth::setP()
{
    #pragma omp simd
    for(int i = 1; i < piLen_-1; i++) 
    {
        p_(i,pjLen_-1) = p_(i,pjLen_-2);
    }
}

void DirichletEast::setUV()
{
    #pragma omp simd
    for(int j = 0; j < ujLen_; j++) 
    {
        u_(uiLen_-1,j) = velX_;
    }
    #pragma omp simd
    for(int j = 0; j < vjLen_; j++)
    {
        v_(viLen_-1,j) = 2*velY_ - v_(viLen_-2,j);
    }
}

void DirichletEast::setFG()
{
    #pragma omp simd
    for(int j = 0; j < ujLen_; j++) 
    {
        f_(uiLen_-1,j) = u_(uiLen_-1,j);
    }
    #pragma omp simd
    for(int j = 0; j < vjLen_; j++)
    {
        g_(viLen_-1,j) = v_(viLen_-1,j);
    }
}

void DirichletEast::setP()
{
    #pragma omp simd
    for(int j = 0; j < pjLen_; j++) 
    {
        p_(piLen_-1,j) = p_(piLen_-2,j);
    }
}

void DirichletSouth::setUV()
{
    #pragma omp simd
    for(int i = 1; i < uiLen_-1; i++) 
    {
        u_(i,0) = 2*velX_ - u_(i,1);
    }
    #pragma omp simd
    for(int i = 1; i < viLen_-1; i++) 
    {
        v_(i,0) = velY_;
    }
}

void DirichletSouth::setFG()
{
    #pragma omp simd
    for(int i = 1; i < uiLen_-1; i++) 
    {
        f_(i,0) = u_(i,0);
    }
    #pragma omp simd
    for(int i = 1; i < viLen_-1; i++) 
    {
        g_(i,0) = v_(i,0);
    }
}

void DirichletSouth::setP()
{
    #pragma omp simd
    for(int i = 1; i < piLen_-1; i++) 
    {
        p_(i,0) = p_(i,1);
    }
}

void DirichletWest::setUV()
{
    #pragma omp simd
    for(int j = 0; j < ujLen_; j++) 
    {
        u_(0,j) = velX_;
    }
    #pragma omp simd
    for(int j = 0; j < vjLen_; j++)
    {
        v_(0,j) = 2*velY_ - v_(1,j);
    }
}

void DirichletWest::setFG()
{
    #pragma omp simd
    for(int j = 0; j < ujLen_; j++) 
    {
        f_(0,j) = u_(0,j);
    }
    #pragma omp simd
    for(int j = 0; j < vjLen_; j++)
    {
        g_(0,j) = v_(0,j);
    }
}

void DirichletWest::setP()
{
    #pragma omp simd
    for(int j=0; j < pjLen_; j++) 
    {
        p_(0,j) = p_(1,j);
    }
}