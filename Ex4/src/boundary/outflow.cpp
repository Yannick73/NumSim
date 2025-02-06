#include "boundary/outflow.h"

// All remain, but the normal is approximated using interpolation
void OutflowTop::setUVW()
{
    #pragma omp simd collapse(2)
    for(int k = 1; k < ukLen_-1; k++)
    {
        for(int i = 1; i < uiLen_-1; i++) 
        {
            u_(i,ujLen_-1,k) = u_(i,ujLen_-2,k);
        }
    }
    #pragma omp simd collapse(2)
    for(int k = 1; k < vkLen_-1; k++)
    {
        for(int i = 1; i < viLen_-1; i++) 
        {
            v_(i,vjLen_-1,k) = (v_(i,vjLen_-1,k)- v_(i,vjLen_-2,k)) / discretization_->dy();
        }
    }
    #pragma omp simd collapse(2)
    for(int k = 1; k < wkLen_-1; k++)
    {
        for(int i = 1; i < wiLen_-1; i++)
        {
            w_(i,wjLen_-1,k) = w_(i,wjLen_-2,k);
        }
    }
}

// All remain (I think)
void OutflowTop::setFGH()
{
    #pragma omp simd collapse(2)
    for(int k = 1; k < ukLen_-1; k++)
    {
        for(int i = 1; i < uiLen_-1; i++) 
        {
            f_(i,ujLen_-1,k) = u_(i,ujLen_-1,k);
        }
    }
    #pragma omp simd collapse(2)
    for(int k = 1; k < vkLen_-1; k++)
    {
        for(int i = 1; i < viLen_-1; i++) 
        {
            g_(i,vjLen_-1,k) = v_(i,vjLen_-1,k);
        }
    }
    #pragma omp simd collapse(2)
    for(int k = 1; k < wkLen_-1; k++)
    {
        for(int i = 1; i < wiLen_-1; i++)
        {
            h_(i,wjLen_-1,k) = w_(i,wjLen_-1,k);
        }
    }
}

void OutflowTop::setP()
{
    #pragma omp simd collapse(2)
    for(int k = 1; k < pkLen_-1; k++)
    {
        for(int i = 1; i < piLen_-1; i++) 
        {
            p_(i,pjLen_-1,k) = p_(i,pjLen_-2,k);
        }
    }
}

void OutflowRight::setUVW()
{
    #pragma omp simd collapse(2)
    for(int k = 0; k < ukLen_; k++)
    {
        for(int j = 0; j < ujLen_; j++) 
        {
            u_(uiLen_-1,j,k) = (u_(uiLen_-2,j,k) - u_(uiLen_-1,j,k)) / discretization_->dx();
        }
    }
    #pragma omp simd collapse(2)
    for(int k = 0; k < vkLen_; k++)
    {
        for(int j = 0; j < vjLen_; j++)
        {
            v_(viLen_-1,j,k) = v_(viLen_-2,j,k);
        }
    }
    #pragma omp simd collapse(2)
    for(int k = 0; k < wkLen_; k++)
    {
        for(int j = 0; j < wjLen_; j++)
        {
            w_(wiLen_-1,j,k) = w_(wiLen_-2,j,k);
        }
    }
}

void OutflowRight::setFGH()
{
    #pragma omp simd collapse(2)
    for(int k = 0; k < ukLen_; k++)
    {
        for(int j = 0; j < ujLen_; j++) 
        {
            f_(uiLen_-1,j,k) = u_(uiLen_-1,j,k);
        }
    }
    #pragma omp simd collapse(2)
    for(int k = 0; k < vkLen_; k++)
    {
        for(int j = 0; j < vjLen_; j++)
        {
            g_(viLen_-1,j,k) = v_(viLen_-1,j,k);
        }
    }
    #pragma omp simd collapse(2)
    for(int k = 0; k < wkLen_; k++)
    {
        for(int j = 0; j < wjLen_; j++)
        {
            h_(wiLen_-1,j,k) = w_(wiLen_-1,j,k);
        }
    }
}

void OutflowRight::setP()
{
    #pragma omp simd collapse(2)
    for(int k = 0; k < pkLen_; k++)
    {
        for(int j = 0; j < pjLen_; j++) 
        {
            p_(piLen_-1,j,k) = p_(piLen_-2,j,k);
        }
    }
}

void OutflowBottom::setUVW()
{
    #pragma omp simd collapse(2)
    for(int k = 1; k < ukLen_-1; k++)
    {
        for(int i = 1; i < uiLen_-1; i++) 
        {
            u_(i,0,k) = u_(i,1,k);
        }
    }
    #pragma omp simd collapse(2)
    for(int k = 1; k < vkLen_-1; k++)
    {
        for(int i = 1; i < viLen_-1; i++) 
        {
            v_(i,0,k) = (v_(i,1,k) - v_(i,0,k)) / discretization_->dy();
        }
    }
    #pragma omp simd collapse(2)
    for(int k = 1; k < wkLen_-1; k++)
    {
        for(int i = 1; i < wiLen_-1; i++) 
        {
            w_(i,0,k) = w_(i,1,k);
        }
    }
}

void OutflowBottom::setFGH()
{
    #pragma omp simd collapse(2)
    for(int k = 1; k < ukLen_-1; k++)
    {
        for(int i = 1; i < uiLen_-1; i++) 
        {
            f_(i,0,k) = u_(i,0,k);
        }
    }
    #pragma omp simd collapse(2)
    for(int k = 1; k < vkLen_-1; k++)
    {
        for(int i = 1; i < viLen_-1; i++) 
        {
            g_(i,0,k) = v_(i,0,k);
        }
    }
    #pragma omp simd collapse(2)
    for(int k = 1; k < wkLen_-1; k++)
    {
        for(int i = 1; i < wiLen_-1; i++)
        {
            h_(i,0,k) = w_(i,0,k);
        }
    }
}

void OutflowBottom::setP()
{
    #pragma omp simd collapse(2)
    for(int k = 1; k < pkLen_-1; k++)
    {
        for(int i = 1; i < piLen_-1; i++) 
        {
            p_(i,0,k) = p_(i,1,k);
        }
    }
}

void OutflowLeft::setUVW()
{
    #pragma omp simd collapse(2)
    for(int k = 0; k < ukLen_; k++)
    {
        for(int j = 0; j < ujLen_; j++) 
        {
            u_(0,j,k) = (u_(1,j,k) - u_(0,j,k)) / discretization_->dx();
        }
    }
    #pragma omp simd collapse(2)
    for(int k = 0; k < vkLen_; k++)
    {
        for(int j = 0; j < vjLen_; j++)
        {
            v_(0,j,k) = v_(1,j,k);
        }
    }
    #pragma omp simd collapse(2)
    for(int k = 0; k < wkLen_; k++)
    {
        for(int j = 0; j < wjLen_; j++)
        {
            w_(0,j,k) = w_(1,j,k);
        }
    }
}

void OutflowLeft::setFGH()
{
    #pragma omp simd collapse(2)
    for(int k = 0; k < ukLen_; k++)
    {
        for(int j = 0; j < ujLen_; j++) 
        {
            f_(0,j,k) = u_(0,j,k);
        }
    }
    #pragma omp simd collapse(2)
    for(int k = 0; k < vkLen_; k++)
    {
        for(int j = 0; j < vjLen_; j++)
        {
            g_(0,j,k) = v_(0,j,k);
        }
    }
    #pragma omp simd collapse(2)
    for(int k = 0; k < wkLen_; k++)
    {
        for(int j = 0; j < wjLen_; j++)
        {
            h_(0,j,k) = w_(0,j,k);
        }
    }
}

void OutflowLeft::setP()
{
    #pragma omp simd collapse(2)
    for(int k = 0; k < pkLen_; k++)
    {
        for(int j = 0; j < pjLen_; j++) 
        {
            p_(0,j,k) = p_(1,j,k);
        }
    }
}

void OutflowFront::setUVW()
{
    #pragma omp simd collapse(2)
    for(int j = 0; j < ujLen_; j++)
    {
        for(int i = 1; i < uiLen_-1; i++)
        {
            u_(i,j,ukLen_-1) = u_(i,j,ukLen_-2);
        }
    }
    #pragma omp simd collapse(2)
    for(int j = 0; j < vjLen_; j++)
    {
        for(int i = 1; i < viLen_-1; i++)
        {
            v_(i,j,vkLen_-1) = v_(i,j,vkLen_-2);
        }
    }
    #pragma omp simd collapse(2)
    for(int j = 0; j < wjLen_; j++)
    {
        for(int i = 1; i < wiLen_; i++)
        {
            w_(i,j,wkLen_-1) = (w_(i,j,wkLen_-2) - w_(i,j,wkLen_-1)) / discretization_->dz();
        }
    }
}

void OutflowFront::setFGH()
{
    #pragma omp simd collapse(2)
    for(int j = 0; j < ujLen_; j++)
    {
        for(int i = 1; i < uiLen_-1; i++)
        {
            // changed to same indices
            f_(i,j,ukLen_-1) = u_(i,j,ukLen_-1);
        }
    }
    #pragma omp simd collapse(2)
    for(int j = 0; j < vjLen_; j++)
    {
        for(int i = 1; i < viLen_-1; i++)
        {
            g_(i,j,vkLen_-1) = v_(i,j,vkLen_-1);
        }
    }
    #pragma omp simd collapse(2)
    for(int j = 0; j < wjLen_; j++)
    {
        for(int i = 1; i < wiLen_; i++)
        {
            h_(i,j,wkLen_-1) = w_(i,j,wkLen_-1);
        }
    }
}

void OutflowFront::setP()
{
    #pragma omp simd collapse(2)
    for(int j = 0; j < pjLen_; j++)
    {
        for(int i = 0; i < piLen_; i++)
        {
            p_(i,j,pkLen_-1) = p_(i,j,pkLen_-2);
        }
    }
}

void OutflowHind::setUVW()
{
    #pragma omp simd collapse(2)
    for(int j = 0; j < ujLen_; j++)
    {
        for(int i = 1; i < uiLen_-1; i++)
        {
            u_(i,j,0) = u_(i,j,1);
        }
    }
    #pragma omp simd collapse(2)
    for(int j = 0; j < vjLen_; j++)
    {
        for(int i = 1; i < viLen_-1; i++)
        {
            v_(i,j,0) = v_(i,j,1);
        }
    }
    #pragma omp simd collapse(2)
    for(int j = 0; j < wjLen_; j++)
    {
        for(int i = 1; i < wiLen_; i++)
        {
            w_(i,j,0) = (w_(i,j,1) - w_(i,j,0)) / discretization_->dz();
        }
    }
}

void OutflowHind::setFGH()
{
    #pragma omp simd collapse(2)
    for(int j = 0; j < ujLen_; j++)
    {
        for(int i = 1; i < uiLen_-1; i++)
        {
            f_(i,j,0) = u_(i,j,0);
        }
    }
    #pragma omp simd collapse(2)
    for(int j = 0; j < vjLen_; j++)
    {
        for(int i = 1; i < viLen_-1; i++)
        {
            g_(i,j,0) = v_(i,j,0);
        }
    }
    #pragma omp simd collapse(2)
    for(int j = 0; j < wjLen_; j++)
    {
        for(int i = 1; i < wiLen_; i++)
        {
            h_(i,j,0) = w_(i,j,0);
        }
    }
}

void OutflowHind::setP()
{
    #pragma omp simd collapse(2)
    for(int j = 0; j < pjLen_; j++)
    {
        for(int i = 0; i < piLen_; i++)
        {
            p_(i,j,0) = p_(i,j,1);
        }
    }
}