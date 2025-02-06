#pragma once

#include <array>
#include <mpi.h>
#include "discretization/discretization.h"
#include "boundary/boundary.h"

class SlipTop : public DomainBoundary
{
public:
    //! Set the edge implicitely
    SlipTop(std::shared_ptr<Discretization> discretization) :
                   DomainBoundary(discretization, BoundaryEdge::TOP) { };

    void setUVW() override;
    void setFGH() override;
    void setP() override;
};

class SlipRight : public DomainBoundary
{
public:
    SlipRight(std::shared_ptr<Discretization> discretization) :
                   DomainBoundary(discretization, BoundaryEdge::TOP) { };

    void setUVW() override;
    void setFGH() override;
    void setP() override;
};

class SlipBottom : public DomainBoundary
{
public:
    SlipBottom(std::shared_ptr<Discretization> discretization) :
                   DomainBoundary(discretization, BoundaryEdge::BOTTOM) { };

    void setUVW() override;
    void setFGH() override;
    void setP() override;
};

class SlipLeft : public DomainBoundary
{
public:
    SlipLeft(std::shared_ptr<Discretization> discretization) :
                   DomainBoundary(discretization, BoundaryEdge::LEFT) { };
                   
    void setUVW() override;
    void setFGH() override;
    void setP() override;
};

class SlipFront : public DomainBoundary
{
public:
    SlipFront(std::shared_ptr<Discretization> discretization) :
                   DomainBoundary(discretization, BoundaryEdge::FRONT) { };
                   
    void setUVW() override;
    void setFGH() override;
    void setP() override;
};

class SlipHind : public DomainBoundary
{
public:
    SlipHind(std::shared_ptr<Discretization> discretization) :
                   DomainBoundary(discretization, BoundaryEdge::HIND) { };
                   
    void setUVW() override;
    void setFGH() override;
    void setP() override;
};