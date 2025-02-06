#pragma once

#include <array>
#include <mpi.h>
#include "discretization/discretization.h"
#include "boundary/boundary.h"

class OutflowTop : public DomainBoundary
{
public:
    //! Set the edge implicitely
    OutflowTop(std::shared_ptr<Discretization> discretization) :
                   DomainBoundary(discretization, BoundaryEdge::TOP) { };

    void setUVW() override;
    void setFGH() override;
    void setP() override;
};

class OutflowRight : public DomainBoundary
{
public:
    OutflowRight(std::shared_ptr<Discretization> discretization) :
                   DomainBoundary(discretization, BoundaryEdge::RIGHT) { };

    void setUVW() override;
    void setFGH() override;
    void setP() override;
};

class OutflowBottom : public DomainBoundary
{
public:
    OutflowBottom(std::shared_ptr<Discretization> discretization) :
                   DomainBoundary(discretization, BoundaryEdge::BOTTOM) { };

    void setUVW() override;
    void setFGH() override;
    void setP() override;
};

class OutflowLeft : public DomainBoundary
{
public:
    OutflowLeft(std::shared_ptr<Discretization> discretization) :
                   DomainBoundary(discretization, BoundaryEdge::LEFT) { };
                   
    void setUVW() override;
    void setFGH() override;
    void setP() override;
};

class OutflowFront : public DomainBoundary
{
public:
    OutflowFront(std::shared_ptr<Discretization> discretization) :
                   DomainBoundary(discretization, BoundaryEdge::FRONT) { };
                   
    void setUVW() override;
    void setFGH() override;
    void setP() override;
};

class OutflowHind : public DomainBoundary
{
public:
    OutflowHind(std::shared_ptr<Discretization> discretization) :
                   DomainBoundary(discretization, BoundaryEdge::HIND) { };
                   
    void setUVW() override;
    void setFGH() override;
    void setP() override;
};