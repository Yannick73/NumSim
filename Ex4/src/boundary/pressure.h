#pragma once

#include <array>
#include <mpi.h>
#include "discretization/discretization.h"
#include "boundary/boundary.h"

class PressureTop : public DomainBoundary
{
public:
    //! Set the edge implicitely
    PressureTop(std::shared_ptr<Discretization> discretization, double pBoundary) :
                   DomainBoundary(discretization, BoundaryEdge::TOP), pBoundary_(pBoundary) { };

    void setUVW() override;
    void setFGH() override;
    void setP() override;

private:
    const double pBoundary_;
};

class PressureRight : public DomainBoundary
{
public:
    PressureRight(std::shared_ptr<Discretization> discretization, double pBoundary) :
                   DomainBoundary(discretization, BoundaryEdge::TOP), pBoundary_(pBoundary) { };

    void setUVW() override;
    void setFGH() override;
    void setP() override;
    
private:
    const double pBoundary_;
};

class PressureBottom : public DomainBoundary
{
public:
    PressureBottom(std::shared_ptr<Discretization> discretization, double pBoundary) :
                   DomainBoundary(discretization, BoundaryEdge::BOTTOM), pBoundary_(pBoundary) { };

    void setUVW() override;
    void setFGH() override;
    void setP() override;
    
private:
    const double pBoundary_;
};

class PressureLeft : public DomainBoundary
{
public:
    PressureLeft(std::shared_ptr<Discretization> discretization, double pBoundary) :
                   DomainBoundary(discretization, BoundaryEdge::LEFT), pBoundary_(pBoundary) { };
                   
    void setUVW() override;
    void setFGH() override;
    void setP() override;
    
private:
    const double pBoundary_;
};

class PressureFront : public DomainBoundary
{
public:
    PressureFront(std::shared_ptr<Discretization> discretization, double pBoundary) :
                   DomainBoundary(discretization, BoundaryEdge::FRONT), pBoundary_(pBoundary) { };
                   
    void setUVW() override;
    void setFGH() override;
    void setP() override;
    
private:
    const double pBoundary_;
};

class PressureHind : public DomainBoundary
{
public:
    PressureHind(std::shared_ptr<Discretization> discretization, double pBoundary) :
                   DomainBoundary(discretization, BoundaryEdge::HIND), pBoundary_(pBoundary) { };
                   
    void setUVW() override;
    void setFGH() override;
    void setP() override;
    
private:
    const double pBoundary_;
};