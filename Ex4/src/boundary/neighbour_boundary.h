#pragma once

#include <array>
#include <sstream>
#include <iostream>
#include <mpi.h>
#include "discretization/discretization.h"
#include "boundary/boundary.h"
#include "storage/field_variable.h"

//! data used by MPI as tags (well, it probably is not required)
enum data
{
    UV = 0,
    FG = 1,
    P  = 2,
    RA = 3
};

//! Virtual implementation for a boundary to neighbouring partition
class NeighbourBoundary : public Boundary
{
public:
    //! Requires edge and length information, child classes may set them implicitely
    NeighbourBoundary(std::shared_ptr<Discretization> discretization, 
                      BoundaryEdge edge, int neighbourRank, int velBufLen, int pBufLen) : 
                      Boundary(discretization, edge), neighbourRank_(neighbourRank),
                      pBufLen_(pBufLen), velBufLen_(velBufLen),
                      r_(discretization->r()), a_(discretization->a()) { }
    
    virtual void setUV() = 0;

    virtual void setFG() = 0;

    virtual void setP()  = 0;
    // exchange variables R and A between ranks for CG
    virtual void setRA() = 0;

protected:
    //! id used for MPI communication
    const int neighbourRank_;

    //! lengths for velocity and pressure fields
    const int pBufLen_;
    const int velBufLen_;
    //! variable access for CG specific variables (need exchange, but no fix-boundary)
    std::shared_ptr<FieldVariable> r_;
    std::shared_ptr<FieldVariable> a_;
};

//! Following classes work the same, but for the respective orientation
class NeighbourNorth : public NeighbourBoundary
{
public:
    //! Using the discretization and beeing North, this sets up the general NeighbourBandary object
    //! Both velocity directions are sent simultaniously, so add them up for the buffer. 
    //! It is north, so they are sliced in x-direction.
    //! north and south neighbours have the quirk, that they may not override the edge values of their partitions' east/west neighbours
    NeighbourNorth(std::shared_ptr<Discretization> discretization, int neighbourRank, bool westNeighbour, bool eastNeighbour) :
                   NeighbourBoundary(discretization, BoundaryEdge::NORTH, neighbourRank,
                   discretization->u().size()[0] + discretization->v().size()[0], discretization->p().size()[0]),
                   westOffset_(westNeighbour ? 1 : 0), eastOffset_(eastNeighbour ? 1 : 0) {  };
                   // ? means, if westNeighbour is true, then assign 1 as offeset index, else assign 0

    void setUV() override;
    void setFG() override;
    void setP()  override;
    void setRA() override;

private:
    const int westOffset_;
    const int eastOffset_;
};

class NeighbourEast : public NeighbourBoundary
{
public:
    NeighbourEast(std::shared_ptr<Discretization> discretization, int neighbourRank) :
                   NeighbourBoundary(discretization, BoundaryEdge::EAST, neighbourRank,
                   discretization->u().size()[1] + discretization->v().size()[1],
                   discretization->p().size()[1]) { };

    void setUV() override;
    void setFG() override;
    void setP()  override;
    void setRA() override;
};

class NeighbourSouth : public NeighbourBoundary
{
public:
    NeighbourSouth(std::shared_ptr<Discretization> discretization, int neighbourRank, bool westNeighbour, bool eastNeighbour) :
                   NeighbourBoundary(discretization, BoundaryEdge::SOUTH, neighbourRank,
                   discretization->u().size()[0] + discretization->v().size()[0], discretization->p().size()[0]),
                   westOffset_(westNeighbour ? 1 : 0), eastOffset_(eastNeighbour ? 1 : 0) {  };

    void setUV() override;
    void setFG() override;
    void setP()  override;
    void setRA() override;

private:
    const int westOffset_;
    const int eastOffset_;
};

class NeighbourWest : public NeighbourBoundary
{
public:
    NeighbourWest(std::shared_ptr<Discretization> discretization, int neighbourRank) :
                   NeighbourBoundary(discretization, BoundaryEdge::WEST, neighbourRank,
                   discretization->u().size()[1] + discretization->v().size()[1],
                   discretization->p().size()[1]) { };
                   
    void setUV() override;
    void setFG() override;
    void setP()  override;
    void setRA() override;
};