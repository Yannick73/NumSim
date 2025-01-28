#pragma once

#include <array>
#include <sstream>
#include <iostream>
#include <mpi.h>
#include "discretization/discretization.h"
#include "boundary/boundary.h"
#include "boundary/mpi_wrapper.h"
#include "storage/field_variable.h"

/*
    In a quick test using Scenario3 with 80x60 cells and the simple exchange
    (non-blocking send and direct receive using arbitrary edge order & release)
    resulted in 35s runtime on my machine. YD
    
*/

//! Virtual implementation for a boundary to neighbouring partition
class AsyncNeighbourBoundary : public Boundary
{
public:
    //! Requires edge and length information, child classes may set them implicitely
    AsyncNeighbourBoundary(std::shared_ptr<Discretization> discretization, 
                      BoundaryEdge edge, int neighbourRank, int velBufLen, int pBufLen) : 
                      Boundary(discretization, edge), neighbourRank_(neighbourRank),
                      pBufLen_(pBufLen), velBufLen_(velBufLen),
                      // 2*pBufLen is used sending R&A with CG (which is not implemented yet)
                      maxBufLen_(std::max(velBufLen, 2*pBufLen)),
                      r_(discretization->r()), a_(discretization->a()), mpiHandler_(neighbourRank),
                      sendBuf_(maxBufLen_, 0.0), recvBuf_(maxBufLen_, 0.0) { }
    
    //! sends UV data and setups receive for it
    virtual void exchangeUV() = 0;
    virtual void exchangeFG() = 0;
    virtual void exchangeP()  = 0;
    // exchange variables R and A between ranks for CG
    virtual void exchangeRA() = 0;

    //! sets the boundary to the values received in the buffer
    virtual void setRecvUV() = 0;
    virtual void setRecvFG() = 0;
    virtual void setRecvP()  = 0;
    virtual void setRecvRA() = 0;

    //! Communication handler
    MPI_Wrapper mpiHandler_;

protected:
    //! id used for MPI communication
    const int neighbourRank_;

    //! lengths for velocity and pressure fields
    const int velBufLen_;
    const int pBufLen_;
    //! length of the longest buffer, probably not used besides the constructor
    const int maxBufLen_;

    //! buffers for exchange, uv and fg exchange the same amount of data, so reuse
    std::vector<double> sendBuf_;
    std::vector<double> recvBuf_;
    //! it may be a good idea, to wrap the mpi comm into a seperate class

    //! variable access for CG specific variables (need exchange, but no fix-boundary)
    std::shared_ptr<FieldVariable> r_;
    std::shared_ptr<FieldVariable> a_;
    // well, unused as of now
};

//! Following classes work the same, but for the respective orientation
class AsyncNeighbourNorth : public AsyncNeighbourBoundary
{
public:
    //! Using the discretization and beeing North, this sets up the general AsyncNeighbourBandary object
    //! Both velocity directions are sent simultaniously, so add them up for the buffer. 
    //! It is north, so they are sliced in x-direction.
    //! north and south neighbours have the quirk, that they may not override the edge values of their partitions' east/west neighbours
    AsyncNeighbourNorth(std::shared_ptr<Discretization> discretization, int neighbourRank, bool westNeighbour, bool eastNeighbour) :
                   AsyncNeighbourBoundary(discretization, BoundaryEdge::TOP, neighbourRank,
                   discretization->u().size()[0] + discretization->v().size()[0], discretization->p().size()[0]),
                   westOffset_(westNeighbour ? 1 : 0), eastOffset_(eastNeighbour ? 1 : 0) {  };
                   // ? means, if westNeighbour is true, then assign 1, else assign 0

    void exchangeUV() override;
    void exchangeFG() override;
    void exchangeP()  override;
    void exchangeRA() override;

    void setRecvUV() override;
    void setRecvFG() override;
    void setRecvP()  override;
    void setRecvRA() override;

private:
    const int westOffset_;
    const int eastOffset_;
};

class AsyncNeighbourEast : public AsyncNeighbourBoundary
{
public:
    AsyncNeighbourEast(std::shared_ptr<Discretization> discretization, int neighbourRank) :
                   AsyncNeighbourBoundary(discretization, BoundaryEdge::RIGHT, neighbourRank,
                   discretization->u().size()[1] + discretization->v().size()[1],
                   discretization->p().size()[1]) { };

    void exchangeUV() override;
    void exchangeFG() override;
    void exchangeP()  override;
    void exchangeRA() override;

    void setRecvUV() override;
    void setRecvFG() override;
    void setRecvP()  override;
    void setRecvRA() override;
};

class AsyncNeighbourSouth : public AsyncNeighbourBoundary
{
public:
    AsyncNeighbourSouth(std::shared_ptr<Discretization> discretization, int neighbourRank, bool westNeighbour, bool eastNeighbour) :
                   AsyncNeighbourBoundary(discretization, BoundaryEdge::BOTTOM, neighbourRank,
                   discretization->u().size()[0] + discretization->v().size()[0], discretization->p().size()[0]),
                   westOffset_(westNeighbour ? 1 : 0), eastOffset_(eastNeighbour ? 1 : 0) {  };

    void exchangeUV() override;
    void exchangeFG() override;
    void exchangeP()  override;
    void exchangeRA() override;

    void setRecvUV() override;
    void setRecvFG() override;
    void setRecvP()  override;
    void setRecvRA() override;

private:
    const int westOffset_;
    const int eastOffset_;
};

class AsyncNeighbourWest : public AsyncNeighbourBoundary
{
public:
    AsyncNeighbourWest(std::shared_ptr<Discretization> discretization, int neighbourRank) :
                   AsyncNeighbourBoundary(discretization, BoundaryEdge::LEFT, neighbourRank,
                   discretization->u().size()[1] + discretization->v().size()[1],
                   discretization->p().size()[1]) { };
                   
    void exchangeUV() override;
    void exchangeFG() override;
    void exchangeP()  override;
    void exchangeRA() override;

    void setRecvUV() override;
    void setRecvFG() override;
    void setRecvP()  override;
    void setRecvRA() override;
};