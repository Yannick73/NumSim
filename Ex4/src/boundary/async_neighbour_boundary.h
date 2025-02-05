#pragma once

#include <array>
#include <sstream>
#include <iostream>
#include <mpi.h>
#include "discretization/discretization.h"
#include "boundary/boundary.h"
#include "boundary/mpi_wrapper.h"
#include "storage/field_variable.h"
#include "storage/array3D.h"
#include "storage/array2D.h"

/*
    In a quick test using Scenario3 with 80x60 cells and the simple exchange
    (non-blocking send and direct receive using arbitrary edge order & release)
    resulted in 35s runtime on my machine. YD
    
*/

enum ind : int
{
    X = 0,
    Y = 1,
    Z = 2
};

//! Virtual implementation for a boundary to neighbouring partition
class AsyncNeighbourBoundary : public Boundary
{
public:
    //! Requires edge and length information, child classes may set them implicitely
    AsyncNeighbourBoundary(std::shared_ptr<Discretization> d, 
                      BoundaryEdge edge, int neighbourRank, std::array<int, 2> velBufLen, std::array<int, 2> pBufLen) : 
                      Boundary(d, edge), neighbourRank_(neighbourRank),
                      pBufLen_(pBufLen), velBufLen_(velBufLen),
                      mpiHandler_(neighbourRank),
                      // 3 for u,v,w
                      velSendBuf_({velBufLen[0], velBufLen[1], 3}), velRecvBuf_({velBufLen[0], velBufLen[1], 3}),
                      pSendBuf_(pBufLen), pRecvBuf_(pBufLen) { }
    
    //! sends UV data and setups receive for it
    virtual void exchangeUVW() = 0;
    virtual void exchangeFGH() = 0;
    virtual void exchangeP()   = 0;

    //! sets the boundary to the values received in the buffer
    virtual void setRecvUVW() = 0;
    virtual void setRecvFGH() = 0;
    virtual void setRecvP()   = 0;

    //! Communication handler
    MPI_Wrapper mpiHandler_;

protected:
    //! id used for MPI communication
    const int neighbourRank_;

    //! lengths for velocity and pressure fields
    const std::array<int, 2> velBufLen_;
    const std::array<int, 2> pBufLen_;

    //! buffers for exchange, uv and fg exchange the same amount of data, so reuse
    Array3D velSendBuf_;
    Array3D velRecvBuf_;
    Array2D pSendBuf_;
    Array2D pRecvBuf_;
    //! it may be a good idea, to wrap the mpi comm into a seperate class
};

//! given the two indices, this function returns the maximum size of the velocity field
inline std::array<int, 2> getVelSize(std::shared_ptr<Discretization> d, int i0, int i1)
{
    const int max_uv0  = std::max(d->u().size()[i0], d->v().size()[i0]);
    const int max_uvw0 = std::max(max_uv0, d->w().size()[i0]);
    const int max_uv1  = std::max(d->u().size()[i1], d->v().size()[i1]);
    const int max_uvw1 = std::max(max_uv1, d->w().size()[i1]);
    return {max_uvw0, max_uvw1};
}

inline std::array<int, 2> getPSize(std::shared_ptr<Discretization> d, int i0, int i1)
{
    return {d->p().size()[i0], d->p().size()[i1]};
}

//! Following classes work the same, but for the respective orientation
class AsyncNeighbourTop : public AsyncNeighbourBoundary
{
public:
    //! Using the d and beeing North, this sets up the general AsyncNeighbourBandary object
    //! Both velocity directions are sent simultaniously, so add them up for the buffer. 
    //! It is north, so they are sliced in x-direction.
    //! north and south neighbours have the quirk, that they may not override the edge values of their partitions' east/west neighbours
    AsyncNeighbourTop(std::shared_ptr<Discretization> d, int neighbourRank, bool leftNeighbour, bool rightNeighbour, bool hindNeighbour, bool frontNeighbour) :
                   AsyncNeighbourBoundary(d, BoundaryEdge::TOP, neighbourRank,
                   getVelSize(d, ind::X, ind::Z), getPSize(d, ind::X, ind::Z)),
                   leftOffset_(leftNeighbour ? 1 : 0), rightOffset_(rightNeighbour ? 1 : 0),
                   frontOffset_(frontNeighbour ? 1 : 0), hindOffset_(hindNeighbour ? 1 : 0) {  };
                   // ? means, if leftNeighbour is true, then assign 1, else assign 0

    void exchangeUVW() override;
    void exchangeFGH() override;
    void exchangeP()  override;

    void setRecvUVW() override;
    void setRecvFGH() override;
    void setRecvP()  override;

private:
    const int leftOffset_;
    const int rightOffset_;
    const int frontOffset_;
    const int hindOffset_;
};

class AsyncNeighbourRight : public AsyncNeighbourBoundary
{
public:
    AsyncNeighbourRight(std::shared_ptr<Discretization> d, int neighbourRank) :
                   AsyncNeighbourBoundary(d, BoundaryEdge::RIGHT, neighbourRank,
                   getVelSize(d, ind::Y, ind::Z), getPSize(d, ind::Y, ind::Z)) { };

    void exchangeUVW() override;
    void exchangeFGH() override;
    void exchangeP()  override;

    void setRecvUVW() override;
    void setRecvFGH() override;
    void setRecvP()  override;
};

class AsyncNeighbourBottom : public AsyncNeighbourBoundary
{
public:
    AsyncNeighbourBottom(std::shared_ptr<Discretization> d, int neighbourRank, bool leftNeighbour, bool rightNeighbour, bool hindNeighbour, bool frontNeighbour) :
                   AsyncNeighbourBoundary(d, BoundaryEdge::BOTTOM, neighbourRank,
                   getVelSize(d, ind::X, ind::Z), getPSize(d, ind::X, ind::Z)),
                   leftOffset_(leftNeighbour ? 1 : 0), rightOffset_(rightNeighbour ? 1 : 0),
                   hindOffset_(hindNeighbour ? 1 : 0), frontOffset_(frontNeighbour ? 1 : 0) { }

    void exchangeUVW() override;
    void exchangeFGH() override;
    void exchangeP()  override;

    void setRecvUVW() override;
    void setRecvFGH() override;
    void setRecvP()  override;

private:
    const int leftOffset_;
    const int rightOffset_;
    const int hindOffset_;
    const int frontOffset_;
};

class AsyncNeighbourLeft : public AsyncNeighbourBoundary
{
public:
    AsyncNeighbourLeft(std::shared_ptr<Discretization> d, int neighbourRank) :
                   AsyncNeighbourBoundary(d, BoundaryEdge::LEFT, neighbourRank,
                   getVelSize(d, ind::X, ind::Z), getPSize(d, ind::X, ind::Z)) { };
                   
    void exchangeUVW() override;
    void exchangeFGH() override;
    void exchangeP()  override;

    void setRecvUVW() override;
    void setRecvFGH() override;
    void setRecvP()  override;
};

class AsyncNeighbourHind : public AsyncNeighbourBoundary
{
public:
    AsyncNeighbourHind(std::shared_ptr<Discretization> d, int neighbourRank, bool leftNeighbour, bool rightNeighbour) :
                   AsyncNeighbourBoundary(d, BoundaryEdge::HIND, neighbourRank,
                   getVelSize(d, ind::X, ind::Y), getPSize(d, ind::X, ind::Y)),
                   leftOffset_(leftNeighbour ? 1 : 0), rightOffset_(rightNeighbour ? 1 : 0) { };
                   
    void exchangeUVW() override;
    void exchangeFGH() override;
    void exchangeP()  override;

    void setRecvUVW() override;
    void setRecvFGH() override;
    void setRecvP()  override;

private:
    int leftOffset_;
    int rightOffset_;
};

class AsyncNeighbourFront : public AsyncNeighbourBoundary
{
public:
    AsyncNeighbourFront(std::shared_ptr<Discretization> d, int neighbourRank, bool leftNeighbour, bool rightNeighbour) :
                   AsyncNeighbourBoundary(d, BoundaryEdge::FRONT, neighbourRank,
                   getVelSize(d, ind::X, ind::Y), getPSize(d, ind::X, ind::Y)),
                   leftOffset_(leftNeighbour ? 1 : 0), rightOffset_(rightNeighbour ? 1 : 0) { };
                   
    void exchangeUVW() override;
    void exchangeFGH() override;
    void exchangeP()  override;

    void setRecvUVW() override;
    void setRecvFGH() override;
    void setRecvP()  override;

private:
    int leftOffset_;
    int rightOffset_;
};