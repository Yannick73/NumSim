#pragma once

#include <array>
#include <mpi.h>
#include "discretization/discretization.h"
#include "boundary/boundary.h"

class Dirichlet : public Boundary
{
public:
    // additionall to the discretization and edge, it also requires the known velocities
    Dirichlet(std::shared_ptr<Discretization> discretization, BoundaryEdge edge, std::array<double, 2> vel) : 
              Boundary(discretization, edge),
              velX_(vel[0]), velY_(vel[1]) { }
    
    virtual void setUV() = 0;

    virtual void setFG() = 0;

    virtual void setP() = 0;

protected:
    //! direction of flow
    const double velX_;
    const double velY_;
};

class DirichletNorth : public Dirichlet
{
public:
    //! Set the edge implicitely
    DirichletNorth(std::shared_ptr<Discretization> discretization, std::array<double, 2> vel) :
                   Dirichlet(discretization, BoundaryEdge::NORTH, vel) { };

    void setUV() override;
    void setFG() override;
    void setP() override;
};

class DirichletEast : public Dirichlet
{
public:
    DirichletEast(std::shared_ptr<Discretization> discretization, std::array<double, 2> vel) :
                   Dirichlet(discretization, BoundaryEdge::EAST, vel) { };

    void setUV() override;
    void setFG() override;
    void setP() override;
};

class DirichletSouth : public Dirichlet
{
public:
    DirichletSouth(std::shared_ptr<Discretization> discretization, std::array<double, 2> vel) :
                   Dirichlet(discretization, BoundaryEdge::SOUTH, vel) { };

    void setUV() override;
    void setFG() override;
    void setP() override;
};

class DirichletWest : public Dirichlet
{
public:
    DirichletWest(std::shared_ptr<Discretization> discretization, std::array<double, 2> vel) :
                   Dirichlet(discretization, BoundaryEdge::WEST, vel) { };
                   
    void setUV() override;
    void setFG() override;
    void setP() override;
};