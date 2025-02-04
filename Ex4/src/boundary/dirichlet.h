#pragma once

#include <array>
#include <mpi.h>
#include "discretization/discretization.h"
#include "boundary/boundary.h"

class Dirichlet : public Boundary
{
public:
    // additionall to the discretization and edge, it also requires the known velocities
    Dirichlet(std::shared_ptr<Discretization> discretization, BoundaryEdge edge, std::array<double, 3> vel) : 
              Boundary(discretization, edge),
              velX_(vel[0]), velY_(vel[1]), velZ_(vel[2]) { }
    
    virtual void setUVW() = 0;

    virtual void setFGH() = 0;

    virtual void setP() = 0;

protected:
    //! direction of flow
    const double velX_;
    const double velY_;
    const double velZ_;
};

class DirichletTop : public Dirichlet
{
public:
    //! Set the edge implicitely
    DirichletTop(std::shared_ptr<Discretization> discretization, std::array<double, 3> vel) :
                   Dirichlet(discretization, BoundaryEdge::TOP, vel) { };

    void setUVW() override;
    void setFGH() override;
    void setP() override;
};

class DirichletRight : public Dirichlet
{
public:
    DirichletRight(std::shared_ptr<Discretization> discretization, std::array<double, 3> vel) :
                   Dirichlet(discretization, BoundaryEdge::RIGHT, vel) { };

    void setUVW() override;
    void setFGH() override;
    void setP() override;
};

class DirichletBottom : public Dirichlet
{
public:
    DirichletBottom(std::shared_ptr<Discretization> discretization, std::array<double, 3> vel) :
                   Dirichlet(discretization, BoundaryEdge::BOTTOM, vel) { };

    void setUVW() override;
    void setFGH() override;
    void setP() override;
};

class DirichletLeft : public Dirichlet
{
public:
    DirichletLeft(std::shared_ptr<Discretization> discretization, std::array<double, 3> vel) :
                   Dirichlet(discretization, BoundaryEdge::LEFT, vel) { };
                   
    void setUVW() override;
    void setFGH() override;
    void setP() override;
};

class DirichletFront : public Dirichlet
{
public:
    DirichletFront(std::shared_ptr<Discretization> discretization, std::array<double, 3> vel) :
                   Dirichlet(discretization, BoundaryEdge::FRONT, vel) { };
                   
    void setUVW() override;
    void setFGH() override;
    void setP() override;
};

class DirichletHind : public Dirichlet
{
public:
    DirichletHind(std::shared_ptr<Discretization> discretization, std::array<double, 3> vel) :
                   Dirichlet(discretization, BoundaryEdge::HIND, vel) { };
                   
    void setUVW() override;
    void setFGH() override;
    void setP() override;
};