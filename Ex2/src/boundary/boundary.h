#pragma once

#include <memory>
#include <exception>
#include "discretization/discretization.h"

//! Enum class to simply encode all four directions. Char means, that if serialized, it is human readible.
enum BoundaryEdge : char
{
    NORTH = 'N',
    EAST  = 'E',
    SOUTH = 'S',
    WEST  = 'W'
};

//! Boundary base-class. This has two layer of abstraction: 
//! First the boundary type itself, then secondly an implemention class for each direction on top of that.
class Boundary
{
public:
    Boundary(std::shared_ptr<Discretization> discretization, BoundaryEdge edge);

    //! used to have virtual setP, setUV and setFG methods respectively,
    //! but with async-neighbours, they are split into send&receive

    //! Descriptor of which edge to manipulate
    const BoundaryEdge edge_;

protected:
    //! The discretization on which the boundaries are to be set
    const std::shared_ptr<Discretization> discretization_;

    //! For simple access, each field-variable is saved as reference
    FieldVariable &u_;
    FieldVariable &v_;
    FieldVariable &f_;
    FieldVariable &g_;
    FieldVariable &p_;

    //! Length information used for indexing
    const int uiLen_;
    const int viLen_;
    const int ujLen_;
    const int vjLen_;
    const int piLen_;
    const int pjLen_;

};
