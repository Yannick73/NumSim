#pragma once

#include <memory>
#include <exception>
#include "discretization/discretization.h"

//! Enum class to simply encode all four directions. Char means, that if serialized, it is human readible.
enum BoundaryEdge : char
{
    TOP    = 'T',
    RIGHT  = 'R',
    BOTTOM = 'B',
    LEFT   = 'L',
    FRONT  = 'F',
    HIND   = 'H'
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
    FieldVariable &w_;
    FieldVariable &f_;
    FieldVariable &g_;
    FieldVariable &h_;
    FieldVariable &p_;

    //! Length information used for indexing
    const int uiLen_;
    const int viLen_;
    const int wiLen_;

    const int ujLen_;
    const int vjLen_;
    const int wjLen_;

    const int ukLen_;
    const int vkLen_;
    const int wkLen_;

    const int piLen_;
    const int pjLen_;
    const int pkLen_;
};
