#include "boundary/boundary.h"

Boundary::Boundary(std::shared_ptr<Discretization> discretization, 
                   BoundaryEdge edge) : 
                   discretization_(discretization),
                   edge_(edge),
                   uiLen_(discretization->u().size()[0]), viLen_(discretization->v().size()[0]),
                   ujLen_(discretization->u().size()[1]), vjLen_(discretization->v().size()[1]),
                   piLen_(discretization->p().size()[0]), pjLen_(discretization->p().size()[1]),
                   u_(discretization->u()), v_(discretization->v()),
                   f_(discretization->f()), g_(discretization->g()),
                   p_(discretization->p())
{
    assert(discretization != nullptr);
}