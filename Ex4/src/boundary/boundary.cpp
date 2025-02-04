#include "boundary/boundary.h"

Boundary::Boundary(std::shared_ptr<Discretization> discretization, 
                   BoundaryEdge edge) : 
                   discretization_(discretization),
                   edge_(edge),
                   uiLen_(discretization->u().size()[0]), viLen_(discretization->v().size()[0]), wiLen_(discretization->w().size()[0]), piLen_(discretization->p().size()[0]),
                   ujLen_(discretization->u().size()[1]), vjLen_(discretization->v().size()[1]), wjLen_(discretization->w().size()[1]), pjLen_(discretization->p().size()[1]),
                   ukLen_(discretization->u().size()[2]), vkLen_(discretization->v().size()[2]), wkLen_(discretization->w().size()[2]), pkLen_(discretization->p().size()[2]),
                   u_(discretization->u()), v_(discretization->v()), w_(discretization->w()),
                   f_(discretization->f()), g_(discretization->g()), h_(discretization->h()),
                   p_(discretization->p())
{
    assert(discretization != nullptr);
}