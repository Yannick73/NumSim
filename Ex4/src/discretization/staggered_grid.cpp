#include "discretization/staggered_grid.h"

StaggeredGrid::StaggeredGrid(std::array<int, 2> nCells,
                             std::array<double, 2> meshWidth)
    : u_({nCells[0] + 1, nCells[1] + 2},
         {0.0 + meshWidth[0], 3 * meshWidth[1] / 2.}, meshWidth, "u"),
      v_({nCells[0] + 2, nCells[1] + 1}, {3 * meshWidth[0] / 2., meshWidth[1]},
         meshWidth, "v"),
      p_({nCells[0] + 2, nCells[1] + 2},
         {3 * meshWidth[0] / 2., 3 * meshWidth[1] / 2.}, meshWidth, "p"),
      f_({nCells[0] + 1, nCells[1] + 2}, {meshWidth[0], 3 * meshWidth[1] / 2.},
         meshWidth, "f"),
      g_({nCells[0] + 2, nCells[1] + 1}, {3 * meshWidth[0] / 2., meshWidth[1]},
         meshWidth, "g"),
      rhs_({nCells[0], nCells[1]},
           {3 * meshWidth[0] / 2., 3 * meshWidth[1] / 2.}, meshWidth, "rhs"),
      meshWidth_(meshWidth), nCells_(nCells),
      ui0_(1), uj0_(1), vi0_(1), vj0_(1), pi0_(1), pj0_(1),
      uGhost_(0), vGhost_(0) {}

StaggeredGrid::StaggeredGrid(PartitionInformation &pi)
    : u_({pi.nCellsLocal()[0] + 1 + pi.uGhostLayer(), pi.nCellsLocal()[1] + 2},
         {0.0 + pi.meshWidth()[0], 3 * pi.meshWidth()[1] / 2.}, pi.meshWidth(), "u"),
      v_({pi.nCellsLocal()[0] + 2, pi.nCellsLocal()[1] + 1 + pi.vGhostLayer()}, 
         {3 * pi.meshWidth()[0] / 2., pi.meshWidth()[1]}, pi.meshWidth(), "v"),
      p_({pi.nCellsLocal()[0] + 2, pi.nCellsLocal()[1] + 2},
         {3 * pi.meshWidth()[0] / 2., 3 * pi.meshWidth()[1] / 2.}, pi.meshWidth(), "p"),
      f_({pi.nCellsLocal()[0] + 1 + pi.uGhostLayer(), pi.nCellsLocal()[1] + 2},
         {pi.meshWidth()[0], 3 * pi.meshWidth()[1] / 2.}, pi.meshWidth(), "f"),
      g_({pi.nCellsLocal()[0] + 2, pi.nCellsLocal()[1] + 1 + pi.vGhostLayer()},
         {3 * pi.meshWidth()[0] / 2., pi.meshWidth()[1]}, pi.meshWidth(), "g"),
      rhs_(pi.nCellsLocal(),
           {3 * pi.meshWidth()[0] / 2., 3 * pi.meshWidth()[1] / 2.}, pi.meshWidth(), "rhs"),
      meshWidth_(pi.meshWidth()), nCells_(pi.nCellsLocal()),
      ui0_(1), uj0_(1), vi0_(1), vj0_(1), pi0_(1), pj0_(1),
      uGhost_(pi.uGhostLayer()), vGhost_(pi.vGhostLayer()) {}

void StaggeredGrid::makeCGFields()
{
   //! this function may only be called once
   if(cgFieldsMade_)
      throw std::runtime_error("makeCGFields may only be called once\n");

   cgFieldsMade_ = true;
   // the origin probably doesn't matter
   r_  = std::make_shared<FieldVariable>(p_.size(),   p_.getOrigin(), meshWidth_, "r");
   a_  = std::make_shared<FieldVariable>(p_.size(),   p_.getOrigin(), meshWidth_, "a");
}