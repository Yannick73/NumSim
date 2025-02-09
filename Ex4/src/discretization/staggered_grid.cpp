#include "discretization/staggered_grid.h"

enum field : char
{
   u = 'u',
   v = 'v',
   w = 'w',
   p = 'p',
   rhs = 'r'
};

inline std::array<int, 3> getCells(PartitionInformation &pi, field ind)
{
   switch (ind)
   {
   case field::u:
      return {pi.nCellsLocal()[0] + 1 + pi.uGhostLayer(),
              pi.nCellsLocal()[1] + 2,
              pi.nCellsLocal()[2] + 2};
   case field::v:
      return {pi.nCellsLocal()[0] + 2,
              pi.nCellsLocal()[1] + 1 + pi.vGhostLayer(),
              pi.nCellsLocal()[2] + 2};
   case field::w:
      return {pi.nCellsLocal()[0] + 2,
              pi.nCellsLocal()[1] + 2,
              pi.nCellsLocal()[2] + 1 + pi.wGhostLayer()};
   case field::p:
      return {pi.nCellsLocal()[0] + 2,
              pi.nCellsLocal()[1] + 2,
              pi.nCellsLocal()[2] + 2};
   case field::rhs:
      return {pi.nCellsLocal()[0],
              pi.nCellsLocal()[1],
              pi.nCellsLocal()[2]};
   default:
      throw std::range_error("Undefined index");
   }
}

inline std::array<double, 3> getOrigin(PartitionInformation &pi, field ind)
{
   switch (ind)
   {
   case field::u:
      return {0.0,
              -pi.meshWidth()[1] / 2.,
              -pi.meshWidth()[2] / 2.};
   case field::v:
      return {-pi.meshWidth()[0] / 2.,
              0.0,
              -pi.meshWidth()[2] / 2.};
   case field::w:
      return {-pi.meshWidth()[0] / 2.,
              -pi.meshWidth()[1] / 2.,
              0.0};
   case field::p:
      return {-pi.meshWidth()[0] / 2.,
              -pi.meshWidth()[1] / 2.,
              -pi.meshWidth()[2] / 2.};
   case field::rhs:
      return {pi.meshWidth()[0] / 2.,
              pi.meshWidth()[1] / 2.,
              pi.meshWidth()[2] / 2.};
   default:
      throw std::range_error("Undefined index");
   }
}

StaggeredGrid::StaggeredGrid(PartitionInformation &pi)
    :   u_(getCells(pi, field::u),   getOrigin(pi, field::u),   pi.meshWidth(), "u"),
        v_(getCells(pi, field::v),   getOrigin(pi, field::v),   pi.meshWidth(), "v"),
        w_(getCells(pi, field::w),   getOrigin(pi, field::w),   pi.meshWidth(), "w"),
        f_(getCells(pi, field::u),   getOrigin(pi, field::u),   pi.meshWidth(), "f"),
        g_(getCells(pi, field::v),   getOrigin(pi, field::v),   pi.meshWidth(), "g"),
        h_(getCells(pi, field::w),   getOrigin(pi, field::w),   pi.meshWidth(), "h"),
        p_(getCells(pi, field::p),   getOrigin(pi, field::p),   pi.meshWidth(), "p"),
      rhs_(getCells(pi, field::rhs), getOrigin(pi, field::rhs), pi.meshWidth(), "rhs"),
      meshWidth_(pi.meshWidth()), nCells_(pi.nCellsLocal()),
      ui0_(1), uj0_(1), uk0_(1),
      vi0_(1), vj0_(1), vk0_(1),
      wi0_(1), wj0_(1), wk0_(1),
      pi0_(1), pj0_(1), pk0_(1),
      uGhost_(pi.uGhostLayer()), vGhost_(pi.vGhostLayer()), wGhost_(pi.wGhostLayer()) {}