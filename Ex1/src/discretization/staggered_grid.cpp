#include "discretization/staggered_grid.h"

StaggeredGrid::StaggeredGrid(std::array<int, 2> nCells,
                             std::array<double, 2> meshWidth)
    : u_({nCells[0] + 4, nCells[1] + 4},
         {0.0 + meshWidth[0], 3 * meshWidth[1] / 2.}, meshWidth),
      v_({nCells[0] + 4, nCells[1] + 4}, {3 * meshWidth[0] / 2., meshWidth[1]},
         meshWidth),
      p_({nCells[0] + 4, nCells[1] + 4},
         {3 * meshWidth[0] / 2., 3 * meshWidth[1] / 2.}, meshWidth),
      f_({nCells[0] + 4, nCells[1] + 4}, {meshWidth[0], 3 * meshWidth[1] / 2.},
         meshWidth),
      g_({nCells[0] + 4, nCells[1] + 4}, {3 * meshWidth[0] / 2., meshWidth[1]},
         meshWidth),
      rhs_({nCells[0] + 4, nCells[1] + 4},
           {3 * meshWidth[0] / 2., 3 * meshWidth[1] / 2.}, meshWidth),
      meshWidth_(meshWidth), nCells_(nCells) {}

void StaggeredGrid::setBoundary(Settings settings) 
{
   dirichletBcBottom = settings.dirichletBcBottom;
   dirichletBcTop = settings.dirichletBcTop;
   dirichletBcLeft = settings.dirichletBcLeft;
   dirichletBcRight = settings.dirichletBcRight;
}

void StaggeredGrid::setBoundaryValues()
{
    const int nx = nCells()[0];
    const int ny = nCells()[1];
   // left and right boundaries take priority, +1 for corner
   for(int j = -1; j < ujEnd()+1; j++) 
   {
      // left edge x-direction
      u(-1,j) = dirichletBcLeft[0];
      f(-1,j) = u(-1,j);
      // right edge x-direction, field is one cell smaller
      u(nx-1,j) = dirichletBcRight[0];
      f(nx-1,j) = u(nx-1,j);
   }

   for(int j = -1; j < vjEnd()+1; j++)
   {
      // left edge y-direction
      v(-1,j) = 2*dirichletBcLeft[1]-v(0,j);
      g(-1,j) = v(-1,j);
      // right edge y-direction
      v(nx,j) = 2*dirichletBcRight[1]-v(nx-1,j);
      g(nx,j) = v(nx,j);
   }
   
   for(int i=0; i < uiEnd(); i++) 
   {
      // lower edge x-direction
      u(i,-1) = 2*dirichletBcBottom[0]-u(i, 0);
      f(i,-1) = u(i, -1);
      // upper edge x-direction
      u(i,ny) = 2*dirichletBcTop[0]-u(i, ny-1);
      f(i,ny) = u(i, ny);
   }
   
   for(int i=0; i < viEnd(); i++) 
   {
      // lower edge y-direction
      v(i,-1) = dirichletBcBottom[1];
      g(i,-1) = v(i, -1);
      // upper edge y-direction, field is one cell smaller
      v(i,ny-1) = dirichletBcTop[1];
      g(i,ny-1) = v(i, ny-1);
   }
}
