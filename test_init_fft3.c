#include "grid/octree.h"
#include "view.h"
#include "initial_conditions_dimonte_fft2.h"

#define r2 (sq(x) + sq(y))
int main()
{
  L0 = 2*pi;
  X0 = Y0 = Z0 = -L0 / 2;
  N = 1 << 5;
  init_grid(N);

#if TREE
  refine(((z < L0/16) && (z > -L0/16)) && level < 9);
#endif

  scalar f[];
  initial_condition_dimonte_fft2(f, amplitude=0.01, NX=512, NY=512, kmin = 10, kmax = 50, isvof=1);

  stats s = statsf (f);
  fprintf (stderr, "# %f %.12f %g %g %g %g\n", t, s.sum, s.min, s.max, s.stddev, s.volume);

  view(camera="top");
  draw_vof ("f", color="z");
  save("test1.png");

  view(camera="front");
  draw_vof ("f", color="z");
  save("test2.png");

  view(camera="iso");
  draw_vof ("f", color="z");
  save("test3.png");

  squares("f", linear = false);
  save("test4.png");
}