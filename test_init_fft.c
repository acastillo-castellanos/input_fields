#include "view.h"
#include "initial_conditions_dimonte_fft2.h"

#define r2 (sq(x) + sq(y))
int main()
{
  L0 = 2*pi;
  X0 = Y0 = Z0 = -L0 / 2;
  N = 1 << 9;
  init_grid(N);

#if TREE
  refine(((r2 < sq(0.25)) && (r2 > sq(0.1))) && level < 9);
#endif

  scalar f[];
  initial_condition_dimonte_fft2(f, 1, kmin = 25, kmax = 75);

  stats s = statsf (f);
  fprintf (stderr, "# %f %.12f %g %g %g %g\n", t, s.sum, s.min, s.max, s.stddev, s.volume);

  box();
  cells();
  save("test1.png");

  squares("f", linear = false);
  save("test2.png");
}