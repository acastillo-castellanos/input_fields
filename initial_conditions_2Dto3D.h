#include "auxiliar_input.h"
#define wallbox(d, extra) intersection((d + extra - y), (-d + extra + y))
void initial_condition_2Dto3D(scalar f, vector u, scalar p, double d1, double d2){
  // Initialize transveral velocity 
	foreach()
	  u.y[] = 0.;
  
  // Now, we match the refinement level
  scalar l[];
  for (int i = 0; i < npe(); i++)
    if (pid() == i)
      read_matrix(file_restart_path, "_l", l);

  int maxlevel = MAXLEVEL;
  for (int li = maxlevel; li >= 4; li--){
    unrefine((l[] < level && ((y < d1) && (y > d2))) && level > li);
  }

  // Then, we read the corresponding fields
  for (int i = 0; i < npe(); i++){
    if (pid() == i){
      read_matrix(file_restart_path, "_f", f);
      read_matrix(file_restart_path, "_u", u.x);
      #if dimension == 2
        read_matrix(file_restart_path, "_v", u.y);
      #else
        read_matrix(file_restart_path, "_v", u.z);
      #endif
      read_matrix(file_restart_path, "_p", p);
    }
  }
}
#undef wallbox
