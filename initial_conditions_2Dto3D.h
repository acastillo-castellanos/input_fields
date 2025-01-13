#include "auxiliar_input.h"
#define cond1(y,d1,d2) ((l[] < level) && ((y < d1) && (y > d2)))
#define cond2(y,d1,d2,del) ((y > d1+del) || (y < d2-del))
#define wallbox(d, extra) intersection((d + extra - y), (-d + extra + y))
void initial_condition_2Dto3D(scalar f, vector u, scalar p, double d1, double d2){
  // Initialize transveral velocity 
	foreach()
	  u.y[] = 0.;
  
  // Now, we match the refinement level
  scalar l[];
  read_matrix(file_restart_path, "_l", l);      

#if dimension == 3
  int maxlevel = MAXLEVEL;
  for (int li = maxlevel; li >= 4; li--){
    unrefine( (cond1(y,d1,d2) || cond2(y,d1,d2,16*_mindel)) && level > li);
  }
#endif

  // Then, we read the corresponding fields
  read_matrix(file_restart_path, "_f", f);
  read_matrix(file_restart_path, "_u", u.x);
  #if dimension == 2
    read_matrix(file_restart_path, "_v", u.y);
  #else
    read_matrix(file_restart_path, "_v", u.z);
  #endif
  read_matrix(file_restart_path, "_p", p);
}
#undef wallbox
#undef cond1
#undef cond2