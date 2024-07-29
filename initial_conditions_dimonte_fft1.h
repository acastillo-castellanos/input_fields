/**
# Initializing an interface using a given spectrum

To initialize the interface we'll use a Fourier transform and some interpolation
routines from GSL.
*/

#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_spline.h>
#pragma autolink -lgsl -lgslcblas
#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])

/**
## Save the initial perturbation in a gnuplot-compatible format 
### save_data_for_gnuplot_complex(): saves a (complex) 1D array as a .dat file
*/

void save_data_for_gnuplot_complex(double *data, int NX, const char *filename){
  FILE *file = fopen(filename, "w");
  if (file == NULL){
    fprintf(stderr, "Error opening file %s\n", filename);
    return;
  }

  for (int i = 0; i < NX; i++){
    double magnitude = sqrt(sq(REAL(data, i)) + sq(IMAG(data, i)));
    fprintf(file, "%d %f %f %f\n", i, REAL(data,i), IMAG(data,i), magnitude);
  }
  fclose(file);
}

/** 
### save_data_for_gnuplot_real(): saves a 1D array as a .dat file
*/
void save_data_for_gnuplot_real(double *data, int NX, const char *filename){
  FILE *file = fopen(filename, "w");
  if (file == NULL){
    fprintf(stderr, "Error opening file %s\n", filename);
    return;
  }

  for (int i = 0; i < NX; i++){
    fprintf(file, "%d %f \n", i, data[i]);
  }
  fclose(file);
}

/** 
### init_1D_complex(): initializes the perturbation in Fourier space
*/

void init_1D_complex(double *data, int N, double kmin, double kmax, double eta0_target=1){
  double *kx = malloc(N * sizeof(double));
  double cst = eta0_target / sqrt((2./kmin)-(2./kmax));

  // Calculate wavenumbers
  for (int i = 0; i <= N / 2; ++i)
    kx[i] = 2 * pi * i / L0;

  for (int i = N / 2 + 1; i < N; ++i)
    kx[i] = 2 * pi * (i - N) / L0;

  // Initialize spectrum in the specified range with magnitude 1/k and random phase
  double dkx = kx[1]-kx[0];
  double eta0 = 0.;
  memset(data, 0, 2 * N * sizeof(double));
  for (int i = 0; i < N; ++i){
    double k = sqrt(sq(kx[i]));
    if ((k >= kmin) && (k < kmax)){
      double magnitude = cst / k;
      double phase = noise() * pi;
      REAL(data, i) = magnitude * cos(phase);
      IMAG(data, i) = magnitude * sin(phase);
      eta0 += sq(magnitude)*dkx;
    }    
  }
  fprintf(stdout, "real eta0 is %g \n", sqrt(eta0));

  free(kx);
}


/** 
### initial_condition_dimonte_fft(): initializes a scalar field f 
*/
void initial_condition_dimonte_fft(vertex scalar phi, double amplitude=1, int NX=N, double kmin=1, double kmax=1){
  
  /** We declare the arrays and initialize the physical space*/
  double *data = malloc(2 * NX * sizeof(double));
  double *xdata = (double *)malloc(NX * sizeof(double));
  double *ydata = (double *)malloc(NX * sizeof(double));

  double dx = L0 / (NX - 2);
  for (int i = 0; i < NX; i++){
    xdata[i] = i * dx + X0;
  }

  /** The perturbation is generated only by the main process */
  if (pid() == 0){
    // Initialize the spectrum
    init_1D_complex(data, NX, kmin, kmax, eta0_target=amplitude);
    save_data_for_gnuplot_complex(data, NX, "initial_spectra.dat");

    // Perform the FFT2D 
    gsl_fft_complex_radix2_backward(data, 1, NX);
    save_data_for_gnuplot_complex(data, NX, "final_deformation.dat");

    // Save the results into a 2D array
    for (int i = 0; i < NX; i++){
      ydata[i] = REAL(data,i);
    }
    save_data_for_gnuplot_real(ydata, NX, "final_deformation2.dat");
  }

  /** and broadcasted to the other processes if MPI */
  @ if _MPI
    MPI_Bcast(ydata, NX, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  @endif

  /** Now, we'll interpolate the perturbation into the mesh. First, we set up
  the accelerations */
  gsl_interp *interp = gsl_interp_alloc(gsl_interp_linear, NX);
  gsl_interp_accel *acc = gsl_interp_accel_alloc();

  /** Initialize the interpolation object */ 
  gsl_interp_init(interp, xdata, ydata, NX);

  /** Apply initial condition to the scalar field. Here, we take care to
  normalize the perturbation using the standard deviation and multiply it by
  some amplitude. Also, we may apply the perturbation directly (`isvof=0`) to a
  field or use it to generate a VOF surface (`isvof=1`)
  */

  
  foreach_vertex()
    phi[] = gsl_interp_eval(interp, xdata, ydata, x, acc) - y;
  

  /** Release interpolation objects */ 
  gsl_interp_free(interp);
  gsl_interp_accel_free(acc);

  /** Free Dynamically Allocated Memory */ 
  free(xdata);
  free(ydata);
  free(data);
}
