/**
# input_matrix(): reads a gnuplot-compatible binary file
These routines read a gnuplot-compatible binary file, i.e., single precision
matrix stored in the following format:

```
 <N+1>  <y0>   <y1>   <y2>  ...  <yN>
 <x0> <z0,0> <z0,1> <z0,2> ... <z0,N>
 <x1> <z1,0> <z1,1> <z1,2> ... <z1,N>
```

For example: using `input_matrix(T, fp, N, X0, Y0, L0)` will read a square
matrix of size `N` from file descriptor `fp` into a scalar field `T`.

## input_matrix(): reads and store the matrix from the binary file */
void input_matrix(scalar s, FILE * fp, int n = N, double ox=X0, double oy=Y0, double width=L0){

  float nfile;
  fread(&nfile, sizeof(float), 1, fp);          // Read the matrix size from file
  n = (int)nfile;                               // Set matrix size

  float yp[n], xp[n];                           // Arrays to store y and x coordinates
  float **v = matrix_new(n, n, sizeof(float));  // Allocate memory for matrix

  fread(&yp, sizeof(float), n, fp);             // Read y coordinates from file
  for (int i = 0; i < n; i++){
    fread(&xp[i], sizeof(float), 1, fp);        // Read x coordinate for each row
    for (int j = 0; j < n; j++){
      fread(&v[i][j], sizeof(float), 1, fp);    // Read matrix values
    }
  }

  /** Loop over the domain and assign matrix values to the scalar field */ 
  foreach () {
    int i = (x - ox) * n / width;
    int j = (y - oy) * n / width;
    if (i >= 0 && i < n && j >= 0 && j < n){
      s[] = v[i][j];
    }
    else{
      s[] = 0.0;
    }
  }
  matrix_free(v);
}

/** 
## input_matrix_double(): reads and store the matrix from a double precision binary file 
*/
void input_matrix_double(scalar s, FILE * fp, int n = N, double ox=X0, double oy=Y0, double width=L0){

  double nfile;
  fread(&nfile, sizeof(double), 1, fp);           // Read the matrix size from file
  n = (int)nfile;                                 // Set matrix size

  double yp[n], xp[n];                            // Arrays to store y and x coordinates
  double **v = matrix_new(n, n, sizeof(double));  // Allocate memory for matrix

  fread(&yp, sizeof(double), n, fp);              // Read y coordinates from file
  for (int i = 0; i < n; i++){
    fread(&xp[i], sizeof(double), 1, fp);         // Read x coordinate for each row
    for (int j = 0; j < n; j++){
      fread(&v[i][j], sizeof(double), 1, fp);     // Read matrix values
    }
  }

  /** Loop over the domain and assign matrix values to the scalar field */ 
  foreach () {
    int i = (x - ox) * n / width;
    int j = (y - oy) * n / width;
    if (i >= 0 && i < n && j >= 0 && j < n){
      s[] = v[i][j];
    }
    else{
      s[] = 0.0;
    }
  }
  matrix_free(v);
}

/** 
## read_matrix(): reads and store the matrix from a double precision binary file 
*/
void read_matrix(const char *prefix, const char *suffix, scalar s) {
  char filename[256];
  snprintf(filename, sizeof(filename), "%s%s.bin", prefix, suffix);
  FILE *fp = fopen(filename, "r");
  if (!fp){
    printf("Binary file %s not found\n", filename);
    exit(EXIT_FAILURE);
  }
  input_matrix_double(s, fp, N, X0, Y0, L0);
  fclose(fp);
}
