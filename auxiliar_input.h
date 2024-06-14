/*
 * input_matrix() reads a gnuplot-compatible binary file, i.e., single precision 
 * matrix stored in the following format:
 *
 *  <N+1>  <y0>   <y1>   <y2>  ...  <yN>
 *  <x0> <z0,0> <z0,1> <z0,2> ... <z0,N>
 *  <x1> <z1,0> <z1,1> <z1,2> ... <z1,N>
 *
 * For example: input_matrix(T, fp, N, X0, Y0, L0);
 * Reads a square matrix of size N from file descriptor "fp" into scalar field 
 * T, defined at points (X,Y).
 */

// Define a structure to hold input matrix parameters
struct InputMatrix {
    scalar s;  // Scalar field to store the matrix values
    FILE * fp; // File pointer for the binary input file
    int n;     // Size of the matrix (optional)
    double ox, oy, width, oz; // Offsets and dimensions (optional)
    double arz; // Additional parameter (optional)
};

// Function to read and store the matrix from the binary file
void input_matrix (struct InputMatrix p) {
	scalar s = p.s; // Assign the scalar field
	if (p.width == 0.0) p.width = L0; // Use default width if not provided

	float width=0;
	fread(&width, sizeof(float), 1, p.fp); // Read the matrix width from file

	if (p.n != width) p.n = (int) width; // Set matrix size

	float yp[p.n], xp[p.n]; // Arrays to store y and x coordinates
	float ** v = matrix_new (p.n, p.n, sizeof(float));  // Allocate memory for matrix
	
	fread(&yp, sizeof(float), p.n, p.fp); // Read y coordinates from file
	for (int i = 0; i < p.n; i++) {
		fread(&xp[i], sizeof(float), 1, p.fp); // Read x coordinate for each row
		for (int j = 0; j < p.n; j++) {
			fread(&v[i][j], sizeof(float), 1, p.fp); // Read matrix values
		}
	}
	
	// Loop over the domain and assign matrix values to the scalar field
	foreach() {
		int i = (x - p.ox) * width / p.width;
        int j = (y - p.oy) * width / p.width;
        if (i >= 0 && i < width && j >= 0 && j < width) {
            s[] = v[i][j];
        } else {
            s[] = 0.0;
        }
	}
	
	// Free the allocated memory
	matrix_free(v);
}

void input_matrix_double (struct InputMatrix p) {
	scalar s = p.s; // Assign the scalar field
	if (p.width == 0.0) p.width = L0; // Use default width if not provided

	double width=0;
	fread(&width, sizeof(double), 1, p.fp); // Read the matrix width from file

	if (p.n != width) p.n = (int) width; // Set matrix size

	double yp[p.n], xp[p.n]; // Arrays to store y and x coordinates
	double ** v = matrix_new (p.n, p.n, sizeof(double));  // Allocate memory for matrix
	
	fread(&yp, sizeof(double), p.n, p.fp); // Read y coordinates from file
	for (int i = 0; i < p.n; i++) {
		fread(&xp[i], sizeof(double), 1, p.fp); // Read x coordinate for each row
		for (int j = 0; j < p.n; j++) {
			fread(&v[i][j], sizeof(double), 1, p.fp); // Read matrix values
		}
	}
	
	// Loop over the domain and assign matrix values to the scalar field
	foreach() {
		int i = (x - p.ox) * width / p.width;
        int j = (y - p.oy) * width / p.width;
        if (i >= 0 && i < width && j >= 0 && j < width) {
            s[] = v[i][j];
        } else {
            s[] = 0.0;
        }
	}
	
	// Free the allocated memory
	matrix_free(v);
}

void read_matrix(const char *prefix, const char *suffix, scalar s) {
    char filename[256];
    snprintf(filename, sizeof(filename), "%s%s.bin", prefix, suffix);
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        printf("Binary file %s not found\n", filename);
        exit(EXIT_FAILURE);
    }
    input_matrix_double(s, fp, N, X0, Y0, L0);
    fclose(fp);
}

