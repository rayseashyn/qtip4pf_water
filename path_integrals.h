#ifndef PATH_INTEGRALS
#define PATH_INTEGRALS

double* get_q_centroid(const double input_array, const int length);
double get_kcv(const double position_array[], const double force_array[]);
double get_conserved(double position_array[], double momentum_array[], const double mass_array[], double energy);

#endif
