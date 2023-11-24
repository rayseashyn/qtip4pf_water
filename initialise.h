#ifndef WATER_INIT
#define WATER_INIT
#include <array>

double* initialise_fcc();
double* initialise_fromfile(std::string filename);
std::array<double, 3> get_dipole_moment(const double position_array[]);
double* initialise_oxygens();
std::vector<std::array<int, 2> > get_oxygen_connectivity();
double* initialise_orth(std::mt19937& rng, std::ofstream& outfile);
double* initialise_positions(std::string cell_type, std::mt19937& rng, std::ofstream& outfile);

#endif
