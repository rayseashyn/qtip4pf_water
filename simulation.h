#ifndef SIMULATION
#define SIMULATION

void run_classical_wspline(double position_array[], std::ofstream& outfile, std::mt19937& rng);
void run_classical_nospline(double position_array[], std::ofstream& outfile, std::mt19937& rng);
void run_sc_wspline(double position_array[], std::ofstream& outfile, std::mt19937& rng);
void run_sc_nospline(double position_array[], std::ofstream& outfile, std::mt19937& rng);

#endif
