//definition of various parameters needed from the qtip4pf potential, as well as system-level
//parameters, such as the temperation and the simulation size.
//some of these will be defined in the input file - include some of the intramolecular
//qtip4pf parameters as these may be changed by the procedure we used to generate potentials
//for the centroid molecular dynamics simulations

#ifndef QTIP4PFPARAMS
#define QTIP4PFPARAMS
#include <math.h>
#include <cmath>
#include <array>

//------------------------------------------
//input parameters for you to change
//------------------------------------------
//simulation temperature in Kelvin
extern double T_K;
//number of beads used in the simulation (should be even for Suzuki-Chin)
extern int n_beads;
//parameter for Suzuki-Chin finite difference (if you want to play around with this)
extern double epsilonsc;
//number of starting fcc unit cells in each direction
//this is the smallest cell that is larger than 2.0*r_cut_LJ in all directions for the current density
extern int n_cell_x;
extern int n_cell_y;
extern int n_cell_z;
//density in g/cm^3
extern double rho_inp;
//simulation timestep in fs
extern double dtfs;
//approximate number of significant figures given accurately by ewald sum
extern double eps;
//number of equilibration steps
extern int n_equil;
//number of sampling steps
extern int n_samp;
const int n_step = 100;
extern bool is_spline;

//cutoff radius for LJ interaction
const double r_cut_angs = 9.0;

//numerical constants it's useful to have
//and unit conversions to atomic units, which this code uses throughout
const double pi = 3.14159265359;
const double angs2bohr = 1.0 / 0.529177210903;
const double me = 9.10938356e-31;
const double kb = 1.38064852e-23;
const double Na = 6.0221409e23;
const double eh = 4.3597447222071e-18;
const double kcal2au = 1.5936e-3;
const double mhgmol = 1.0080;
const double mogmol = 15.9994;
//masses of oxygen and hydrogen in atomic units
const double mo = 0.001 * mogmol / (Na * me);
const double mh = 0.001 * mhgmol / (Na * me);
const double K2au = kb / eh;
const double hbar_SI = 1.054571817e-34;
const double hbar = 1.0;
const double au2fs = hbar_SI * 1e15 / eh;
const double fs2au = 1.0 / au2fs;

//store basic fractions and sqrts for ease of use in generating starting lattice if not reading from file
const double one_third = 1.0 / 3.0;
const double two_thirds = 2.0 * one_third;
const double one_sixth = 0.5 * one_third;
const double five_sixths = 5.0 / 6.0;
const double sqrt3 = sqrt(3.0);
//and in atomic units
extern double T_red;
//inverse temperature
extern double beta;
//various quantities needed for path integrals
extern double piovern;
extern double beta_n;
extern double omega_p;
extern double omega_psq;

//and total number of such unit cells
extern int n_cells;
//number of molecules in the simulation
//and related stuff like number of atoms and degrees of freedom in each replica
//four molecules in each fcc cell
extern int n_mols;
extern int n_atoms;
extern int n_dof;
//dimensions of the orthorhombic unit cell in atomic units
const double mh2o = (mhgmol*2.0+mogmol) / Na;
extern double a_cell_angs;
extern double a_cell;
extern double b_cell;
extern double c_cell;
//parameter used in starting oxygen positions if not reading from file
//1/16 is the "idealised" number - in real ice I believe it is close to but not exactly this
const double zo = 1.0 / 16.0;
extern double rnn;
//dimensions of the total simulation box stored as a 3D array
//for ease of access in implementing periodic boundary conditions
extern std::array<double, 3> boxl;
//volume of each fcc unit cell
extern double v_cell;
//and total volume of the simulation box
extern double vol;
//cutoff radius for LJ oxygen interactions
//this could be a little short - have seen 10.5 \AA used in some papers
const double r_cut = r_cut_angs * angs2bohr;
const double r_cut_sq = r_cut * r_cut;
//"skin" radius for Verlet neighbour list
const double r_skin = 1.30 * r_cut;
const double r_skin_sq = r_skin * r_skin;
//density of oxygen atoms in atomic units
extern double rho_ox;
//approximate size of neighbour list
extern int list_size;

//various parameters for defining the qTIP4P/F potential - the meaning of each should be self-explanatory from David's original paper
const double gamma_vs = 0.73612;
const double epslj = 0.1852 * kcal2au;
const double sigmalj = 3.1589 * angs2bohr;
const double sigmalj6 = std::pow(sigmalj, 6.0);
const double sigmalj12 = sigmalj6 * sigmalj6;
//energy of LJ interaction at the chosen cutoff radius
const double e_cut = 4.0 * epslj * (std::pow((sigmalj/r_cut), 12.0) - std::pow((sigmalj/r_cut), 6.0));
const double roheq = 0.9419 * angs2bohr;
const double ar = 2.287 / angs2bohr;
const double ar2 = ar * ar;
const double ar3 = ar2 * ar;
const double ar4fac = 7.0 * ar3 * ar / 12.0;
const double Dr = 116.09 * kcal2au;
const double theta_eq = 107.4 * pi / 180.0;
const double ktheta = 87.85 * kcal2au;
const double qo = -1.1128;
const double qh = -0.5*qo;

//standard deviations for Boltzmann distribution of momenta when initialising for oxygen and hydrogen atoms
extern double sd_ox;
extern double sd_h;


//simulation parameters - eg temperature, timestep etc - go here
//-----------------------------------------------------
extern double dt;
extern double halfdt;

//calculation of rwald and kwald cutoff radii
//-----------------------------------------------------
extern double acc;
extern double rho;
//we need to find the minimum box side length out of x, y and z
//the only way I can get this to work on Europa is rather ugly - find the minimum out of x and y and then check that against z in the next line
//checking all three in one line works fine on my laptop though, just not on Europa
extern double min_xy;
extern double min_box;
extern double kwald_prefactor;
extern double max_xy;
extern double max_box;

#endif
