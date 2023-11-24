//main file, reads in the desired input data from file and then runs the necessary simulation

#include <stdio.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <math.h>
#include <array>
#include <string>
#include <sstream>
#include <fftw3.h>
#include "omp.h"
#include <boost/math/interpolators/cubic_b_spline.hpp>
#include "initialise.h"
#include "qtip4pfparams.h"
#include "neighbour.h"
#include "forces.h"
#include "thermostat.h"
#include "update.h"
#include "indexing.h"
#include "path_integrals.h"
#include "simulation.h"

//define parameters to be read in from config file and derived parameters in file scope so that it can be used in the qtip4pfparams header (and so accessed by all other files that include that header)
double T_K, epsilonsc, rho_inp, dtfs, eps, a_cell_angs, a_cell, b_cell, c_cell, rnn;
int n_beads, n_cell_x, n_cell_y, n_cell_z, n_cells;
int n_mols, n_atoms, n_dof, list_size;
double T_red, beta, piovern, beta_n, omega_p, omega_psq;
double v_cell, vol, rho_ox, sd_ox, sd_h, dt, halfdt;
double acc, rho, min_xy, min_box, kwald_prefactor, max_xy, max_box;
std::array<double, 3> boxl;
int n_equil, n_samp;
bool is_spline;
std::string cell_type, from_file, filename;

//main function for setting up and running the simulation
int main () {
    
    //read in parameters from config file
    std::ifstream config;
    std::string string_name;
    config.open("config.txt");
    while (config.is_open()) {
        config >> string_name >> n_beads;
        config >> string_name >> T_K;
        config >> string_name >> epsilonsc;
        config >> string_name >> n_cell_x;
        config >> string_name >> n_cell_y;
        config >> string_name >> n_cell_z;
        config >> string_name >> rho_inp;
        config >> string_name >> dtfs;
        config >> string_name >> eps;
        config >> string_name >> n_equil;
        config >> string_name >> n_samp;
        config >> string_name >> is_spline;
        config >> string_name >> cell_type;
        config >> string_name >> from_file;
        if (from_file=="yes") {
            config >> filename;
        }
        config.close();
    }
    
    //check if n_beads is a bad number or not
    if (n_beads != 1 && n_beads%2 == 1) {
        std::cerr << "n_beads should be even for Suzuki-Chin calculations, or set to 1 for classical calculations \n";
        exit(1);
    }
    
    //find derived parameters from input parameters
    //reduced temperature in atomic units
    T_red = K2au * T_K;
    //inverse temperature
    beta = 1.0 / T_red;
    //various quantities needed for path integrals
    piovern = pi / n_beads;
    beta_n = beta / n_beads;
    omega_p = 1.0 / (beta_n*hbar);
    omega_psq = omega_p * omega_p;
    
    //total number of simulation cells
    n_cells = n_cell_x * n_cell_y * n_cell_z;
    
    //number of molecules in the simulation
    //and related stuff like number of atoms and degrees of freedom in each replica
    if (cell_type == "fcc") {
        //four molecules in each fcc cell
        n_mols = 4 * n_cells;
        //unit cell parameter
        a_cell_angs = 1e8 * std::cbrt(4 * mh2o / rho_inp);
        a_cell = a_cell_angs * angs2bohr;
        b_cell = a_cell;
        c_cell = a_cell;
    }
    
    else if (cell_type == "orth") {
        //eight molecules in each orthorhombic unit cell
        n_mols = 8 * n_cells;
        
        a_cell_angs = 1e8 * std::cbrt(4 * mh2o / (sqrt(2.0)*rho_inp));
        a_cell = a_cell_angs * angs2bohr;
        b_cell = sqrt3 * a_cell;
        c_cell = 2.0 * sqrt(2.0/3.0) * a_cell;
        
        rnn = (0.5 - 2.0 * zo) * c_cell;
    }
    
    else {
        std::cerr << "cell type not supported, should be \"orth\" or \"fcc\" \n";
        exit(1);
    }
    n_atoms = 3 * n_mols;
    n_dof = 3 * n_atoms;
    
    //box length stored as an array
    boxl = {n_cell_x * a_cell, n_cell_y * b_cell, n_cell_z * c_cell};
    
    //volume of each fcc unit cell
    v_cell = a_cell * b_cell * c_cell;
    //and total volume of the simulation box
    vol = n_cells * v_cell;
    
    //density of oxygen atoms in atomic units
    if (cell_type == "fcc") {
        rho_ox = 4.0 / v_cell;
    }
    else if (cell_type == "orth") {
        rho_ox = 8.0 / v_cell;
    }
    //approximate size of neighbour list
    list_size = std::round(4.0 * r_skin * r_skin * r_skin * n_mols * rho_ox);
    
    //standard deviations for Boltzmann distribution of momenta when initialising for oxygen and hydrogen atoms
    sd_ox = sqrt(mo/beta_n);
    sd_h = sqrt(mh/beta_n);
    
    //simulation parameters - eg temperature, timestep etc - go here
    //-----------------------------------------------------
    dt = dtfs * fs2au;
    halfdt = 0.5 * dt;

    //calculation of rwald and kwald cutoff radii
    //-----------------------------------------------------
    acc = sqrt(-log(eps));
    rho = n_atoms / vol;
    //we need to find the minimum box side length out of x, y and z
    //the only way I can get this to work on Europa is rather ugly - find the minimum out of x and y and then check that against z in the next line
    //checking all three in one line works fine on my laptop though, just not on Europa
    min_xy = std::min(boxl[0], boxl[1]);
    min_box = std::min(min_xy, boxl[2]);
    kwald_prefactor = 2.0 * pi / vol;
    max_xy = std::max(boxl[0], boxl[1]);
    max_box = std::max(max_xy, boxl[2]);

    std::random_device rd;
    //save seed so can reproduce run again later if needed
    //the implementation of either mt19937 or gamma_distribution/normal_distribution appears to be different on the c++ installations on my machine and dirac/europa
    //but regardless, using the same seed on europa reproduces the same trajectory on europa
    //and the same for dirac and my machine
    //so if you want to reproduce trajectories, they need to be run on the same cluster/machine
    unsigned int seed = rd();
    std::ofstream outfile("output.txt");
    if (outfile.is_open()) {
        outfile << "rd output " << seed << "\n";
        outfile.flush();
    }
    else {
        std::cerr << "couldn't write, outfile appears to be closed \n";
    }
    std::mt19937 rng(seed);
    
    if (2.0*r_cut > min_box) {
        std::cerr << "LJ cutoff is greater than half the smallest dimension of the simulation box! \n";
    }

    //different ways to initialise positions - comment out as needed
    //if you want to initialise on an ordinary proton disordered ice lattice or as a fcc lattice for water simulation
    double* position_array = nullptr;
    if (from_file == "no") {
        position_array = initialise_positions(cell_type, rng, outfile);
    }
    else if (from_file == "yes") {
        position_array = initialise_fromfile(filename);
    }
    
    int i, bead;
    
    std::ofstream initfile("initial_positions.txt");
    if (initfile.is_open()) {
        for (bead=0; bead < n_beads; ++bead) {
            for (i=0; i < n_atoms; ++i) {
                initfile << position_array[bindex(bead, 3*i)] << " " << position_array[bindex(bead, 3*i+1)] << " " << position_array[bindex(bead, 3*i+2)] << "\n";
                initfile.flush();
            }
        }
    }
    else {
        std::cerr << "couldn't write, initfile appears to be closed \n";
    }
    initfile.close();
    
    if (n_beads == 1) {
        if (is_spline) {
            run_classical_wspline(position_array, outfile, rng);
        }
        else {
            run_classical_nospline(position_array, outfile, rng);
        }
    }
    else {
        if (is_spline) {
            run_sc_wspline(position_array, outfile, rng);
        }
        else {
            run_sc_nospline(position_array, outfile, rng);
        }
    }
    
    outfile.close();
    
    //output final positions to file
    std::ofstream finfile("final_positions.txt");
    if (finfile.is_open()) {
        for (bead=0; bead < n_beads; ++bead) {
            for (i=0; i < n_atoms; ++i) {
                finfile << position_array[bindex(bead, 3*i)] << " " << position_array[bindex(bead, 3*i+1)] << " " << position_array[bindex(bead, 3*i+2)] << "\n";
                finfile.flush();
            }
        }
    }
    else {
        std::cerr << "couldn't write, finfile appears to be closed \n";
    }
    finfile.close();
    
    delete[] position_array;
    return 0;
}
