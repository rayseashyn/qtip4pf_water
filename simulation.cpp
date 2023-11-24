//functions to run classical, CMD or PIMD simulations
//these just calculate the relevant average kinetic and potential energies and
//output to a file, as this is just a basic simulation setup
//calculation of more complex quantities like the thermal conductivity or the 
//radial distribution functions are easily included by including their calculation 
//as an additional file and by additionally calculating the stress tensor in the force calculation

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <random>
#include <boost/math/interpolators/cubic_b_spline.hpp>
#include "update.h"
#include "thermostat.h"
#include "qtip4pfparams.h"
#include "neighbour.h"
#include "forces.h"
#include "path_integrals.h"
#include "indexing.h"

void run_classical_wspline(double position_array[], std::ofstream& outfile, std::mt19937& rng) {
    //create arrays to store necessary data
    double* momentum_array = new double[n_dof];
    double* force_array = new double[n_dof];
    double* disp_array = new double[3*n_mols];
    double* mass_array = new double[n_dof];
    std::array<double, 3> p_cofm;
    int* list = new int[list_size];
    int* point = new int[n_mols];
    double* charge_array = new double[n_atoms];
    int i, j, dim;
    double energy;
    double ke;
    
    //generate cubic spline for O-O
    //first data point to read in for O-O spline
    double oo_start = 0.0;
    double oo_stepsize = 0.01;
    //initialise empty vector to fill in with potential values from file
    std::vector<double> oo_pot;
    //filename to read in from - should be in the directory that the executable is run from
    //change this to whatever naming you'd like to use
    std::ifstream oo_infile("oo_spline.txt");
    std::string line;
    if (oo_infile.is_open()) {
        while (!oo_infile.eof()) {
            std::getline(oo_infile, line);
            if (line != "") {
                double entry = std::stod(line);
                oo_pot.push_back(entry);
            }
        }
    }
    else {
        std::cerr << "couldn't read, oo_infile appears to be closed \n";
    }
    oo_infile.close();
    double oo_spline_cutsq = oo_start + oo_stepsize * (oo_pot.size() - 1);
    oo_spline_cutsq = oo_spline_cutsq * oo_spline_cutsq;
    
    boost::math::cubic_b_spline<double> oo_spline(oo_pot.begin(), oo_pot.end(), oo_start, oo_stepsize);
    
    //generate cubic spline for O-H
    //first data point to read in for O-H spline
    double oh_start = 0.0;
    double oh_stepsize = 0.01;
    //initialise empty vector to fill in with potential values from file
    std::vector<double> oh_pot;
    //filename to read in from - should be in the directory that the executable is run from
    //change this to whatever naming you'd like to use
    std::ifstream oh_infile("oh_spline.txt");
    if (oh_infile.is_open()) {
        while (!oh_infile.eof()) {
            std::getline(oh_infile, line);
            if (line != "") {
                double entry = std::stod(line);
                oh_pot.push_back(entry);
            }
        }
    }
    else {
        std::cerr << "couldn't read, oo_infile appears to be closed \n";
    }
    oh_infile.close();
    
    boost::math::cubic_b_spline<double> oh_spline(oh_pot.begin(), oh_pot.end(), oh_start, oh_stepsize);
    double oh_spline_cutsq = oh_start + oh_stepsize * (oh_pot.size() - 1);
    oh_spline_cutsq = oh_spline_cutsq * oh_spline_cutsq;
    
    //generate cubic spline for H-H
    //first data point to read in for H-H spline
    double hh_start = 0.0;
    double hh_stepsize = 0.01;
    //initialise empty vector to fill in with potential values from file
    std::vector<double> hh_pot;
    //filename to read in from - should be in the directory that the executable is run from
    //change this to whatever naming you'd like to use
    std::ifstream hh_infile("hh_spline.txt");
    if (hh_infile.is_open()) {
        while (!hh_infile.eof()) {
            std::getline(hh_infile, line);
            if (line != "") {
                double entry = std::stod(line);
                hh_pot.push_back(entry);
            }
        }
    }
    else {
        std::cerr << "couldn't read, oo_infile appears to be closed \n";
    }
    hh_infile.close();
    
    boost::math::cubic_b_spline<double> hh_spline(hh_pot.begin(), hh_pot.end(), hh_start, hh_stepsize);
    double hh_spline_cutsq = hh_start + hh_stepsize * (hh_pot.size() - 1);
    hh_spline_cutsq = hh_spline_cutsq * hh_spline_cutsq;
    
    for (i=0; i < n_mols; ++i) {
        //store the charges on each atom
        charge_array[3*i] = qo;
        charge_array[3*i+1] = qh;
        charge_array[3*i+2] = qh;
        //store the masses of each atom (differently sized array)
        for (dim=0; dim < 3; ++dim) {
            mass_array[9*i+dim] = mo;
            mass_array[9*i+dim+3] = mh;
            mass_array[9*i+dim+6] = mh;
        }
    }

    //set up gaussian distribution
    //needed for generating momenta from boltzmann distribution
    std::normal_distribution<double> boltz_dev_ox(0.0, sd_ox);
    std::normal_distribution<double> boltz_dev_h(0.0, sd_h);

    //initialise displacements of oxygen atoms as zero
    for (i=0; i < n_mols; ++i) {
        for (dim=0; dim < 3; ++dim) {
            disp_array[3*i+dim] = 0.0;
        }
    }
    
    //set up the neighbour list
    neighbour_list(position_array, point, list, disp_array);
    
    for (dim=0; dim < 3; ++dim) {
        p_cofm[dim] = 0.0;
    }
    
    //initialise momenta
    for (i=0; i < n_mols; ++i) {
        for (dim=0; dim < 3; ++dim) {
            momentum_array[9*i+dim] = boltz_dev_ox(rng);
            momentum_array[9*i+dim+3] = boltz_dev_h(rng);
            momentum_array[9*i+dim+6] = boltz_dev_h(rng);
            p_cofm[dim] += momentum_array[9*i+dim] + momentum_array[9*i+dim+3] + momentum_array[9*i+dim+6];
        }
    }
    
    //calculate centre of mass momentum of the entire extended system
    for (dim=0; dim < 3; ++dim) {
        p_cofm[dim] /= n_atoms;
    }
    
    //subtract any centre of mass momentum from the entire extended system
    for (i=0; i < n_atoms; ++i) {
        for (dim=0; dim < 3; ++dim) {
            momentum_array[3*i+dim] -= p_cofm[dim];
        }
    }
    
    //fill in force array to begin verlet integration
    energy = get_forces_wspline(position_array, force_array, charge_array, point, list, oo_spline, oh_spline, hh_spline, oo_spline_cutsq, oh_spline_cutsq, hh_spline_cutsq);
    
    std::ofstream consfile("cons.txt");
    
    int n_out = n_equil / 10;
    
    for (i=0; i < n_equil; ++i) {
        energy = update_trajectory_class_wspline(position_array, momentum_array, mass_array, force_array, charge_array, point, list, disp_array, oo_spline, oh_spline, hh_spline, oo_spline_cutsq, oh_spline_cutsq, hh_spline_cutsq);
        //energy -= bussi(momentum_array, mass_array, rng);
        
        if (i%n_out==0) {
            if (outfile.is_open()) {
                outfile << "Equilibration is " << 10*i / n_out << "% done \n";
                outfile.flush();
            }
            else {
                std::cerr << "couldn't write, outfile appears to be closed \n";
            }
        }
        
        if (i%n_step==0) {
            ke = 0.0;
            for (j=0; j < n_atoms; ++j) {
                ke += 0.5 * (momentum_array[3*j]*momentum_array[3*j] + momentum_array[3*j+1]*momentum_array[3*j+1] + momentum_array[3*j+2]*momentum_array[3*j+2]) / mass_array[3*j];
            }
            if (consfile.is_open()) {
                consfile << energy + ke << "\n";
                consfile.flush();
            }
            else {
                std::cerr << "couldn't write, consfile appears to be closed \n";
            }
        }
    }
    
    //resample momenta
    for (dim=0; dim < 3; ++dim) {
        p_cofm[dim] = 0.0;
    }
    
    for (i=0; i < n_mols; ++i) {
        for (dim=0; dim < 3; ++dim) {
            momentum_array[9*i+dim] = boltz_dev_ox(rng);
            momentum_array[9*i+dim+3] = boltz_dev_h(rng);
            momentum_array[9*i+dim+6] = boltz_dev_h(rng);
            p_cofm[dim] += momentum_array[9*i+dim] + momentum_array[9*i+dim+3] + momentum_array[9*i+dim+6];
        }
    }
    
    for (dim=0; dim < 3; ++dim) {
        p_cofm[dim] /= n_atoms;
    }
    
    //subtract any centre of mass momentum from the system
    for (i=0; i < n_atoms; ++i) {
        for (dim=0; dim < 3; ++dim) {
            momentum_array[3*i+dim] -= p_cofm[dim];
        }
    }
    
    double avpot = 0.0;
    double avke = 0.0;
    
    n_out = n_samp / 10;
    
    for (i=0; i < n_samp; ++i) {
        energy = update_trajectory_class_wspline(position_array, momentum_array, mass_array, force_array, charge_array, point, list, disp_array, oo_spline, oh_spline, hh_spline, oo_spline_cutsq, oh_spline_cutsq, hh_spline_cutsq);
        //energy -= bussi(momentum_array, mass_array, rng);
        //calculate potential and centroid virial kinetic energy
        avpot += energy;
        ke = 0.0;
        for (j=0; j < n_atoms; ++j) {
            ke += 0.5 * (momentum_array[3*j]*momentum_array[3*j] + momentum_array[3*j+1]*momentum_array[3*j+1] + momentum_array[3*j+2]*momentum_array[3*j+2]) / mass_array[3*j];
        }
        avke += ke;
        
        if (i%n_out==0) {
            if (outfile.is_open()) {
                outfile << "Sampling is " << 10*i / n_out << "% done \n";
                outfile.flush();
            }
            else {
                std::cerr << "couldn't write, outfile appears to be closed \n";
            }
        }
        
        if (i%n_step==0) {
            if (consfile.is_open()) {
                consfile << energy + ke << "\n";
                consfile.flush();
            }
            else {
                std::cerr << "couldn't write, consfile appears to be closed \n";
            }
        }
    }
    
    //calculate average energy and centroid virial kinetic energy
    avpot /= n_samp;
    avke /= n_samp;
    
    //output average potential energy to file
    std::ofstream potfile("pot.txt");
    if (potfile.is_open()) {
        potfile << avpot << "\n";
        potfile.flush();
    }
    else {
        std::cerr << "couldn't write, potfile appears to be closed \n";
    }
    potfile.close();
    
    //output kinetic energy to file
    std::ofstream kinfile("kin.txt");
    if (kinfile.is_open()) {
        kinfile << avke << "\n";
        kinfile.flush();
    }
    else {
        std::cerr << "couldn't write, kinfile appears to be closed \n";
    }
    kinfile.close();
    //close any open files and finish up
    consfile.close();
    delete[] momentum_array;
    delete[] force_array;
    delete[] disp_array;
    delete[] mass_array;
    delete[] list;
    delete[] point;
    delete[] charge_array;
}

void run_classical_nospline(double position_array[], std::ofstream& outfile, std::mt19937& rng) {
    //create arrays to store necessary data
    double* momentum_array = new double[n_dof];
    double* force_array = new double[n_dof];
    double* disp_array = new double[3*n_mols];
    double* mass_array = new double[n_dof];
    std::array<double, 3> p_cofm;
    int* list = new int[list_size];
    int* point = new int[n_mols];
    double* charge_array = new double[n_atoms];
    int i, j, dim;
    double energy;
    double ke;
    
    for (i=0; i < n_mols; ++i) {
        //store the charges on each atom
        charge_array[3*i] = qo;
        charge_array[3*i+1] = qh;
        charge_array[3*i+2] = qh;
        //store the masses of each atom (differently sized array)
        for (dim=0; dim < 3; ++dim) {
            mass_array[9*i+dim] = mo;
            mass_array[9*i+dim+3] = mh;
            mass_array[9*i+dim+6] = mh;
        }
    }

    //set up gaussian distribution
    //needed for generating momenta from boltzmann distribution
    std::normal_distribution<double> boltz_dev_ox(0.0, sd_ox);
    std::normal_distribution<double> boltz_dev_h(0.0, sd_h);

    //initialise displacements of oxygen atoms as zero
    for (i=0; i < n_mols; ++i) {
        for (dim=0; dim < 3; ++dim) {
            disp_array[3*i+dim] = 0.0;
        }
    }
    
    //set up the neighbour list
    neighbour_list(position_array, point, list, disp_array);
    
    for (dim=0; dim < 3; ++dim) {
        p_cofm[dim] = 0.0;
    }
    
    //initialise momenta
    for (i=0; i < n_mols; ++i) {
        for (dim=0; dim < 3; ++dim) {
            momentum_array[9*i+dim] = boltz_dev_ox(rng);
            momentum_array[9*i+dim+3] = boltz_dev_h(rng);
            momentum_array[9*i+dim+6] = boltz_dev_h(rng);
            p_cofm[dim] += momentum_array[9*i+dim] + momentum_array[9*i+dim+3] + momentum_array[9*i+dim+6];
        }
    }
    
    //calculate centre of mass momentum of the entire extended system
    for (dim=0; dim < 3; ++dim) {
        p_cofm[dim] /= n_atoms;
    }
    
    //subtract any centre of mass momentum from the entire extended system
    for (i=0; i < n_atoms; ++i) {
        for (dim=0; dim < 3; ++dim) {
            momentum_array[3*i+dim] -= p_cofm[dim];
        }
    }
    
    //fill in force array to begin verlet integration
    energy = get_forces_nospline(position_array, force_array, charge_array, point, list);
    
    std::ofstream consfile("cons.txt");
    
    int n_out = n_equil / 10;
    
    for (i=0; i < n_equil; ++i) {
        energy = update_trajectory_class_nospline(position_array, momentum_array, mass_array, force_array, charge_array, point, list, disp_array);
        //energy -= bussi(momentum_array, mass_array, rng);
        
        if (i%n_out==0) {
            if (outfile.is_open()) {
                outfile << "Equilibration is " << 10*i / n_out << "% done \n";
                outfile.flush();
            }
            else {
                std::cerr << "couldn't write, outfile appears to be closed \n";
            }
        }
        
        if (i%n_step==0) {
            ke = 0.0;
            for (j=0; j < n_atoms; ++j) {
                ke += 0.5 * (momentum_array[3*j]*momentum_array[3*j] + momentum_array[3*j+1]*momentum_array[3*j+1] + momentum_array[3*j+2]*momentum_array[3*j+2]) / mass_array[3*j];
            }
            if (consfile.is_open()) {
                consfile << energy + ke << "\n";
                consfile.flush();
            }
            else {
                std::cerr << "couldn't write, consfile appears to be closed \n";
            }
        }
    }
    
    //resample momenta
    for (dim=0; dim < 3; ++dim) {
        p_cofm[dim] = 0.0;
    }
    
    for (i=0; i < n_mols; ++i) {
        for (dim=0; dim < 3; ++dim) {
            momentum_array[9*i+dim] = boltz_dev_ox(rng);
            momentum_array[9*i+dim+3] = boltz_dev_h(rng);
            momentum_array[9*i+dim+6] = boltz_dev_h(rng);
            p_cofm[dim] += momentum_array[9*i+dim] + momentum_array[9*i+dim+3] + momentum_array[9*i+dim+6];
        }
    }
    
    for (dim=0; dim < 3; ++dim) {
        p_cofm[dim] /= n_atoms;
    }
    
    //subtract any centre of mass momentum from the system
    for (i=0; i < n_atoms; ++i) {
        for (dim=0; dim < 3; ++dim) {
            momentum_array[3*i+dim] -= p_cofm[dim];
        }
    }
    
    double avpot = 0.0;
    double avke = 0.0;
    
    n_out = n_samp / 10;
    
    for (i=0; i < n_samp; ++i) {
        energy = update_trajectory_class_nospline(position_array, momentum_array, mass_array, force_array, charge_array, point, list, disp_array);
        //energy -= bussi(momentum_array, mass_array, rng);
        //calculate potential and centroid virial kinetic energy
        avpot += energy;
        ke = 0.0;
        for (j=0; j < n_atoms; ++j) {
            ke += 0.5 * (momentum_array[3*j]*momentum_array[3*j] + momentum_array[3*j+1]*momentum_array[3*j+1] + momentum_array[3*j+2]*momentum_array[3*j+2]) / mass_array[3*j];
        }
        avke += ke;
        
        if (i%n_out==0) {
            if (outfile.is_open()) {
                outfile << "Sampling is " << 10*i / n_out << "% done \n";
                outfile.flush();
            }
            else {
                std::cerr << "couldn't write, outfile appears to be closed \n";
            }
        }
        
        if (i%n_step==0) {
            if (consfile.is_open()) {
                consfile << energy + ke << "\n";
                consfile.flush();
            }
            else {
                std::cerr << "couldn't write, consfile appears to be closed \n";
            }
        }
    }
    
    //calculate average energy and centroid virial kinetic energy
    avpot /= n_samp;
    avke /= n_samp;
    
    //output average potential energy to file
    std::ofstream potfile("pot.txt");
    if (potfile.is_open()) {
        potfile << avpot << "\n";
        potfile.flush();
    }
    else {
        std::cerr << "couldn't write, potfile appears to be closed \n";
    }
    potfile.close();
    
    //output kinetic energy to file
    std::ofstream kinfile("kin.txt");
    if (kinfile.is_open()) {
        kinfile << avke << "\n";
        kinfile.flush();
    }
    else {
        std::cerr << "couldn't write, kinfile appears to be closed \n";
    }
    kinfile.close();
    //close any open files and finish up
    consfile.close();
    delete[] momentum_array;
    delete[] force_array;
    delete[] disp_array;
    delete[] mass_array;
    delete[] list;
    delete[] point;
    delete[] charge_array;
}

void run_sc_wspline(double position_array[], std::ofstream& outfile, std::mt19937& rng) {
    //create arrays to store necessary data
    double* momentum_array = new double[n_beads*n_dof];
    double* force_array = new double[n_beads*n_dof];
    double* fphys = new double[n_beads*n_dof];
    double* disp_array = new double[n_beads*3*n_mols];
    double* mass_array = new double[n_dof];
    std::array<double, 3> p_cofm;
    int* list = new int[n_beads*list_size];
    int* point = new int[n_beads*n_mols];
    double* charge_array = new double[n_atoms];
    int i, dim, bead;
    double energy;
    double heat = 0.0;
    double sc_pot = 0.0;
    
    //generate cubic spline for O-O
    //first data point to read in for O-O spline
    double oo_start = 0.0;
    double oo_stepsize = 0.01;
    //initialise empty vector to fill in with potential values from file
    std::vector<double> oo_pot;
    //filename to read in from - should be in the directory that the executable is run from
    //change this to whatever naming you'd like to use
    std::ifstream oo_infile("oo_spline.txt");
    std::string line;
    if (oo_infile.is_open()) {
        while (!oo_infile.eof()) {
            std::getline(oo_infile, line);
            if (line != "") {
                double entry = std::stod(line);
                oo_pot.push_back(entry);
            }
        }
    }
    else {
        std::cerr << "couldn't read, oo_infile appears to be closed \n";
    }
    oo_infile.close();
    double oo_spline_cutsq = oo_start + oo_stepsize * (oo_pot.size() - 1);
    oo_spline_cutsq = oo_spline_cutsq * oo_spline_cutsq;
    
    boost::math::cubic_b_spline<double> oo_spline(oo_pot.begin(), oo_pot.end(), oo_start, oo_stepsize);
    
    //generate cubic spline for O-H
    //first data point to read in for O-H spline
    double oh_start = 0.0;
    double oh_stepsize = 0.01;
    //initialise empty vector to fill in with potential values from file
    std::vector<double> oh_pot;
    //filename to read in from - should be in the directory that the executable is run from
    //change this to whatever naming you'd like to use
    std::ifstream oh_infile("oh_spline.txt");
    if (oh_infile.is_open()) {
        while (!oh_infile.eof()) {
            std::getline(oh_infile, line);
            if (line != "") {
                double entry = std::stod(line);
                oh_pot.push_back(entry);
            }
        }
    }
    else {
        std::cerr << "couldn't read, oo_infile appears to be closed \n";
    }
    oh_infile.close();
    
    boost::math::cubic_b_spline<double> oh_spline(oh_pot.begin(), oh_pot.end(), oh_start, oh_stepsize);
    double oh_spline_cutsq = oh_start + oh_stepsize * (oh_pot.size() - 1);
    oh_spline_cutsq = oh_spline_cutsq * oh_spline_cutsq;
    
    //generate cubic spline for H-H
    //first data point to read in for H-H spline
    double hh_start = 0.0;
    double hh_stepsize = 0.01;
    //initialise empty vector to fill in with potential values from file
    std::vector<double> hh_pot;
    //filename to read in from - should be in the directory that the executable is run from
    //change this to whatever naming you'd like to use
    std::ifstream hh_infile("hh_spline.txt");
    if (hh_infile.is_open()) {
        while (!hh_infile.eof()) {
            std::getline(hh_infile, line);
            if (line != "") {
                double entry = std::stod(line);
                hh_pot.push_back(entry);
            }
        }
    }
    else {
        std::cerr << "couldn't read, oo_infile appears to be closed \n";
    }
    hh_infile.close();
    
    boost::math::cubic_b_spline<double> hh_spline(hh_pot.begin(), hh_pot.end(), hh_start, hh_stepsize);
    double hh_spline_cutsq = hh_start + hh_stepsize * (hh_pot.size() - 1);
    hh_spline_cutsq = hh_spline_cutsq * hh_spline_cutsq;

    for (i=0; i < n_mols; ++i) {
        //store the charges on each atom
        charge_array[3*i] = qo;
        charge_array[3*i+1] = qh;
        charge_array[3*i+2] = qh;
        //store the masses of each atom (differently sized array)
        for (dim=0; dim < 3; ++dim) {
            mass_array[9*i+dim] = mo;
            mass_array[9*i+dim+3] = mh;
            mass_array[9*i+dim+6] = mh;
        }
    }

    //set up gaussian distribution
    //needed for generating momenta from boltzmann distribution
    std::normal_distribution<double> boltz_dev_ox(0.0, sd_ox);
    std::normal_distribution<double> boltz_dev_h(0.0, sd_h);

    //initialise displacements of oxygen atoms in each replica as zero
    for (bead=0; bead < n_beads; ++bead) {
        for (i=0; i < n_mols; ++i) {
            for (dim=0; dim < 3; ++dim) {
                disp_array[dindex(bead, 3*i+dim)] = 0.0;
            }
        }
    }
    
    //set up the neighbour list for each replica
    for (bead=0; bead < n_beads; ++bead) {
        neighbour_list(position_array, point, list, disp_array, bead);
    }
    
    for (dim=0; dim < 3; ++dim) {
        p_cofm[dim] = 0.0;
    }
    
    //initialise momenta
    for (bead=0; bead < n_beads; ++bead) {
        for (i=0; i < n_mols; ++i) {
            for (dim=0; dim < 3; ++dim) {
                momentum_array[bindex(bead, 9*i+dim)] = boltz_dev_ox(rng);
                momentum_array[bindex(bead, 9*i+dim+3)] = boltz_dev_h(rng);
                momentum_array[bindex(bead, 9*i+dim+6)] = boltz_dev_h(rng);
                p_cofm[dim] += momentum_array[bindex(bead, 9*i+dim)] + momentum_array[bindex(bead, 9*i+dim+3)] + momentum_array[bindex(bead, 9*i+dim+6)];
            }
        }
    }
    
    //calculate centre of mass momentum of the entire extended system
    for (dim=0; dim < 3; ++dim) {
        p_cofm[dim] /= n_atoms * n_beads;
    }
    
    //subtract any centre of mass momentum from the entire extended system
    for (bead=0; bead < n_beads; ++bead) {
        for (i=0; i < n_atoms; ++i) {
            for (dim=0; dim < 3; ++dim) {
                momentum_array[bindex(bead, 3*i+dim)] -= p_cofm[dim];
            }
        }
    }
    
    //fill in force array to begin verlet integration
    energy = force_eval_sc_wspline(position_array, force_array, charge_array, mass_array, point, list, oo_spline, oh_spline, hh_spline, oo_spline_cutsq, oh_spline_cutsq, hh_spline_cutsq, sc_pot, fphys);
    
    std::ofstream consfile("cons.txt");
    
    int n_out = n_equil / 10;
    //equilibrate system
    for (i=0; i < n_equil; ++i) {
        energy = update_sc_wspline(position_array, momentum_array, mass_array, force_array, fphys, charge_array, point, list, disp_array, rng, 0, heat, sc_pot, oo_spline, oh_spline, hh_spline, oo_spline_cutsq, oh_spline_cutsq, hh_spline_cutsq);
        
        if (i%n_out==0) {
            if (outfile.is_open()) {
                outfile << "Equilibration is " << 10*i / n_out << "% done \n";
                outfile.flush();
            }
            else {
                std::cerr << "couldn't write, outfile appears to be closed \n";
            }
        }
        
        if (i%n_step == 0) {
            if (consfile.is_open()) {
                consfile << get_conserved(position_array, momentum_array, mass_array, energy) + heat << "\n";
                consfile.flush();
            }
            else {
                std::cerr << "couldn't write, consfile appears to be closed \n";
            }
        }
    }
    
    //resample momenta
    for (dim=0; dim < 3; ++dim) {
        p_cofm[dim] = 0.0;
    }
    
    for (bead=0; bead < n_beads; ++bead) {
        for (i=0; i < n_mols; ++i) {
            for (dim=0; dim < 3; ++dim) {
                momentum_array[bindex(bead, 9*i+dim)] = boltz_dev_ox(rng);
                momentum_array[bindex(bead, 9*i+dim+3)] = boltz_dev_h(rng);
                momentum_array[bindex(bead, 9*i+dim+6)] = boltz_dev_h(rng);
                p_cofm[dim] += momentum_array[bindex(bead, 9*i+dim)] + momentum_array[bindex(bead, 9*i+dim+3)] + momentum_array[bindex(bead, 9*i+dim+6)];
            }
        }
    }
    
    for (dim=0; dim < 3; ++dim) {
        p_cofm[dim] /= n_atoms * n_beads;
    }
    
    //subtract any centre of mass momentum from the entire extended system
    for (bead=0; bead < n_beads; ++bead) {
        for (i=0; i < n_atoms; ++i) {
            for (dim=0; dim < 3; ++dim) {
                momentum_array[bindex(bead, 3*i+dim)] -= p_cofm[dim];
            }
        }
    }
    
    double pot = 0.0;
    double kcv = 0.0;
    n_out = n_samp / 10;
    for (i=0; i < n_samp; ++i) {
        energy = update_sc_wspline(position_array, momentum_array, mass_array, force_array, fphys, charge_array, point, list, disp_array, rng, 1, heat, sc_pot, oo_spline, oh_spline, hh_spline, oo_spline_cutsq, oh_spline_cutsq, hh_spline_cutsq);
        
        if (i%n_out==0) {
            if (outfile.is_open()) {
                outfile << "Sampling is " << 10*i / n_out << "% done \n";
                outfile.flush();
            }
            else {
                std::cerr << "couldn't write, outfile appears to be closed \n";
            }
        }
        
        if (i%n_step == 0) {
            if (consfile.is_open()) {
                consfile << get_conserved(position_array, momentum_array, mass_array, energy) + heat << "\n";
                consfile.flush();
            }
            else {
                std::cerr << "couldn't write, consfile appears to be closed \n";
            }
        }
        //calculate potential and centroid virial kinetic energy
        pot += sc_pot;
        kcv += get_kcv(position_array, fphys);
    }
    
    //calculate average energy and centroid virial kinetic energy
    pot /= n_samp;
    kcv /= n_samp;
    
    //output average potential energy to file
    std::ofstream potfile("pot_sc.txt");
    if (potfile.is_open()) {
        potfile << pot << "\n";
        potfile.flush();
    }
    else {
        std::cerr << "couldn't write, potfile appears to be closed \n";
    }
    potfile.close();
    
    //output kinetic energy to file
    std::ofstream kcvfile("kcv.txt");
    if (kcvfile.is_open()) {
        kcvfile << kcv << "\n";
        kcvfile.flush();
    }
    else {
        std::cerr << "couldn't write, kcvfile appears to be closed \n";
    }
    kcvfile.close();
    //close any open files and finish up
    consfile.close();
    delete[] momentum_array;
    delete[] force_array;
    delete[] fphys;
    delete[] disp_array;
    delete[] mass_array;
    delete[] list;
    delete[] point;
    delete[] charge_array;
}

void run_sc_nospline(double position_array[], std::ofstream& outfile, std::mt19937& rng) {
    //create arrays to store necessary data
    double* momentum_array = new double[n_beads*n_dof];
    double* force_array = new double[n_beads*n_dof];
    double* fphys = new double[n_beads*n_dof];
    double* disp_array = new double[n_beads*3*n_mols];
    double* mass_array = new double[n_dof];
    std::array<double, 3> p_cofm;
    int* list = new int[n_beads*list_size];
    int* point = new int[n_beads*n_mols];
    double* charge_array = new double[n_atoms];
    int i, dim, bead;
    double energy;
    double heat = 0.0;
    double sc_pot = 0.0;

    for (i=0; i < n_mols; ++i) {
        //store the charges on each atom
        charge_array[3*i] = qo;
        charge_array[3*i+1] = qh;
        charge_array[3*i+2] = qh;
        //store the masses of each atom (differently sized array)
        for (dim=0; dim < 3; ++dim) {
            mass_array[9*i+dim] = mo;
            mass_array[9*i+dim+3] = mh;
            mass_array[9*i+dim+6] = mh;
        }
    }

    //set up gaussian distribution
    //needed for generating momenta from boltzmann distribution
    std::normal_distribution<double> boltz_dev_ox(0.0, sd_ox);
    std::normal_distribution<double> boltz_dev_h(0.0, sd_h);

    //initialise displacements of oxygen atoms in each replica as zero
    for (bead=0; bead < n_beads; ++bead) {
        for (i=0; i < n_mols; ++i) {
            for (dim=0; dim < 3; ++dim) {
                disp_array[dindex(bead, 3*i+dim)] = 0.0;
            }
        }
    }
    
    //set up the neighbour list for each replica
    for (bead=0; bead < n_beads; ++bead) {
        neighbour_list(position_array, point, list, disp_array, bead);
    }
    
    for (dim=0; dim < 3; ++dim) {
        p_cofm[dim] = 0.0;
    }
    
    //initialise momenta
    for (bead=0; bead < n_beads; ++bead) {
        for (i=0; i < n_mols; ++i) {
            for (dim=0; dim < 3; ++dim) {
                momentum_array[bindex(bead, 9*i+dim)] = boltz_dev_ox(rng);
                momentum_array[bindex(bead, 9*i+dim+3)] = boltz_dev_h(rng);
                momentum_array[bindex(bead, 9*i+dim+6)] = boltz_dev_h(rng);
                p_cofm[dim] += momentum_array[bindex(bead, 9*i+dim)] + momentum_array[bindex(bead, 9*i+dim+3)] + momentum_array[bindex(bead, 9*i+dim+6)];
            }
        }
    }
    
    //calculate centre of mass momentum of the entire extended system
    for (dim=0; dim < 3; ++dim) {
        p_cofm[dim] /= n_atoms * n_beads;
    }
    
    //subtract any centre of mass momentum from the entire extended system
    for (bead=0; bead < n_beads; ++bead) {
        for (i=0; i < n_atoms; ++i) {
            for (dim=0; dim < 3; ++dim) {
                momentum_array[bindex(bead, 3*i+dim)] -= p_cofm[dim];
            }
        }
    }
    
    //fill in force array to begin verlet integration
    energy = force_eval_sc_nospline(position_array, force_array, charge_array, mass_array, point, list,  sc_pot, fphys);
    
    std::ofstream consfile("cons.txt");
    
    int n_out = n_equil / 10;
    //equilibrate system
    for (i=0; i < n_equil; ++i) {
        energy = update_sc_nospline(position_array, momentum_array, mass_array, force_array, fphys, charge_array, point, list, disp_array, rng, 1, heat, sc_pot);
        
        if (i%n_out==0) {
            if (outfile.is_open()) {
                outfile << "Equilibration is " << 10*i / n_out << "% done \n";
                outfile.flush();
            }
            else {
                std::cerr << "couldn't write, outfile appears to be closed \n";
            }
        }
        
        if (i%n_step == 0) {
            if (consfile.is_open()) {
                consfile << get_conserved(position_array, momentum_array, mass_array, energy) + heat << "\n";
                consfile.flush();
            }
            else {
                std::cerr << "couldn't write, consfile appears to be closed \n";
            }
        }
    }
    
    //resample momenta
    for (dim=0; dim < 3; ++dim) {
        p_cofm[dim] = 0.0;
    }
    
    for (bead=0; bead < n_beads; ++bead) {
        for (i=0; i < n_mols; ++i) {
            for (dim=0; dim < 3; ++dim) {
                momentum_array[bindex(bead, 9*i+dim)] = boltz_dev_ox(rng);
                momentum_array[bindex(bead, 9*i+dim+3)] = boltz_dev_h(rng);
                momentum_array[bindex(bead, 9*i+dim+6)] = boltz_dev_h(rng);
                p_cofm[dim] += momentum_array[bindex(bead, 9*i+dim)] + momentum_array[bindex(bead, 9*i+dim+3)] + momentum_array[bindex(bead, 9*i+dim+6)];
            }
        }
    }
    
    for (dim=0; dim < 3; ++dim) {
        p_cofm[dim] /= n_atoms * n_beads;
    }
    
    //subtract any centre of mass momentum from the entire extended system
    for (bead=0; bead < n_beads; ++bead) {
        for (i=0; i < n_atoms; ++i) {
            for (dim=0; dim < 3; ++dim) {
                momentum_array[bindex(bead, 3*i+dim)] -= p_cofm[dim];
            }
        }
    }
    
    double pot = 0.0;
    double kcv = 0.0;
    n_out = n_samp / 10;
    for (i=0; i < n_samp; ++i) {
        energy = update_sc_nospline(position_array, momentum_array, mass_array, force_array, fphys, charge_array, point, list, disp_array, rng, 0, heat, sc_pot);
        
        if (i%n_out==0) {
            if (outfile.is_open()) {
                outfile << "Sampling is " << 10*i / n_out << "% done \n";
                outfile.flush();
            }
            else {
                std::cerr << "couldn't write, outfile appears to be closed \n";
            }
        }
        
        if (i%n_step == 0) {
            if (consfile.is_open()) {
                consfile << get_conserved(position_array, momentum_array, mass_array, energy) + heat << "\n";
                consfile.flush();
            }
            else {
                std::cerr << "couldn't write, consfile appears to be closed \n";
            }
        }
        //calculate potential and centroid virial kinetic energy
        pot += sc_pot;
        kcv += get_kcv(position_array, fphys);
    }
    
    //calculate average energy and centroid virial kinetic energy
    pot /= n_samp;
    kcv /= n_samp;
    
    //output average potential energy to file
    std::ofstream potfile("pot_sc.txt");
    if (potfile.is_open()) {
        potfile << pot << "\n";
        potfile.flush();
    }
    else {
        std::cerr << "couldn't write, potfile appears to be closed \n";
    }
    potfile.close();
    
    //output kinetic energy to file
    std::ofstream kcvfile("kcv.txt");
    if (kcvfile.is_open()) {
        kcvfile << kcv << "\n";
        kcvfile.flush();
    }
    else {
        std::cerr << "couldn't write, kcvfile appears to be closed \n";
    }
    kcvfile.close();
    //close any open files and finish up
    consfile.close();
    delete[] momentum_array;
    delete[] force_array;
    delete[] fphys;
    delete[] disp_array;
    delete[] mass_array;
    delete[] list;
    delete[] point;
    delete[] charge_array;
}

