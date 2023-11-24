//definition of thermostats used to run classical, CMD and PIMD simulations in
//the NVT ensemble

#include <vector>
#include <random>
#include "thermostat.h"
#include "qtip4pfparams.h"
#include "indexing.h"

//thermostat system to desired temperature for running in the NVT ensemble
//using a bussi thermostat
double bussi(double momentum_array[], const double mass_array[], std::mt19937& gen) {
    
    //keep track of the heat absorbed by the thermostat, initialised as zero
    int i, dof, dim;
    double R, alpha_sq, g2, R_add, alpha;
    double p2 = 0.0;
    std::normal_distribution<double> gauss_dev(0.0, 1.0);
    double heat = 0.0;
    
    //add initial kinetic energy to heat
    for (i=0; i < n_atoms; ++i) {
        for (dim=0; dim < 3; ++dim) {
            heat += 0.5 * momentum_array[3*i+dim] * momentum_array[3*i+dim] / mass_array[3*i];
        }
    }
    
    for (i=0; i < n_atoms; ++i) {
        for (dim=0; dim < 3; ++dim) {
            p2 += momentum_array[3*i+dim] * momentum_array[3*i+dim] / mass_array[3*i];
        }
    }
    
    R = gauss_dev(gen);
    dof = 3*n_atoms;
    if (dof%2 == 0) {
        std::gamma_distribution<double> gamma_dev((dof-2)/2, 1.0);
        R_add = gauss_dev(gen);
        g2 = 2.0 * gamma_dev(gen) + R_add * R_add;
    }
    else {
        std::gamma_distribution<double> gamma_dev((dof-1)/2, 1.0);
        g2 = 2.0 * gamma_dev(gen);
    }
    alpha_sq = (R*R + g2) / (beta * p2);
    alpha = sqrt(alpha_sq);
    
    for (i=0; i < n_dof; ++i) {
        momentum_array[i] *= alpha;
    }
    
    //subtract new kinetic energy from heat
    for (i=0; i < n_atoms; ++i) {
        for (dim=0; dim < 3; ++dim) {
            heat -= 0.5 * momentum_array[3*i+dim] * momentum_array[3*i+dim] / mass_array[3*i];
        }
    }
    
    return heat;
}

//thermostat system to desired temperature for running in the NVT ensemble
//using a PILE-G thermostat
double pile_g(double momentum_array[], const double mass_array[], std::mt19937& gen) {
    
    //keep track of the heat absorbed by the thermostat, initialised as zero
    int i, dof, dim, bead;
    double ke_centroid, R, alpha_sq, g2, R_add, alpha, sqrt_m_betapm1;
    double c1, c2, gamma_k;
    double prefactor = 4.0 / (beta_n * hbar);
    std::normal_distribution<double> gauss_dev(0.0, 1.0);
    double heat = 0.0;
    
    //we should already be in the normal mode representation
    //attach a global rescaling thermostat to the centroid
    //find the kinetic energy of the centroid mode
    
    ke_centroid = 0.0;
    for (i=0; i < n_atoms; ++i) {
        for (dim=0; dim < 3; ++dim) {
            ke_centroid += 0.5 * momentum_array[3*i+dim]*momentum_array[3*i+dim]/ mass_array[3*i+dim];
        }
    }
    
    //add initial kinetic energy of centroid to heat
    heat += ke_centroid;
    
    R = gauss_dev(gen);
    dof = 3*n_atoms;
    if (dof%2 == 0) {
        std::gamma_distribution<double> gamma_dev((dof-2)/2, 1.0);
        R_add = gauss_dev(gen);
        g2 = 2.0 * gamma_dev(gen) + R_add * R_add;
    }
    else {
        std::gamma_distribution<double> gamma_dev((dof-1)/2, 1.0);
        g2 = 2.0 * gamma_dev(gen);
    }
    alpha_sq = (R*R + g2) / (2.0 * beta_n * ke_centroid);
    alpha = sqrt(alpha_sq);
    
    for (i=0; i < n_dof; ++i) {
        momentum_array[i] *= alpha;
    }
    
    //subtract new kinetic energy from heat
    for (i=0; i < n_atoms; ++i) {
        for (dim=0; dim < 3; ++dim) {
            heat -= 0.5 * momentum_array[3*i+dim] * momentum_array[3*i+dim] / mass_array[3*i];
        }
    }
    
    //now apply a langevin thermostat to the other normal modes of the ring polymer
    for (bead=1; bead < n_beads; ++bead) {
        //add on original kinetic energy of this normal mode to the heat
        for (i=0; i < n_dof; ++i) {
            heat += 0.5 * momentum_array[bindex(bead, i)] * momentum_array[bindex(bead, i)] / mass_array[i];
        }
        gamma_k = prefactor * sin(bead * pi / n_beads);
        c1 = exp(-dt * gamma_k);
        c2 = sqrt(1.0 - c1 * c1);
        for (i=0; i < n_dof; ++i) {
            sqrt_m_betapm1 = sqrt(mass_array[i] / beta_n);
            momentum_array[bindex(bead, i)] = c1 * momentum_array[bindex(bead, i)] + sqrt_m_betapm1 * c2 * gauss_dev(gen);
        }
        //and subtract new kinetic energy from the heat
        for (i=0; i < n_dof; ++i) {
            heat -= 0.5 * momentum_array[bindex(bead, i)] * momentum_array[bindex(bead, i)] / mass_array[i];
        }
    }
    
    return heat;
}
