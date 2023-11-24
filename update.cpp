//functions to update the configurations of classical, CMD and PIMD simulations
//for a timestep

#include <boost/math/interpolators/cubic_b_spline.hpp>
#include <random>
#include <iostream>
#include "update.h"
#include "neighbour.h"
#include "forces.h"
#include "qtip4pfparams.h"
#include "indexing.h"
#include "fourier.h"
#include "thermostat.h"

//update positions and velocities using a velocity verlet algorithm
double update_trajectory_class_wspline(double position_array[], double momentum_array[], const double mass_array[], double force_array[], double charge_array[], int point[], int list[], double disp_array[], const boost::math::cubic_b_spline<double>& oo_spline, const boost::math::cubic_b_spline<double>& oh_spline, const boost::math::cubic_b_spline<double>& hh_spline, const double oo_spline_cutsq, const double oh_spline_cutsq, const double hh_spline_cutsq) {
    
    //indexing parameters
    int i, dim;
    //potential energy
    double energy, disp_sq, dispmax, dispmax2;

    for (i=0; i < n_atoms; ++i) {
        for (dim=0; dim < 3; ++dim) {
            momentum_array[3*i+dim] += halfdt * force_array[3*i+dim];
            position_array[3*i+dim] += momentum_array[3*i+dim] * dt / mass_array[3*i+dim];
            if (i%3==0) {
                disp_array[i+dim] += momentum_array[3*i+dim] * dt / mass_array[3*i+dim];
            }
        }
    }
    
    //check if the neighbour list needs updating
    //initialise the two largest displacements in this replica as zero
    dispmax = 0.0;
    dispmax2 = 0.0;
    //find the squared distance travelled by each atom in this replica since the list was last updated
    for (i=0; i < n_mols; ++i) {
        disp_sq = disp_array[3*i]*disp_array[3*i] + disp_array[3*i+1]*disp_array[3*i+1] + disp_array[3*i+2]*disp_array[3*i+2];
        //check if this is larger than the currently stored largest squared displacement
        if (disp_sq > dispmax) {
            //if so, set the old largest as the new second largest
            dispmax2 = dispmax;
            //and store this as the new largest squared displacement
            dispmax = disp_sq;
        }
        //if not, check if it is larger than the second largest currently stored squared displacement
        else if (disp_sq > dispmax2) {
            //store this as the new second largest squared displacement
            dispmax2 = disp_sq;
        }
    }
    //compute sqrt of the squared displacements stored
    dispmax = sqrt(dispmax);
    dispmax2 = sqrt(dispmax2);
    
    //check if the two largest displacements in this replica is larger than the "skin" of the neighbour list
    if (dispmax + dispmax2 > r_skin - r_cut) {
        //if so, the neighbour list for this replica needs updating
        neighbour_list(position_array, point, list, disp_array);
    }
    
    //calculate forces and energy as well as stress tensor/site energies at the updated positions
    energy = get_forces_wspline(position_array, force_array, charge_array, point, list, oo_spline, oh_spline, hh_spline, oo_spline_cutsq, oh_spline_cutsq, hh_spline_cutsq);

    //final update of momenta via classical forces for half a timestep
    for (i=0; i < n_atoms; ++i) {
        for (dim=0; dim < 3; ++dim) {
            momentum_array[3*i+dim] += halfdt * force_array[3*i+dim];
        }
    }
    
    //return the potential energy of the system
    return energy;
}

//update positions and velocities using a velocity verlet algorithm
double update_trajectory_class_nospline(double position_array[], double momentum_array[], const double mass_array[], double force_array[], double charge_array[], int point[], int list[], double disp_array[]) {
    
    //indexing parameters
    int i, dim;
    //potential energy
    double energy, disp_sq, dispmax, dispmax2;

    for (i=0; i < n_atoms; ++i) {
        for (dim=0; dim < 3; ++dim) {
            momentum_array[3*i+dim] += halfdt * force_array[3*i+dim];
            position_array[3*i+dim] += momentum_array[3*i+dim] * dt / mass_array[3*i+dim];
            if (i%3==0) {
                disp_array[i+dim] += momentum_array[3*i+dim] * dt / mass_array[3*i+dim];
            }
        }
    }
    
    //check if the neighbour list needs updating
    //initialise the two largest displacements in this replica as zero
    dispmax = 0.0;
    dispmax2 = 0.0;
    //find the squared distance travelled by each atom in this replica since the list was last updated
    for (i=0; i < n_mols; ++i) {
        disp_sq = disp_array[3*i]*disp_array[3*i] + disp_array[3*i+1]*disp_array[3*i+1] + disp_array[3*i+2]*disp_array[3*i+2];
        //check if this is larger than the currently stored largest squared displacement
        if (disp_sq > dispmax) {
            //if so, set the old largest as the new second largest
            dispmax2 = dispmax;
            //and store this as the new largest squared displacement
            dispmax = disp_sq;
        }
        //if not, check if it is larger than the second largest currently stored squared displacement
        else if (disp_sq > dispmax2) {
            //store this as the new second largest squared displacement
            dispmax2 = disp_sq;
        }
    }
    //compute sqrt of the squared displacements stored
    dispmax = sqrt(dispmax);
    dispmax2 = sqrt(dispmax2);
    
    //check if the two largest displacements in this replica is larger than the "skin" of the neighbour list
    if (dispmax + dispmax2 > r_skin - r_cut) {
        //if so, the neighbour list for this replica needs updating
        neighbour_list(position_array, point, list, disp_array);
    }
    
    //calculate forces and energy as well as stress tensor/site energies at the updated positions
    energy = get_forces_nospline(position_array, force_array, charge_array, point, list);

    //final update of momenta via classical forces for half a timestep
    for (i=0; i < n_atoms; ++i) {
        for (dim=0; dim < 3; ++dim) {
            momentum_array[3*i+dim] += halfdt * force_array[3*i+dim];
        }
    }
    
    //return the potential energy of the system
    return energy;
}

double update_sc_wspline(double position_array[], double momentum_array[], const double mass_array[], double force_array[], double fphys[], double charge_array[], int point[], int list[], double disp_array[], std::mt19937& rng, bool is_nvt, double& heat, double& sc_pot, const boost::math::cubic_b_spline<double>& oo_spline, const boost::math::cubic_b_spline<double>& oh_spline, const boost::math::cubic_b_spline<double>& hh_spline, const double oo_spline_cutsq, const double oh_spline_cutsq, const double hh_spline_cutsq) {
    
    //indexing parameters
    int i, dim, bead;
    //potential energy
    double energy, disp_sq, omega_jsq, scale, s, temp, dispmax, dispmax2;
    double prefactor = 2.0 / (hbar * beta_n);
    double prefactor_sq = prefactor * prefactor;
    
    //propagation of momenta by classical forces for one half time step
    for (bead=0; bead < n_beads; ++bead) {
        for (i=0; i < n_atoms; ++i) {
            for (dim=0; dim < 3; ++dim) {
                momentum_array[bindex(bead, 3*i+dim)] += halfdt * force_array[bindex(bead, 3*i+dim)];
            }
        }
    }
    
    //convert arrays to the normal mode representation
    ftransform(position_array, n_dof, 1);
    ftransform(momentum_array, n_dof, 1);
    ftransform(disp_array, 3*n_mols, 1);
    
    //propagate positions and momenta using the cayley transformed propagator
    //this allows for better exploration of phase space and larger time steps to be used
    //https://aip.scitation.org/doi/pdf/10.1063/1.5120282
    //since we are using the BCOCB splitting of the integrator, we are only propagating for half a time step
    //we actually use the square root of this evolution, as discussed in
    //https://aip.scitation.org/doi/pdf/10.1063/1.5134810
    //in order to maintain "strong stability" - as sqrt(cay(A*dt)) != cay(A*dt/2), which is not an issue for the "normal" exponential propagators
    for (bead=0; bead < n_beads; ++bead) {
        s = sin(bead * piovern);
        omega_jsq = prefactor_sq * s * s;
        scale = 1.0 / sqrt(1.0 + 0.25 * omega_jsq * dt * dt);
        for (i=0; i < n_atoms; ++i) {
            for (dim=0; dim < 3; ++dim) {
                temp = halfdt * momentum_array[bindex(bead, 3*i+dim)] / mass_array[3*i] + position_array[bindex(bead, 3*i+dim)];
                temp *= scale;
                momentum_array[bindex(bead, 3*i+dim)] = momentum_array[bindex(bead, 3*i+dim)] - halfdt * omega_jsq * position_array[bindex(bead, 3*i+dim)] * mass_array[3*i];
                momentum_array[bindex(bead, 3*i+dim)] *= scale;
                //check if we are propagating an oxygen atom
                if (i%3==0) {
                    //track changes in positions
                    //this will be converted to the bead representation later on
                    disp_array[dindex(bead, i+dim)] += temp - position_array[bindex(bead, 3*i+dim)];
                }
                position_array[bindex(bead, 3*i+dim)] = temp;
            }
        }
    }
    
    //check if we are supposed to be in the NVT ensemble
    if (is_nvt) {
        //if so, apply thermostat in the middle
        heat -= pile_g(momentum_array, mass_array, rng);
    }
    
    //propagate positions and momenta using the cayley transformed propagator for half a time step
    for (bead=0; bead < n_beads; ++bead) {
        s = sin(bead * piovern);
        omega_jsq = prefactor_sq * s * s;
        scale = 1.0 / sqrt(1.0 + 0.25 * omega_jsq * dt * dt);
        for (i=0; i < n_atoms; ++i) {
            for (dim=0; dim < 3; ++dim) {
                temp = halfdt * momentum_array[bindex(bead, 3*i+dim)] / mass_array[3*i] + position_array[bindex(bead, 3*i+dim)];
                temp *= scale;
                momentum_array[bindex(bead, 3*i+dim)] = momentum_array[bindex(bead, 3*i+dim)] - halfdt * omega_jsq * position_array[bindex(bead, 3*i+dim)] * mass_array[3*i];
                momentum_array[bindex(bead, 3*i+dim)] *= scale;
                //check if we are propagating an oxygen atom
                if (i%3==0) {
                    //track changes in positions
                    //this will be converted to the bead representation later on
                    disp_array[dindex(bead, i+dim)] += temp - position_array[bindex(bead, 3*i+dim)];
                }
                position_array[bindex(bead, 3*i+dim)] = temp;
            }
        }
    }
    //convert back to the bead representation
    ftransform(position_array, n_dof, 0);
    ftransform(momentum_array, n_dof, 0);
    ftransform(disp_array, 3*n_mols, 0);
    
    //check if the neighbour list for any of the replicas needs updating
    //loop over all replicas in the extended system
    for (bead=0; bead < n_beads; ++bead) {
        //initialise the two largest displacements in this replica as zero
        dispmax = 0.0;
        dispmax2 = 0.0;
        //find the squared distance travelled by each atom in this replica since the list was last updated
        for (i=0; i < n_mols; ++i) {
            disp_sq = disp_array[dindex(bead, 3*i)]*disp_array[dindex(bead, 3*i)] + disp_array[dindex(bead, 3*i+1)]*disp_array[dindex(bead, 3*i+1)] + disp_array[dindex(bead, 3*i+2)]*disp_array[dindex(bead, 3*i+2)];
            //check if this is larger than the currently stored largest squared displacement
            if (disp_sq > dispmax) {
                //if so, set the old largest as the new second largest
                dispmax2 = dispmax;
                //and store this as the new largest squared displacement
                dispmax = disp_sq;
            }
            //if not, check if it is larger than the second largest currently stored squared displacement
            else if (disp_sq > dispmax2) {
                //store this as the new second largest squared displacement
                dispmax2 = disp_sq;
            }
        }
        //compute sqrt of the squared displacements stored
        dispmax = sqrt(dispmax);
        dispmax2 = sqrt(dispmax2);
        
        //check if the two largest displacements in this replica is larger than the "skin" of the neighbour list
        if (dispmax + dispmax2 > r_skin - r_cut) {
            //if so, the neighbour list for this replica needs updating
            neighbour_list(position_array, point, list, disp_array, bead);
        }
    }
    
    //calculate forces and energy as well as stress tensor/site energies at the updated positions
    energy = force_eval_sc_wspline(position_array, force_array, charge_array, mass_array, point, list, oo_spline, oh_spline, hh_spline, oo_spline_cutsq, oh_spline_cutsq, hh_spline_cutsq, sc_pot, fphys);

    //final update of momenta via classical forces for half a timestep
    for (bead=0; bead < n_beads; ++bead) {
        for (i=0; i < n_atoms; ++i) {
            for (dim=0; dim < 3; ++dim) {
                momentum_array[bindex(bead, 3*i+dim)] += halfdt * force_array[bindex(bead, 3*i+dim)];
            }
        }
    }
    
    //return the potential energy of the system
    return energy;
}

double update_sc_nospline(double position_array[], double momentum_array[], const double mass_array[], double force_array[], double fphys[], double charge_array[], int point[], int list[], double disp_array[], std::mt19937& rng, bool is_nvt, double& heat, double& sc_pot) {
    
    //indexing parameters
    int i, dim, bead;
    //potential energy
    double energy, disp_sq, omega_jsq, scale, s, temp, dispmax, dispmax2;
    double prefactor = 2.0 / (hbar * beta_n);
    double prefactor_sq = prefactor * prefactor;
    
    //propagation of momenta by classical forces for one half time step
    for (bead=0; bead < n_beads; ++bead) {
        for (i=0; i < n_atoms; ++i) {
            for (dim=0; dim < 3; ++dim) {
                momentum_array[bindex(bead, 3*i+dim)] += halfdt * force_array[bindex(bead, 3*i+dim)];
            }
        }
    }
    
    //convert arrays to the normal mode representation
    ftransform(position_array, n_dof, 1);
    ftransform(momentum_array, n_dof, 1);
    ftransform(disp_array, 3*n_mols, 1);
    
    //propagate positions and momenta using the cayley transformed propagator
    //this allows for better exploration of phase space and larger time steps to be used
    //https://aip.scitation.org/doi/pdf/10.1063/1.5120282
    //since we are using the BCOCB splitting of the integrator, we are only propagating for half a time step
    //we actually use the square root of this evolution, as discussed in
    //https://aip.scitation.org/doi/pdf/10.1063/1.5134810
    //in order to maintain "strong stability" - as sqrt(cay(A*dt)) != cay(A*dt/2), which is not an issue for the "normal" exponential propagators
    for (bead=0; bead < n_beads; ++bead) {
        s = sin(bead * piovern);
        omega_jsq = prefactor_sq * s * s;
        scale = 1.0 / sqrt(1.0 + 0.25 * omega_jsq * dt * dt);
        for (i=0; i < n_atoms; ++i) {
            for (dim=0; dim < 3; ++dim) {
                temp = halfdt * momentum_array[bindex(bead, 3*i+dim)] / mass_array[3*i] + position_array[bindex(bead, 3*i+dim)];
                temp *= scale;
                momentum_array[bindex(bead, 3*i+dim)] = momentum_array[bindex(bead, 3*i+dim)] - halfdt * omega_jsq * position_array[bindex(bead, 3*i+dim)] * mass_array[3*i];
                momentum_array[bindex(bead, 3*i+dim)] *= scale;
                //check if we are propagating an oxygen atom
                if (i%3==0) {
                    //track changes in positions
                    //this will be converted to the bead representation later on
                    disp_array[dindex(bead, i+dim)] += temp - position_array[bindex(bead, 3*i+dim)];
                }
                position_array[bindex(bead, 3*i+dim)] = temp;
            }
        }
    }
    
    //check if we are supposed to be in the NVT ensemble
    if (is_nvt) {
        //if so, apply thermostat in the middle
        heat -= pile_g(momentum_array, mass_array, rng);
    }
    
    //propagate positions and momenta using the cayley transformed propagator for half a time step
    for (bead=0; bead < n_beads; ++bead) {
        s = sin(bead * piovern);
        omega_jsq = prefactor_sq * s * s;
        scale = 1.0 / sqrt(1.0 + 0.25 * omega_jsq * dt * dt);
        for (i=0; i < n_atoms; ++i) {
            for (dim=0; dim < 3; ++dim) {
                temp = halfdt * momentum_array[bindex(bead, 3*i+dim)] / mass_array[3*i] + position_array[bindex(bead, 3*i+dim)];
                temp *= scale;
                momentum_array[bindex(bead, 3*i+dim)] = momentum_array[bindex(bead, 3*i+dim)] - halfdt * omega_jsq * position_array[bindex(bead, 3*i+dim)] * mass_array[3*i];
                momentum_array[bindex(bead, 3*i+dim)] *= scale;
                //check if we are propagating an oxygen atom
                if (i%3==0) {
                    //track changes in positions
                    //this will be converted to the bead representation later on
                    disp_array[dindex(bead, i+dim)] += temp - position_array[bindex(bead, 3*i+dim)];
                }
                position_array[bindex(bead, 3*i+dim)] = temp;
            }
        }
    }
    //convert back to the bead representation
    ftransform(position_array, n_dof, 0);
    ftransform(momentum_array, n_dof, 0);
    ftransform(disp_array, 3*n_mols, 0);
    
    //check if the neighbour list for any of the replicas needs updating
    //loop over all replicas in the extended system
    for (bead=0; bead < n_beads; ++bead) {
        //initialise the two largest displacements in this replica as zero
        dispmax = 0.0;
        dispmax2 = 0.0;
        //find the squared distance travelled by each atom in this replica since the list was last updated
        for (i=0; i < n_mols; ++i) {
            disp_sq = disp_array[dindex(bead, 3*i)]*disp_array[dindex(bead, 3*i)] + disp_array[dindex(bead, 3*i+1)]*disp_array[dindex(bead, 3*i+1)] + disp_array[dindex(bead, 3*i+2)]*disp_array[dindex(bead, 3*i+2)];
            //check if this is larger than the currently stored largest squared displacement
            if (disp_sq > dispmax) {
                //if so, set the old largest as the new second largest
                dispmax2 = dispmax;
                //and store this as the new largest squared displacement
                dispmax = disp_sq;
            }
            //if not, check if it is larger than the second largest currently stored squared displacement
            else if (disp_sq > dispmax2) {
                //store this as the new second largest squared displacement
                dispmax2 = disp_sq;
            }
        }
        //compute sqrt of the squared displacements stored
        dispmax = sqrt(dispmax);
        dispmax2 = sqrt(dispmax2);
        
        //check if the two largest displacements in this replica is larger than the "skin" of the neighbour list
        if (dispmax + dispmax2 > r_skin - r_cut) {
            //if so, the neighbour list for this replica needs updating
            neighbour_list(position_array, point, list, disp_array, bead);
        }
    }
    
    //calculate forces and energy as well as stress tensor/site energies at the updated positions
    energy = force_eval_sc_nospline(position_array, force_array, charge_array, mass_array, point, list, sc_pot, fphys);

    //final update of momenta via classical forces for half a timestep
    for (bead=0; bead < n_beads; ++bead) {
        for (i=0; i < n_atoms; ++i) {
            for (dim=0; dim < 3; ++dim) {
                momentum_array[bindex(bead, 3*i+dim)] += halfdt * force_array[bindex(bead, 3*i+dim)];
            }
        }
    }
    
    //return the potential energy of the system
    return energy;
}
