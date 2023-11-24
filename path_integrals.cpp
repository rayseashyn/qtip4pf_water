//functions to calculate useful quantities from path integral simulations
//such as the centroid position, quantum kinetic energy estimator and
//the conserved quantity from the ring polymer hamiltonian in a NVE path integral simulation

#include "path_integrals.h"
#include "qtip4pfparams.h"
#include "indexing.h"
#include "fourier.h"

double* get_q_centroid(const double input_array[], const int length) {
    double* q_centroid = new double[length];
    int bead, i;
    double qc_dof;
    
    for (i=0; i < length; ++i) {
        qc_dof = 0.0;
        for (bead=0; bead < n_beads; ++bead) {
            qc_dof += input_array[bead*length + i];
        }
        
        qc_dof /= static_cast <float> (n_beads);
        q_centroid[i] = qc_dof;
    }
    
    return q_centroid;
}

//returns the centroid virial estimator for the kinetic energy
double get_kcv(const double position_array[], const double force_array[]) {
    
    int bead, j;
    double tq, del_q, tv;
    double* q_centroid = get_q_centroid(position_array, n_dof);
    
    tq = 0.0;
    //increment loop by two as only even beads contribute for the OP estimator
    for (bead=0; bead < n_beads; bead+=2) {
        for (j=0; j < n_dof; ++j) {
            del_q = position_array[bindex(bead, j)] - q_centroid[j];
            tq -= del_q * force_array[bindex(bead, j)];
        }
    }
    tq /= n_beads;
    tv = 0.5 * n_dof / beta + tq;
    
    delete[] q_centroid;
    return tv;
}

//returns the ring polymer hamiltonian conserved energy quantity
//takes in the positions to compute the spring contribution in the normal mode basis
//and the velocities to determine the kinetic energy in the normal mode basis
//the energy input to this function should be the total classical potential energy from all
//of the replicas
double get_conserved(double position_array[], double momentum_array[], const double mass_array[], double energy) {
    
    int bead, j;
    double nm_kin, v_spring, q2, prefactor, s, omega_jsq;
    double piovern = pi / n_beads;
    
    ftransform(position_array, n_dof, 1);
    ftransform(momentum_array, n_dof, 1);
    nm_kin = 0.0;
    v_spring = 0.0;
    for (bead=0; bead < n_beads; bead++) {
        q2 = 0.0;
        prefactor = 2.0 / (hbar * beta_n);
        s = sin(bead * piovern);
        omega_jsq = prefactor * prefactor * s * s;
        for (j=0; j < n_dof; ++j) {
            nm_kin += 0.5 * momentum_array[bindex(bead, j)] * momentum_array[bindex(bead, j)] / mass_array[j];
            q2 += mass_array[j] * position_array[bindex(bead, j)] * position_array[bindex(bead, j)];
        }
        q2 *= 0.5 * omega_jsq;
        v_spring += q2;
    }
    ftransform(position_array, n_dof, 0);
    ftransform(momentum_array, n_dof, 0);
    
    return energy + nm_kin + v_spring;
}
