//a set of functions to calculate the classical forces and energy
//for qtip4pf water.
//also included are functions to include an additional splined potential
//to allow for "fast" centroid molecular dynamics (CMD) calculations
//using boost to generate the splines
//finally, also needed are additional functions to calculate the objects require
//to run path integral molecular dynamics simulations using the Suzuki-Chin propagator
//which are used to calculate quantum mechanical static properties and to generate the
//potentials used to run the CMD simulations
//these use OMP to parallelise the loops over the different system replicas

#include <vector>
#include <fstream>
#include <iostream>
#include <array>
#include <boost/math/interpolators/cubic_b_spline.hpp>
#include "omp.h"
#include "qtip4pfparams.h"
#include "forces.h"
#include "indexing.h"

//calculate qtip4pf energy and forces

//lennard-jones interaction between oxygen atoms
double lj_oxygen(const std::array<double, 3>& r_ij, double rijsq, double force_array[], int i, int j, int bead) {
    //indexing parameters
    int dim;
    //we need to compute rij^-2, rij^-6 and rij^-12
    double rijpm2, rijpm6, rijpm12;
    //energy and pair force
    double energy, pair_force;
    
    //various distance related properties we need to calculate
    rijpm2 = 1.0 / rijsq;
    rijpm6 = rijpm2 * rijpm2 * rijpm2;
    rijpm12 = rijpm6 * rijpm6;
    //scaling as we are in atomic rather than LJ units
    rijpm6 *= sigmalj6;
    rijpm12 *= sigmalj12;
    
    //calculate energy (to be scaled by 4*epsilon later)
    energy = rijpm12 - rijpm6;
    //we can also use the energy in this form to calculate the pair force (reduces computations)
    pair_force = energy + rijpm12;
    //scale energy
    energy *= 4.0 * epslj;
    //subtract energy at cutoff radius (maintains energy conservation as molecules move across the cutoff radius in the NVE ensemble)
    energy -= e_cut;
    //and scale forces (this approach only requires us to calculate |rij|^2 and NOT |rij| - so we don't need to use sqrt())
    pair_force *= 24.0 * epslj * rijpm2;
    
    for (dim=0; dim < 3; ++dim) {
        //add on the force contributions to the oxygen atoms
        force_array[bindex(bead, 9*i+dim)] += pair_force * r_ij[dim];
        force_array[bindex(bead, 9*j+dim)] -= pair_force * r_ij[dim];
    }
    
    return energy;
}

//intramolecular quartic expansion of morse potential to allow bond flexibility
//here "i" is the index of the molecule in which we are calculating the potential for
//bond_step is used to keep track of which hydrogen we are calculating the potential for
//ie bond_step=1 calculates the bond potential between the oxygen and the "first" labelled hydrogen, and bond_step=2 for the second
double intra_morse_bond(const std::array<double, 3>& r_oh, double roh, double force_array[], int i, int bond_step, int bead) {
    
    int dim;
    double energy, pair_force;
    double rijpm1, rdisp, rdispp2, rdispp3, rdispp4;
    
    //calculate parameters related to the distance needed later
    rijpm1 = 1.0 / roh;
    rdisp = roh - roheq;
    rdispp2 = rdisp * rdisp;
    rdispp3 = rdispp2 * rdisp;
    rdispp4 = rdispp2 * rdispp2;
    
    //calculate contribution to the energy
    energy = ar2*rdispp2 - ar3*rdispp3 + ar4fac*rdispp4;
    energy *= Dr;
    
    //and 1/r dU/dr
    pair_force = 2.0*ar2*rdisp - 3.0*ar3*rdispp2 + 4.0*ar4fac*rdispp3;
    pair_force *= Dr * rijpm1;
    
    //then calculate force contribution on this oxygen and hydrogen
    for (dim=0; dim < 3; ++dim) {
        force_array[bindex(bead, 9*i+dim)] -= pair_force * r_oh[dim];
        force_array[bindex(bead, 9*i+3*bond_step+dim)] += pair_force * r_oh[dim];
    }
    
    return energy;
}

//harmonic potential based on intramolecular H-O-H bond angle
double intra_harm_angle(const std::array<double, 3>& r_oh1, const std::array<double, 3>& r_oh2, double roh1, double roh2, double force_array[], int i, int bead) {
    
    int dim;
    double costheta, costhetasq, theta, thetadisp, dotprod;
    double r1pm1, r2pm1, r1pm2, r2pm2, r1r2pm1;
    double dxdcosx, c2r1, c2r2;
    double energy, pair_f;
    
    //compute H-O-H angle from dot product
    dotprod = r_oh1[0]*r_oh2[0] + r_oh1[1]*r_oh2[1] + r_oh1[2]*r_oh2[2];
    //calculate different inverse distances that are needed and compute some useful products
    r1pm1 = 1.0 / roh1;
    r2pm1 = 1.0 / roh2;
    r1pm2 = r1pm1 * r1pm1;
    r2pm2 = r2pm1 * r2pm1;
    r1r2pm1 = r1pm1 * r2pm1;
    
    //cosine of the bond angle
    costheta = dotprod*r1r2pm1;
    costhetasq = costheta * costheta;
    //various quantities need to compute the derivatives
    dxdcosx = 1.0 / sqrt(1.0 - costhetasq);
    c2r1 = dotprod*r1pm2 * r1r2pm1;
    c2r2 = dotprod*r2pm2 * r1r2pm1;
    //bond angle in radians
    theta = acos(costheta);
    thetadisp = theta - theta_eq;
    
    energy = 0.5 * ktheta * thetadisp * thetadisp;
    //precompute a shared term in all the different derivatives
    pair_f = ktheta * thetadisp * dxdcosx;

    for (dim=0; dim < 3; ++dim) {
        //calculate the forces on each atom in the molecule
        //this is fairly messy, but of course so are the expressions themselves
        //perhaps this can be simplified somewhat, but as long as it works i'm not going to worry about it too much
        force_array[bindex(bead, 9*i+dim)] += pair_f * (r1r2pm1 * (r_oh1[dim] + r_oh2[dim]) - c2r1 * r_oh1[dim] - c2r2 * r_oh2[dim]);
        force_array[bindex(bead, 9*i+dim+3)] -= pair_f * (r1r2pm1 * r_oh2[dim] - c2r1 * r_oh1[dim]);
        force_array[bindex(bead, 9*i+dim+6)] -= pair_f * (r1r2pm1 * r_oh1[dim] - c2r2 * r_oh2[dim]);
    }
    
    return energy;
}

//real space Ewald interaction
double rwald(const std::array<double, 3>& r_ij, double rijsq, double force_array[], const double charge_array[], int i, int j, double ewald_kappa, int bead) {
    
    int dim;
    double energy, pair_force;
    double rij, rijpm1, rijpm2, qiqj;
    double ewald_kappasq = ewald_kappa * ewald_kappa;
    double ewald_prefactor = 2.0 * ewald_kappa / sqrt(pi);
    
    rij = sqrt(rijsq);
    rijpm1 = 1.0 / rij;
    rijpm2 = rijpm1 * rijpm1;
    qiqj = charge_array[i] * charge_array[j];
    
    //currently using the erfc implemented in c++ rather than any expansion
    //i have tested the expansion used in David's old code and this gave very little speed-up, if any
    double erfc_val = erfc(ewald_kappa*rij);
    //if the two atoms are in the same molecule, then we need to account for this
    if (i/3 == j/3) {
        //remove self interaction
        //the case when i==j will be taken care of later, since this term is not included in the double sum
        //so this checks whether i and j are in the same molecule, but where i!=j
        //in this case we need to use -erf() rather than erfc()
        //(Belhadj, Halper, Levy, 1990)
        //-erf(x) = erfc(x) - 1
        erfc_val -= 1.0;
    }
    //compute rwald energy and forces
    energy = qiqj * rijpm1 * erfc_val;
    pair_force = energy + qiqj * ewald_prefactor * exp(-ewald_kappasq*rijsq);
    pair_force *= rijpm2;
    
    for (dim=0; dim < 3; ++dim) {
        force_array[bindex(bead, 3*i+dim)] += pair_force * r_ij[dim];
        force_array[bindex(bead, 3*j+dim)] -= pair_force * r_ij[dim];
    }
    return energy;
}

//kwald interaction
double kwald(const double position_array[], double force_array[], const double charge_array[], int k_max, double ewald_kappa, double k_cut_sq, int bead) {
    
    //indexing variables
    int i, dim, j;
    int kx, ky, kz;
    double ewald_kappasq = ewald_kappa * ewald_kappa;
    double quartalphasqpm1 = 0.25 / ewald_kappasq;
    //variable to hold the reciprocal space component to the energy
    //initialise as zero for now
    double energy = 0.0;
    //bk is a 3D array holding the reciprocal lattice vectors of the simulation cell
    //3D arrays to hold the k_1 \cdot r_i and cos() and sin() of this
    //as well as the components of each reciprocal lattice vector in the sum
    std::array<double, 3> bk, kr, ck, sk, rk;
    //vectors to store the real and imaginary parts of e^(i*k*r) for each reciprocal lattice vector in the sum
    std::vector<std::vector<double> > cos_k(n_dof, std::vector<double>(k_max+1));
    std::vector<std::vector<double> > sin_k(n_dof, std::vector<double>(k_max+1));
    //various things needed to compute the energy and force contribution of each reciprocal lattice vector
    double ksq, exp_fac, sum_real, sum_imag, pair_force;
    double sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8;
    double coscos, cossin, sincos, sinsin, cca, csa, sca, ssa, ccb, csb, scb, ssb;
    double term1, term2, term3, term4, xterm, yterm, zterm;
    std::vector<std::array<double, 4> > xysum(n_atoms);
    std::vector<std::array<double, 8> > xyzsum(n_atoms);
    
    //find prefactors for reciprocal lattice vectors of the entire simulation cell
    for (dim=0; dim < 3; ++dim) {
        bk[dim] = 2.0 * pi / boxl[dim];
    }
    
    //ewald summation, reciprocal space part
    //loop over all atoms
    for (i=0; i < n_atoms; ++i) {
        //(e^ikx)^0 = 1.0
        for (dim=0; dim < 3; ++dim) {
            //stores the real part
            //cos_k[3*i+dim][0] = 1.0;
            cos_k[3*i+dim][0] = 1.0;
            //and the imaginary part
            //sin_k[3*i+dim][0] = 0.0;
            sin_k[3*i+dim][0] = 0.0;
        }
        //calculate reciprocal lattice vectors bx, by, and bz dotted with the position of charge i
        for (dim=0; dim < 3; ++dim) {
            kr[dim] = bk[dim] * position_array[bindex(bead, 3*i+dim)];
        }

        //and compute cos and sin of these
        for (dim=0; dim < 3; ++dim) {
            ck[dim] = cos(kr[dim]);
            sk[dim] = sin(kr[dim]);
        }
        //now, loop over k vectors
        for (j=1; j < k_max + 1; ++j) {
            //find each power by successive complex multiplication
            //cos_kx stores the real part of this power, and sin_kx the imaginary part
            for (dim=0; dim < 3; ++dim) {
                cos_k[3*i+dim][j] = ck[dim] * cos_k[3*i+dim][j-1] - sk[dim] * sin_k[3*i+dim][j-1];
                sin_k[3*i+dim][j] = sk[dim] * cos_k[3*i+dim][j-1] + ck[dim] * sin_k[3*i+dim][j-1];
            }
        }
        for (j=0; j < k_max + 1; ++j) {
            cos_k[3*i][j] *= charge_array[i];
            sin_k[3*i][j] *= charge_array[i];
        }
    }
    
    //now evaluate the sum
    //loop over reciprocal lattice vectors in the x direction
    //only loop over vectors with positive (or zero) x component as we can exploit symmetry later on
    //you can exploit more symmetry for a little extra speed-up but this leads to the code becoming (even more) messy - on my todo list as well for optimising
    for (kx=0; kx < k_max+1; ++kx) {
        //store the lengths of this reciprocal lattice vector as an array for ease of use later on
        //no need for ugly indexing here as only positive kx vectors needed
        rk[0] = bk[0] * kx;
        ksq = rk[0] * rk[0];
        
        //case where |ky|=0 and |kz|=0
        if (ksq > 0.0 && ksq < k_cut_sq) {
            sum_real = 0.0;
            sum_imag = 0.0;
            
            exp_fac = exp(-ksq * quartalphasqpm1) / ksq;
            
            if (kx != 0) {
                exp_fac *= 2.0;
            }
            
            for (i=0; i < n_atoms; ++i) {
                sum_real += cos_k[3*i][kx];
                sum_imag += sin_k[3*i][kx];
            }
            
            //increment energy from this k-vector
            energy += kwald_prefactor * exp_fac * (sum_real * sum_real + sum_imag * sum_imag);
            //compute forces from this k-vector
            for (i=0; i < n_atoms; ++i) {
                pair_force = (sin_k[3*i][kx]*sum_real - cos_k[3*i][kx]*sum_imag) * 2.0 * kwald_prefactor * exp_fac * rk[0];
                force_array[bindex(bead, 3*i)] += pair_force;
            }
        }
        
        //ky=0 and |kz| > 0
        for (kz=1; kz < k_max+1; ++kz) {
            rk[2] = bk[2] * kz;
            ksq = rk[0]*rk[0] + rk[2]*rk[2];
            if (ksq > k_cut_sq) {
                break;
            }
            sum1 = 0.0;
            sum2 = 0.0;
            sum3 = 0.0;
            sum4 = 0.0;
            for (i=0; i < n_atoms; ++i) {
                coscos = cos_k[3*i][kx]*cos_k[3*i+2][kz];
                cossin = cos_k[3*i][kx]*sin_k[3*i+2][kz];
                sincos = sin_k[3*i][kx]*cos_k[3*i+2][kz];
                sinsin = sin_k[3*i][kx]*sin_k[3*i+2][kz];
                xysum[i][0] = coscos - sinsin;
                xysum[i][1] = cossin + sincos;
                xysum[i][2] = coscos + sinsin;
                xysum[i][3] = sincos - cossin;
                sum1 += xysum[i][0];
                sum2 += xysum[i][1];
                sum3 += xysum[i][2];
                sum4 += xysum[i][3];
            }
            
            exp_fac = kwald_prefactor * exp(-ksq * quartalphasqpm1) / ksq;
            
            if (kx != 0) {
                exp_fac *= 2.0;
            }
            
            //increment energy from this k-vector
            energy += exp_fac * (sum1 * sum1 + sum2 * sum2 + sum3 * sum3 + sum4 * sum4);
            xterm = 2.0 * exp_fac * rk[0];
            zterm = 2.0 * exp_fac * rk[2];
            //compute forces from this k-vector
            for (i=0; i < n_atoms; ++i) {
                term1 = sum2 * xysum[i][0] - sum1 * xysum[i][1];
                term2 = sum4 * xysum[i][2] - sum3 * xysum[i][3];
                force_array[bindex(bead, 3*i)] -= xterm * (term1 + term2);
                force_array[bindex(bead, 3*i+2)] -= zterm * (term1 - term2);
            }
        }
        
        //|ky| > 0 and kz = 0
        
        for (ky=1; ky < k_max + 1; ++ky) {
            rk[1] = bk[1] * ky;
            ksq = rk[0]*rk[0] + rk[1]*rk[1];
            
            if (ksq > k_cut_sq) {
                break;
            }
            
            sum1 = 0.0;
            sum2 = 0.0;
            sum3 = 0.0;
            sum4 = 0.0;
            
            for (i=0; i < n_atoms; ++i) {
                coscos = cos_k[3*i][kx] * cos_k[3*i+1][ky];
                cossin = cos_k[3*i][kx] * sin_k[3*i+1][ky];
                sincos = sin_k[3*i][kx] * cos_k[3*i+1][ky];
                sinsin = sin_k[3*i][kx] * sin_k[3*i+1][ky];
                xysum[i][0] = coscos - sinsin;
                xysum[i][1] = cossin + sincos;
                xysum[i][2] = coscos + sinsin;
                xysum[i][3] = sincos - cossin;
                sum1 += xysum[i][0];
                sum2 += xysum[i][1];
                sum3 += xysum[i][2];
                sum4 += xysum[i][3];
            }
            
            exp_fac = kwald_prefactor * exp(-ksq * quartalphasqpm1) / ksq;
            
            if (kx != 0) {
                exp_fac *= 2.0;
            }
            
            //increment energy from this k-vector
            energy += exp_fac * (sum1 * sum1 + sum2 * sum2 + sum3 * sum3 + sum4 * sum4);
            xterm = 2.0 * exp_fac * rk[0];
            yterm = 2.0 * exp_fac * rk[1];
            //compute forces from this k-vector
            for (i=0; i < n_atoms; ++i) {
                term1 = sum2 * xysum[i][0] - sum1 * xysum[i][1];
                term2 = sum4 * xysum[i][2] - sum3 * xysum[i][3];
                force_array[bindex(bead, 3*i)] -= xterm * (term1 + term2);
                force_array[bindex(bead, 3*i+1)] -= yterm * (term1 - term2);
            }
            
            //|ky| > 0 and |kz| > 0
            for (kz=1; kz < k_max + 1; ++kz) {
                rk[2] = bk[2] * kz;
                ksq = rk[0]*rk[0] + rk[1]*rk[1] + rk[2]*rk[2];
                
                if (ksq > k_cut_sq) {
                    break;
                }
                
                sum1 = 0.0;
                sum2 = 0.0;
                sum3 = 0.0;
                sum4 = 0.0;
                sum5 = 0.0;
                sum6 = 0.0;
                sum7 = 0.0;
                sum8 = 0.0;
                
                for (i=0; i < n_atoms; ++i) {
                    cca = xysum[i][0] * cos_k[3*i+2][kz];
                    csa = xysum[i][0] * sin_k[3*i+2][kz];
                    sca = xysum[i][1] * cos_k[3*i+2][kz];
                    ssa = xysum[i][1] * sin_k[3*i+2][kz];
                    ccb = xysum[i][2] * cos_k[3*i+2][kz];
                    csb = xysum[i][2] * sin_k[3*i+2][kz];
                    scb = xysum[i][3] * cos_k[3*i+2][kz];
                    ssb = xysum[i][3] * sin_k[3*i+2][kz];
                    xyzsum[i][0] = cca - ssa;
                    xyzsum[i][1] = sca + csa;
                    xyzsum[i][2] = cca + ssa;
                    xyzsum[i][3] = sca - csa;
                    xyzsum[i][4] = ccb - ssb;
                    xyzsum[i][5] = scb + csb;
                    xyzsum[i][6] = ccb + ssb;
                    xyzsum[i][7] = scb - csb;
                    
                    sum1 += xyzsum[i][0];
                    sum2 += xyzsum[i][1];
                    sum3 += xyzsum[i][2];
                    sum4 += xyzsum[i][3];
                    sum5 += xyzsum[i][4];
                    sum6 += xyzsum[i][5];
                    sum7 += xyzsum[i][6];
                    sum8 += xyzsum[i][7];
                }
                
                exp_fac = kwald_prefactor * exp(-ksq * quartalphasqpm1) / ksq;
                
                if (kx != 0) {
                    exp_fac *= 2.0;
                }
                
                //increment energy from this k-vector
                energy += exp_fac * (sum1*sum1 + sum2*sum2 + sum3*sum3 + sum4*sum4 + sum5*sum5 + sum6*sum6 + sum7*sum7 + sum8*sum8);
                
                xterm = 2.0 * exp_fac * rk[0];
                yterm = 2.0 * exp_fac * rk[1];
                zterm = 2.0 * exp_fac * rk[2];
                //compute forces from this k-vector
                for (i=0; i < n_atoms; ++i) {
                    term1 = sum2 * xyzsum[i][0] - sum1 * xyzsum[i][1];
                    term2 = sum4 * xyzsum[i][2] - sum3 * xyzsum[i][3];
                    term3 = sum6 * xyzsum[i][4] - sum5 * xyzsum[i][5];
                    term4 = sum8 * xyzsum[i][6] - sum7 * xyzsum[i][7];
                    force_array[bindex(bead, 3*i)] -= xterm * (term1 + term2 + term3 + term4);
                    force_array[bindex(bead, 3*i+1)] -= yterm * (term1 + term2 - term3 - term4);
                    force_array[bindex(bead, 3*i+2)] -= zterm * (term1 - term2 + term3 - term4);
                }
            }
        }
    }
    
    return energy;
}

//a function to find the forces and energy from a generic splined potential
//"spline" will be one of the O-O, O-H or H-H splines
double get_spline(const std::array<double, 3> r_ij, const double rij, double force_array[], const int i, const int j, const boost::math::cubic_b_spline<double>& spline, int bead) {
    
    int dim;
    double rijpm1 = 1.0 / rij;
    double energy = spline(rij);
    double fp = spline.prime(rij);
    fp *= rijpm1;
    
    for (dim=0; dim < 3; ++dim) {
        force_array[bindex(bead, 3*i+dim)] -= fp * r_ij[dim];
        force_array[bindex(bead, 3*j+dim)] += fp * r_ij[dim];
    }
    
    return energy;
}

double get_forces_nospline(const double position_array[], double force_array[], const double charge_array[], const int point[], const int list[], int bead) {
    
    int i, j, dim, neighbour_i;
    //array to hold displacement vectors
    std::array<double, 3> r_ij, r_oh1, r_oh2;
    //squared distance
    double rijsq, roh1sq, roh1, roh2sq, roh2;
    double* ewald_force = new double[n_beads*n_dof];
    double* charge_pos = new double[n_beads*n_dof];
    double energy;
    //initialise energy for this replica as zero
    energy = 0.0;
    //set all forces in each replica to zero
    for (i=0; i < n_atoms; ++i) {
        for (dim=0; dim < 3; ++dim) {
            force_array[bindex(bead, 3*i+dim)] = 0.0;
        }
    }
    //loop over all oxygen atoms in the system
    //only go up to penultimate oxygen as point[n_mols] is undefined memory
    //where we can miss out the last oxygen due to Newton's third law
    for (i=0; i < n_mols-1; ++i) {
        //calculate potential energy contribution here
        //loop over neighbouring oxygens
        for (neighbour_i=point[pindex(bead, i)]; neighbour_i < point[pindex(bead, i+1)]; ++neighbour_i) {
            //find index of this neighbour
            j = list[lindex(bead, neighbour_i)];
            //compute distance vector between the two oxygens
            for (dim=0; dim < 3; ++dim) {
                //only calculate using oxygen positions (why there's a factor of 9 instead of 3)
                r_ij[dim] = position_array[bindex(bead, 9*i+dim)] - position_array[bindex(bead, 9*j+dim)];
            }
            
            //minimum image separation in cell units
            for (dim=0; dim < 3; ++dim) {
                r_ij[dim] -= boxl[dim] * std::round(r_ij[dim] / boxl[dim]);
            }
            
            //now compute squared distance
            rijsq = r_ij[0]*r_ij[0] + r_ij[1]*r_ij[1] + r_ij[2]*r_ij[2];
            
            //check if oxygens are close enough to be interacting
            if (rijsq < r_cut_sq) {
                //and if they are, compute the energy and forces between them
                energy += lj_oxygen(r_ij, rijsq, force_array, i, j, bead);
            }
        }
    }
    
    //intramolecular component
    //now we need to loop over all molecules to include the intramolecular components in the final molecule
    for (i=0; i < n_mols; ++i) {
        //calculate distance vector to first bonded hydrogen
        for (dim=0; dim < 3; ++dim) {
            r_oh1[dim] = position_array[bindex(bead, 9*i+dim)] - position_array[bindex(bead, 9*i+dim+3)];
        }
        
        //minimum image separation in cell units
        for (dim=0; dim < 3; ++dim) {
            r_oh1[dim] -= boxl[dim] * std::round(r_oh1[dim] / boxl[dim]);
        }
        
        //squared distance and distance to this hydrogen
        roh1sq = r_oh1[0]*r_oh1[0] + r_oh1[1]*r_oh1[1] + r_oh1[2]*r_oh1[2];
        roh1 = sqrt(roh1sq);
        //bond energy and forces between oxygen and first hydrogen
        energy += intra_morse_bond(r_oh1, roh1, force_array, i, 1, bead);
        
        //distance vector to second hydrogen in this molecule
        for (dim=0; dim < 3; ++dim) {
            r_oh2[dim] = position_array[bindex(bead, 9*i+dim)] - position_array[bindex(bead, 9*i+dim+6)];
        }
        
        //minimum image separation in cell units
        for (dim=0; dim < 3; ++dim) {
            r_oh2[dim] -= boxl[dim] * std::round(r_oh2[dim] / boxl[dim]);
        }
        
        //squared distance and distance to the second hydrogen
        roh2sq = r_oh2[0]*r_oh2[0] + r_oh2[1]*r_oh2[1] + r_oh2[2]*r_oh2[2];
        roh2 = sqrt(roh2sq);
        //bond energy and forces between oxygen and second hydrogen
        energy += intra_morse_bond(r_oh2, roh2, force_array, i, 2, bead);
        
        //and energy and forces from angular term
        energy += intra_harm_angle(r_oh1, r_oh2, roh1, roh2, force_array, i, bead);
    }
    
    //calculate point charge contribution to the energies and forces via Ewald summation
    //calculate position of the m-sites
    for (i=0; i < n_mols; ++i) {
        for (dim=0; dim < 3; ++dim) {
            charge_pos[bindex(bead, 9*i+dim)] = gamma_vs*position_array[bindex(bead, 9*i+dim)] + 0.5 * (1.0 - gamma_vs) * (position_array[bindex(bead, 9*i+dim+3)] + position_array[bindex(bead, 9*i+dim+6)]);
            charge_pos[bindex(bead, 9*i+dim+3)] = position_array[bindex(bead, 9*i+dim+3)];
            charge_pos[bindex(bead, 9*i+dim+6)] = position_array[bindex(bead, 9*i+dim+6)];
        }
    }
    
    //we need a new force_array to keep track of the ewald forces, as the forces on the virtual sites ARE NOT the same as the forces on the oxygen atoms
    for (i=0; i < n_atoms; ++i) {
        for (dim=0; dim < 3; ++dim) {
            ewald_force[bindex(bead, 3*i+dim)] = 0.0;
        }
    }
    
    //parameters for ewald calculation
    double ewald_kappa = 2.6 * std::pow(rho/vol, 1.0/6.0);
    double rwald_cut = acc / ewald_kappa;
    if (rwald_cut > 0.5 * min_box) {
        rwald_cut = 0.5 * min_box;
        ewald_kappa = acc / rwald_cut;
    }
    double rwald_cut_sq = rwald_cut * rwald_cut;
    double k_cut = 2.0 * acc * ewald_kappa;
    double k_cut_sq = k_cut * k_cut;
    int k_max = static_cast <int> (k_cut * max_box / (2.0 * pi));
    double ewald_prefactor = 2.0 * ewald_kappa / sqrt(pi);
    
    //ewald summation, real space part
    //cutoff for real space interactions currently set as half the smallest box length
    //loop over all atoms in the system
    for (i=0; i < n_atoms-1; ++i) {
        //loop over all other atoms
        for (j=i+1; j < n_atoms; ++j) {
            //calculate distance between each atom pair
            for (dim=0; dim < 3; ++dim) {
                r_ij[dim] = charge_pos[bindex(bead, 3*i+dim)] - charge_pos[bindex(bead, 3*j+dim)];
            }
            
            //minimum image separation in cell units
            for (dim=0; dim < 3; ++dim) {
                r_ij[dim] -= boxl[dim] * std::round(r_ij[dim] / boxl[dim]);
            }
            
            rijsq = r_ij[0]*r_ij[0] + r_ij[1]*r_ij[1] + r_ij[2]*r_ij[2];
            
            if (rijsq < rwald_cut_sq) {
                energy += rwald(r_ij, rijsq, ewald_force, charge_array, i, j, ewald_kappa, bead);
            }
        }
    }
    
    //calculate reciprocal space charge interaction energy
    energy += kwald(charge_pos, ewald_force, charge_array, k_max, ewald_kappa, k_cut_sq, bead);
    //convert ewald forces on virtual site to forces on oxygen/hydrogen via the chain rule
    for (i=0; i < n_mols; ++i) {
        for (dim=0; dim < 3; ++dim) {
            force_array[bindex(bead, 9*i+dim)] += gamma_vs * ewald_force[bindex(bead, 9*i+dim)];
            force_array[bindex(bead, 9*i+dim+3)] += ewald_force[bindex(bead, 9*i+dim+3)] + 0.5 * (1.0 - gamma_vs) * ewald_force[bindex(bead, 9*i+dim)];
            force_array[bindex(bead, 9*i+dim+6)] += ewald_force[bindex(bead, 9*i+dim+6)] + 0.5 * (1.0 - gamma_vs) * ewald_force[bindex(bead, 9*i+dim)];
        }
    }
    
    //subtract self term from energy
    for (i=0; i < n_atoms; ++i) {
        energy -= 0.5 * ewald_prefactor * charge_array[i] * charge_array[i];
    }
    
    delete[] charge_pos;
    delete[] ewald_force;
    
    return energy;
}

//calculates forces on the beads from the potential for a configuration
//best to have this as its own function as we will need to call it multiple times each step to find the Suzuki-Chin forces
double get_forces_wspline(const double position_array[], double force_array[], const double charge_array[], const int point[], const int list[], const boost::math::cubic_b_spline<double>& oo_spline, const boost::math::cubic_b_spline<double>& oh_spline, const boost::math::cubic_b_spline<double>& hh_spline, const double oo_spline_cutsq, const double oh_spline_cutsq, const double hh_spline_cutsq, int bead) {
    
    
    double energy = get_forces_nospline(position_array, force_array, charge_array, point, list, bead);
    
    int i, j, dim;
    std::array<double, 3> r_ij;
    double rijsq, rij;
    
    //if we're including the splined potential, add in a cubic spline between all O-O, O-H and H-H pairs
    for (i=0; i < n_atoms; ++i) {
        for (j=i+1; j < n_atoms; ++j) {
            //compute distance vector
            for (dim=0; dim < 3; ++dim) {
                r_ij[dim] = position_array[bindex(bead, 3*i+dim)] - position_array[bindex(bead, 3*j+dim)];
            }
            
            //minimum image separation in cell units
            for (dim=0; dim < 3; ++dim) {
                r_ij[dim] -= boxl[dim] * std::round(r_ij[dim] / boxl[dim]);
            }
            rijsq = r_ij[0]*r_ij[0] + r_ij[1]*r_ij[1] + r_ij[2]*r_ij[2];
            //check which type of spline is needed
            if (i%3==0) {
                if (j%3==0) {
                    if (rijsq < oo_spline_cutsq) {
                        rij = sqrt(rijsq);
                        energy += get_spline(r_ij, rij, force_array, i, j, oo_spline, bead);
                    }
                }
                
                else {
                    if (rijsq < oh_spline_cutsq) {
                        rij = sqrt(rijsq);
                        energy += get_spline(r_ij, rij, force_array, i, j, oh_spline, bead);
                    }
                }
            }
            
            else {
                if (j%3==0) {
                    if (rijsq < oh_spline_cutsq) {
                        rij = sqrt(rijsq);
                        energy += get_spline(r_ij, rij, force_array, i, j, oh_spline, bead);
                    }
                }
                else {
                    if (rijsq < hh_spline_cutsq) {
                        rij = sqrt(rijsq);
                        energy += get_spline(r_ij, rij, force_array, i, j, hh_spline, bead);
                    }
                }
            }
        }
    }
    
    return energy;
}

//function to calculate the Suzuki-Chin forces on the atoms
//force_array will store the Suzuki-Chin forces that are used to propagate the dynamics
//whilst fphys will store the "ordinary" forces, which are needed to calculate the centroid-virial kinetic energy later on
//the energy returned by this function gives the value associated with the potential part of the Suzuki-Chin conserved quantity
//whereas sc_pot stores the value used in the potential energy estimator
double force_eval_sc_nospline(const double position_array[], double force_array[], const double charge_array[], const double mass_array[], const int point[], const int list[], double& sc_pot, double fphys[]) {
    
    int i, bead;
    double fsq, del_term, delta, sc_scale, fscale, dq, eplus, eminus;
    double* qplus = new double[n_beads*n_dof];
    double* qminus = new double[n_beads*n_dof];
    double* fplus = new double[n_beads*n_dof];
    double* fminus = new double[n_beads*n_dof];
    
    //set initial potential energy to zero to be incremented
    //tot_energy is the potential contribution to the Suzuki-Chin conserved quantity
    double tot_energy = 0.0;
    //energy will keep track of the potential energy of each replica
    double energy;
    //sc_pot is the Suzuki-Chin potential energy estimator
    double vop = 0.0;
    
    //loop over all the replicas in the extended system
#pragma omp parallel for default(shared) private(i, energy, fsq, del_term, delta, sc_scale, fscale, dq, eplus, eminus) reduction(+:tot_energy, vop)
    for (bead=0; bead < n_beads; ++bead) {
        //compute classical forces and potential for this replica
        //fphys keeps track of the classical forces (-dV/dq) as this is needed for the centroid virial kinetic energy estimator
        //whereas force_array will be used to update the momenta in the update_trajectory step
        energy = get_forces_nospline(position_array, fphys, charge_array, point, list, bead);
        
        //contribution to Suzuki-Chin conserved quantity from the potential on the even beads
        if (bead % 2 == 0) {
            tot_energy += 2.0 * energy;
            //we can also calculate the Suzuki-Chin forces on the even beads here as well as no other information is needed
            for (i=0; i < n_dof; ++i) {
                force_array[bindex(bead, i)] = 2.0 * fphys[bindex(bead, i)] / 3.0;
            }
            //only calculate the Suzuki-Chin potential energy on the even beads
            vop += energy;
        }
        //contribution to Suzuki-Chin conserved quantity from the potential on the odd beads
        else {
            tot_energy += 4.0 * energy;
        
            //accumulate sum of squared forces (weighted by mass) for the odd beads (Suzuki-Chin term w/ alpha=0 so even beads don't contribute)
            fsq = 0.0;
            for (i=0; i < n_dof; ++i) {
                fsq += fphys[bindex(bead, i)] * fphys[bindex(bead, i)] / mass_array[i];
            }
            fsq /= 3.0 * omega_psq;
            //add on this correction term to the Suzuki-Chin conserved quantity (to be scaled by 1/3rd later)
            tot_energy += fsq;
            
            //calculate finite difference scaling parameter for this odd indexed replica
            delta = 0.0;
            for (i=0; i < n_dof; ++i) {
                del_term = fphys[bindex(bead, i)] / mass_array[i];
                delta += del_term * del_term;
            }
            delta /= n_atoms * n_beads;
            delta = 1.0 / sqrt(delta);
            sc_scale = epsilonsc * delta;
            fscale = 1.0 / (9.0 * omega_psq * sc_scale);
            
            //set up positions projections for approximating the Hessian via finite difference
            for (i=0; i < n_dof; ++i) {
                dq = sc_scale * fphys[bindex(bead, i)] / mass_array[i];
                qplus[bindex(bead, i)] = position_array[bindex(bead, i)] + dq;
                qminus[bindex(bead, i)] = position_array[bindex(bead, i)] - dq;
            }
            
            //fill in forces for displaced positions
            //assumes that the neighbour list for the displaced systems are the same as the neighbour list for the original system
            //which should be true, since any displacement is small - and so any error introduced will be negligible (and most of the time actually zero)
            eplus = get_forces_nospline(qplus, fplus, charge_array, point, list, bead);
            eminus = get_forces_nospline(qminus, fminus, charge_array, point, list, bead);
            
            //now we can calculate the Suzuki-Chin forces on the odd beads as well
            for (i=0; i < n_dof; ++i) {
                force_array[bindex(bead, i)] = 4.0 * fphys[bindex(bead, i)] / 3.0;
                //NB. there appears to be an error in the HO Path Integrals made easy paper - this term has had the sign changed as a result
                force_array[bindex(bead, i)] -= fscale * (fplus[bindex(bead, i)] - fminus[bindex(bead, i)]);
            }
        }
        
    }
    //return the total potential energy of all the replicas
    tot_energy /= 3.0;
    sc_pot = 2.0 * vop / n_beads;
    
    delete[] qplus;
    delete[] qminus;
    delete[] fplus;
    delete[] fminus;
    return tot_energy;
    
}

double force_eval_sc_wspline(const double position_array[], double force_array[], const double charge_array[], const double mass_array[], const int point[], const int list[], const boost::math::cubic_b_spline<double>& oo_spline, const boost::math::cubic_b_spline<double>& oh_spline, const boost::math::cubic_b_spline<double>& hh_spline, const double oo_spline_cutsq, const double oh_spline_cutsq, const double hh_spline_cutsq, double& sc_pot, double fphys[]) {
    
    int i, bead;
    double fsq, del_term, delta, sc_scale, fscale, dq, eplus, eminus;
    double* qplus = new double[n_beads*n_dof];
    double* qminus = new double[n_beads*n_dof];
    double* fplus = new double[n_beads*n_dof];
    double* fminus = new double[n_beads*n_dof];
    
    //set initial potential energy to zero to be incremented
    //tot_energy is the potential contribution to the Suzuki-Chin conserved quantity
    double tot_energy = 0.0;
    //energy will keep track of the potential energy of each replica
    double energy;
    //sc_pot is the Suzuki-Chin potential energy estimator
    double vop = 0.0;
    
    //loop over all the replicas in the extended system
#pragma omp parallel for default(shared) private(i, energy, fsq, del_term, delta, sc_scale, fscale, dq, eplus, eminus) reduction(+:tot_energy, vop)
    for (bead=0; bead < n_beads; ++bead) {
        //compute classical forces and potential for this replica
        //fphys keeps track of the classical forces (-dV/dq) as this is needed for the centroid virial kinetic energy estimator
        //whereas force_array will be used to update the momenta in the update_trajectory step
        energy = get_forces_wspline(position_array, fphys, charge_array, point, list, oo_spline, oh_spline, hh_spline, oo_spline_cutsq, oh_spline_cutsq, hh_spline_cutsq, bead);
        
        //contribution to Suzuki-Chin conserved quantity from the potential on the even beads
        if (bead % 2 == 0) {
            tot_energy += 2.0 * energy;
            //we can also calculate the Suzuki-Chin forces on the even beads here as well as no other information is needed
            for (i=0; i < n_dof; ++i) {
                force_array[bindex(bead, i)] = 2.0 * fphys[bindex(bead, i)] / 3.0;
            }
            //only calculate the Suzuki-Chin potential energy on the even beads
            vop += energy;
        }
        //contribution to Suzuki-Chin conserved quantity from the potential on the odd beads
        else {
            tot_energy += 4.0 * energy;
        
            //accumulate sum of squared forces (weighted by mass) for the odd beads (Suzuki-Chin term w/ alpha=0 so even beads don't contribute)
            fsq = 0.0;
            for (i=0; i < n_dof; ++i) {
                fsq += fphys[bindex(bead, i)] * fphys[bindex(bead, i)] / mass_array[i];
            }
            fsq /= 3.0 * omega_psq;
            //add on this correction term to the Suzuki-Chin conserved quantity (to be scaled by 1/3rd later)
            tot_energy += fsq;
            
            //calculate finite difference scaling parameter for this odd indexed replica
            delta = 0.0;
            for (i=0; i < n_dof; ++i) {
                del_term = fphys[bindex(bead, i)] / mass_array[i];
                delta += del_term * del_term;
            }
            delta /= n_atoms * n_beads;
            delta = 1.0 / sqrt(delta);
            sc_scale = epsilonsc * delta;
            fscale = 1.0 / (9.0 * omega_psq * sc_scale);
            
            //set up positions projections for approximating the Hessian via finite difference
            for (i=0; i < n_dof; ++i) {
                dq = sc_scale * fphys[bindex(bead, i)] / mass_array[i];
                qplus[bindex(bead, i)] = position_array[bindex(bead, i)] + dq;
                qminus[bindex(bead, i)] = position_array[bindex(bead, i)] - dq;
            }
            
            //fill in forces for displaced positions
            //assumes that the neighbour list for the displaced systems are the same as the neighbour list for the original system
            //which should be true, since any displacement is small - and so any error introduced will be negligible (and most of the time actually zero)
            eplus = get_forces_wspline(qplus, fplus, charge_array, point, list, oo_spline, oh_spline, hh_spline, oo_spline_cutsq, oh_spline_cutsq, hh_spline_cutsq, bead);
            eminus = get_forces_wspline(qminus, fminus, charge_array, point, list, oo_spline, oh_spline, hh_spline, oo_spline_cutsq, oh_spline_cutsq, hh_spline_cutsq, bead);
            
            //now we can calculate the Suzuki-Chin forces on the odd beads as well
            for (i=0; i < n_dof; ++i) {
                force_array[bindex(bead, i)] = 4.0 * fphys[bindex(bead, i)] / 3.0;
                //NB. there appears to be an error in the HO Path Integrals made easy paper - this term has had the sign changed as a result
                force_array[bindex(bead, i)] -= fscale * (fplus[bindex(bead, i)] - fminus[bindex(bead, i)]);
            }
        }
        
    }
    //return the total potential energy of all the replicas
    tot_energy /= 3.0;
    sc_pot = 2.0 * vop / n_beads;
    
    delete[] qplus;
    delete[] qminus;
    delete[] fplus;
    delete[] fminus;
    return tot_energy;
    
}

