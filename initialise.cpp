//a set of functions to either initialise water molecules on a face centred cubic lattice
//for liquid simulations (after equilibration and melting)
//or from an orthorhombic lattice, constrained to have zero dipole moment and to satisfy
//the "ice rules", to allow for simulations of ice
//also included is functionality to start from a configuration provided in a separate file
//helpful to continue simulations if required

#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include <random>
#include "initialise.h"
#include "qtip4pfparams.h"
#include "indexing.h"

//function to initialise starting atomic positions for water simulation
//generate water molecules on a starting fcc lattice
double* initialise_fcc() {
    //array to hold positions
    //3*n_atoms dimensional array,
    // with structure O_x, O_y, O_z, H1_x, H1_y, H1_z, H2_x, H2_y, H2_z
    double* position_array = new double[n_beads*n_dof];
    //indexing parameters for loops
    int i, j, k, ci, bead;
    double xdiff, zdiff;
    
    xdiff = roheq * cos(theta_eq/2.0) / sqrt(2.0);
    zdiff = roheq * sin(theta_eq / 2.0);
    
    //loop over fcc unit cells
    for (i=0; i < n_cell_x; ++i) {
        for (j=0; j < n_cell_y; ++j) {
            for (k=0; k < n_cell_z; ++k) {
                ci = 36 * (k + n_cell_z*j + n_cell_z*n_cell_y*i);
                //fill in oxygen positions
                //first oxygen
                position_array[ci] = (i+0.25)*a_cell;
                position_array[ci+1] = (j+0.25)*a_cell;
                position_array[ci+2] = (k+0.25)*a_cell;
                //hydrogens bonded to first oxygen
                position_array[ci+3] = position_array[ci] + xdiff;
                position_array[ci+4] = position_array[ci+1] + xdiff;
                position_array[ci+5] = position_array[ci+2] + zdiff;
                position_array[ci+6] = position_array[ci] + xdiff;
                position_array[ci+7] = position_array[ci+1] + xdiff;
                position_array[ci+8] = position_array[ci+2] - zdiff;
                
                //second oxygen
                position_array[ci+9] = (i+0.25)*a_cell;
                position_array[ci+10] = (j+0.75)*a_cell;
                position_array[ci+11] = (k+0.75)*a_cell;
                //hydrogens bonded to second oxygen
                position_array[ci+12] = position_array[ci+9] + xdiff;
                position_array[ci+13] = position_array[ci+10] + xdiff;
                position_array[ci+14] = position_array[ci+11] + zdiff;
                position_array[ci+15] = position_array[ci+9] + xdiff;
                position_array[ci+16] = position_array[ci+10] + xdiff;
                position_array[ci+17] = position_array[ci+11] - zdiff;
                
                //third oxygen
                position_array[ci+18] = (i+0.75)*a_cell;
                position_array[ci+19] = (j+0.25)*a_cell;
                position_array[ci+20] = (k+0.75)*a_cell;
                //hydrogens bonded to third oxygen
                position_array[ci+21] = position_array[ci+18] + xdiff;
                position_array[ci+22] = position_array[ci+19] + xdiff;
                position_array[ci+23] = position_array[ci+20] + zdiff;
                position_array[ci+24] = position_array[ci+18] + xdiff;
                position_array[ci+25] = position_array[ci+19] + xdiff;
                position_array[ci+26] = position_array[ci+20] - zdiff;
                
                //fourth oxygen
                position_array[ci+27] = (i+0.75)*a_cell;
                position_array[ci+28] = (j+0.75)*a_cell;
                position_array[ci+29] = (k+0.25)*a_cell;
                //hydrogens bonded to fourth oxygen
                position_array[ci+30] = position_array[ci+27] + xdiff;
                position_array[ci+31] = position_array[ci+28] + xdiff;
                position_array[ci+32] = position_array[ci+29] + zdiff;
                position_array[ci+33] = position_array[ci+27] + xdiff;
                position_array[ci+34] = position_array[ci+28] + xdiff;
                position_array[ci+35] = position_array[ci+29] - zdiff;
            }
        }
    }
    
    //set initial positions in all other replicas to those in the first replica
    for (bead=1; bead < n_beads; ++bead) {
        for (i=0; i < n_dof; ++i) {
            position_array[bindex(bead, i)] = position_array[i];
        }
    }
    
    //return the initial positions
    return position_array;
}

//function to read in positions from some other file
//if you want to read in positions from another file
//assumes the input file will be laid out like
//O1x O1y O1z
//H1x H1y H1z
//H2x H2y H2z
//where H1 and H2 are the hydrogens bonded to O1
double* initialise_fromfile(std::string filename) {
    double* position_array = new double[n_beads*n_dof];
    int i, dim, bead;
    int count = 0;
    
    std::ifstream infile;
    infile.open(filename);
    while (infile.is_open()) {
        std::string line;
        while (std::getline(infile, line)) {
            if (line != "") {
                std::stringstream ss(line);
                for (dim=0; dim < 3; ++dim) {
                    ss >> position_array[3*count+dim];
                }
                count += 1;
            }
        }
        infile.close();
    }
    
    for (bead=1; bead < n_beads; ++bead) {
        for (i=0; i < n_dof; ++i) {
            position_array[bindex(bead, i)] = position_array[i];
        }
    }
    
    return position_array;
}

std::array<double, 3> get_dipole_moment(const double position_array[]) {
    std::array<double, 3> dipole_moment;
    std::array<double, 3> r_oh1, r_oh2, r_bisect;
    double rbsq, rb;
    int i, dim;
    
    for (dim=0; dim < 3; ++dim) {
        dipole_moment[dim] = 0.0;
    }
    
    //loop over all molecules in the simulation
    for (i=0; i < n_mols; ++i) {
        for (dim=0; dim < 3; ++dim) {
            //find separation vector between oxygen and first hydrogen in this molecule
            r_oh1[dim] = position_array[9*i+dim+3] - position_array[9*i+dim];
            //periodic boundary conditions
            //this shouldn't actually be needed, I don't think
            //but just in case a O-H goes over one of the box sides
            if (r_oh1[dim] <= -0.5 * boxl[dim]) {
                r_oh1[dim] += boxl[dim];
            }
            else if (r_oh1[dim] > 0.5 * boxl[dim]) {
                r_oh1[dim] -= boxl[dim];
            }
            
            //find separation vector between oxygen and second hydrogen in this molecule
            r_oh2[dim] = position_array[9*i+dim+6] - position_array[9*i+dim];
            
            if (r_oh2[dim] <= -0.5 * boxl[dim]) {
                r_oh2[dim] += boxl[dim];
            }
            else if (r_oh2[dim] > 0.5 * boxl[dim]) {
                r_oh2[dim] -= boxl[dim];
            }
            
            //find vector that bisects the H-O-H angle
            r_bisect[dim] = r_oh1[dim] + r_oh2[dim];
        }
        //find length of current bisection vector
        rbsq = r_bisect[0]*r_bisect[0] + r_bisect[1]*r_bisect[1] + r_bisect[2]*r_bisect[2];
        rb = sqrt(rbsq);
        
        //and normalise bisection vector and add this to overall dipole moment vector
        for (dim=0; dim < 3; ++dim) {
            r_bisect[dim] /= rb;
            dipole_moment[dim] += r_bisect[dim];
        }
    }
    
    return dipole_moment;
}

double* initialise_oxygens() {
    //array to hold positions
    //3*n_atoms dimensional array,
    // with structure O_x, O_y, O_z, H1_x, H1_y, H1_z, H2_x, H2_y, H2_z
    //this function will only fill in the oxygen positions
    double* position_array = new double[n_beads*n_dof];
    //indexing parameters for loops
    int i, j, k, index, bead;
    
    for (bead=0; bead < n_beads; ++bead) {
        for (i=0; i < n_dof; ++i) {
            position_array[bindex(bead, i)] = 0.0;
        }
    }
    
    //set up oxygen positions on two interleaved hexagonal lattices in cell units
    for (i=0; i < n_cell_x; ++i) {
        for (j=0; j < n_cell_y; ++j) {
            for (k=0; k < n_cell_z; ++k) {
                index = 72 * k + 72 * n_cell_z * j + 72 * n_cell_y * n_cell_z * i;
                //position of first oxygen in the cell
                position_array[index] = (i + 0.5) * a_cell;
                position_array[index+1] = (j + one_sixth) * b_cell;
                position_array[index+2] = (k+zo) * c_cell;
                
                //position of second oxygen in the cell
                position_array[index+9] = (i+0.5) * a_cell;
                position_array[index+10] = (j + one_sixth) * b_cell;
                position_array[index+11] = (k + 0.5 - zo) * c_cell;
                
                //position of the third oxygen in the cell
                position_array[index+18] = i * a_cell;
                position_array[index+19] = (j + one_third) * b_cell;
                position_array[index+20] = (k + 0.5 + zo) * c_cell;
                
                //position of the fourth oxygen in the cell
                position_array[index+27] = i * a_cell;
                position_array[index+28] = (j + one_third) * b_cell;
                position_array[index+29] = (k + 1.0 - zo) * c_cell;
                
                //position of the fifth oxygen in the cell
                position_array[index+36] = i * a_cell;
                position_array[index+37] = (j + two_thirds) * b_cell;
                position_array[index+38] = (k + zo) * c_cell;
                
                //position of the sixth oxygen in the cell
                position_array[index+45] = i * a_cell;
                position_array[index+46] = (j + two_thirds) * b_cell;
                position_array[index+47] = (k + 0.5 - zo) * c_cell;
                
                //position of the seventh oxygen in the cell
                position_array[index+54] = (i + 0.5) * a_cell;
                position_array[index+55] = (j + five_sixths) * b_cell;
                position_array[index+56] = (k + 0.5 + zo) * c_cell;
               
                //position of the final oxygen in the orthorhombic unit cell
                position_array[index+63] = (i + 0.5) * a_cell;
                position_array[index+64] = (j + five_sixths) * b_cell;
                position_array[index+65] = (k + 1.0 - zo) * c_cell;
            }
        }
    }
    return position_array;
}

//a function to get a list of distinct O-O nearest neighbours - of which there are 2*n_mols
//will return a 2*n_mols by 2 2D array where ox_connect[i][0] gives the first oxygen index in the i^th nearest neighbour pair and ox_connect[i][1] the second oxygen index
std::vector<std::array<int, 2> > get_oxygen_connectivity() {
    std::vector<std::array<int, 2> > ox_connect(2*n_mols);
    int i, j, k, cell_index, bond_index, atom_index, next_cell_index;
    
    //loop over all cells
    for (i=0; i < n_cell_x; ++i) {
        for (j=0; j < n_cell_y; ++j) {
            for (k=0; k < n_cell_z; ++k) {
                cell_index = k + n_cell_z * j + n_cell_z * n_cell_y * i;
                bond_index = 16 * cell_index;
                atom_index = 8 * cell_index;
                //oxygen 0 in this cell is next to oxygen 1 in the same cell
                ox_connect[bond_index][0] = atom_index;
                ox_connect[bond_index][1] = atom_index+1;
                //oxygen 0 in this cell is next to oxygen 3 in the cell below
                next_cell_index = (k-1+n_cell_z) % n_cell_z + n_cell_z * j + n_cell_z * n_cell_y * i;
                ox_connect[bond_index+1][0] = atom_index;
                ox_connect[bond_index+1][1] = 8*next_cell_index + 3;
                //oxygen 0 in this cell is also next to oxygen 3 in the cell below and one along x
                next_cell_index = (k-1+n_cell_z)%n_cell_z + n_cell_z*j + n_cell_z*n_cell_y*((i+1)%n_cell_x);
                ox_connect[bond_index+2][0] = atom_index;
                ox_connect[bond_index+2][1] = 8 * next_cell_index + 3;
                //oxygen 0 in this cell is next to oxygen 7 in the cell below and one along -y
                next_cell_index = (k-1+n_cell_z)%n_cell_z + n_cell_z*((j-1+n_cell_y)%n_cell_y) + n_cell_z*n_cell_y*i;
                ox_connect[bond_index+3][0] = atom_index;
                ox_connect[bond_index+3][1] = 8*next_cell_index + 7;
                //oxygen 1 in this cell is next to oxygen 2 in the same cell
                ox_connect[bond_index+4][0] = atom_index+1;
                ox_connect[bond_index+4][1] = atom_index+2;
                //oxygen 1 in this cell is also next to oxygen 2 in the cell one along x
                next_cell_index = k + n_cell_z*j + n_cell_z*n_cell_y*((i+1)%n_cell_x);
                ox_connect[bond_index+5][0] = atom_index+1;
                ox_connect[bond_index+5][1] = 8*next_cell_index + 2;
                //and oxygen 1 in this cell is next to oxygen 6 in the cell one along -y
                next_cell_index = k + n_cell_z*((j-1+n_cell_y)%n_cell_y) + n_cell_z*n_cell_y*i;
                ox_connect[bond_index+6][0] = atom_index+1;
                ox_connect[bond_index+6][1] = 8*next_cell_index + 6;
                //oxygen 2 in this cell is next to oxygen 3 in the same cell
                ox_connect[bond_index+7][0] = atom_index+2;
                ox_connect[bond_index+7][1] = atom_index+3;
                //and oxygen 2 in this cell is also next to oxygen 5 in the same cell
                ox_connect[bond_index+8][0] = atom_index+2;
                ox_connect[bond_index+8][1] = atom_index+5;
                //oxygen 3 in this cell is next to oxygen 4 in the cell above
                next_cell_index = (k+1)%n_cell_z + n_cell_z * j + n_cell_z * n_cell_y * i;
                ox_connect[bond_index+9][0] = atom_index + 3;
                ox_connect[bond_index+9][1] = 8*next_cell_index + 4;
                //oxygen 4 in this cell connects to oxygen 5 in the same cell
                ox_connect[bond_index+10][0] = atom_index + 4;
                ox_connect[bond_index+10][1] = atom_index + 5;
                //oxygen 4 in this cell is also connected to oxygen 7 in the cell below
                next_cell_index = (k-1+n_cell_z) % n_cell_z + n_cell_z * j + n_cell_z * n_cell_y * i;
                ox_connect[bond_index+11][0] = atom_index + 4;
                ox_connect[bond_index+11][1] = 8*next_cell_index + 7;
                //oxygen 4 in this cell is connected to oxygen 7 in the cell below and one along -x
                next_cell_index = (k-1+n_cell_z) % n_cell_z + n_cell_z * j + n_cell_z * n_cell_y * ((i-1+n_cell_x)%n_cell_x);
                ox_connect[bond_index+12][0] = atom_index + 4;
                ox_connect[bond_index+12][1] = 8*next_cell_index + 7;
                //oxygen 5 in this cell is connected to oxygen 6 in the same cell
                ox_connect[bond_index+13][0] = atom_index + 5;
                ox_connect[bond_index+13][1] = atom_index + 6;
                //oxygen 5 in this cell is also connected to oxygen 6 in the cell one along -x
                next_cell_index = k + n_cell_z * j + n_cell_z * n_cell_y * ((i-1+n_cell_x)%n_cell_x);
                ox_connect[bond_index+14][0] = atom_index + 5;
                ox_connect[bond_index+14][1] = 8 * next_cell_index + 6;
                //and finally, oxygen 6 in this cell is connected to oxygen 7 in the same cell
                ox_connect[bond_index+15][0] = atom_index + 6;
                ox_connect[bond_index+15][1] = atom_index + 7;
            }
        }
    }
    
    return ox_connect;
}

//function to initialise starting atomic positions for ice simulation
//generate oxygen atoms on a hexagonal ice lattice
//uses the "hydrogen-ordered" structure of Luo, Yu - 2020
double* initialise_orth(std::mt19937& rng, std::ofstream& outfile) {
    //array to hold positions
    //3*n_atoms dimensional array,
    // with structure O_x, O_y, O_z, H1_x, H1_y, H1_z, H2_x, H2_y, H2_z
    //indexing parameters for loops
    int i, dim, bead;
    int first_oxygen, second_oxygen, old_diff, new_diff, bond_swap;
    int bond_oxygen, pos_index;
    bool is_coord, is_dipole_zero;
    //vector to keep track of what cell we are in
    std::array<double, 3> r_oo, dipole_moment;
    std::vector<int> coord_number(n_mols);
    std::vector<std::array<int, 2> > ox_connect(2*n_mols);
    //store the hydrogen positions in a separate array so that we can rearrange them into a sensible order once the ice rules have been satisfied
    std::vector<double> hyd_pos(6*n_mols);
    //a vector to store which oxygen each hydrogen is currently bonded to
    std::vector<int> hyd_connect(2*n_mols);
    ox_connect = get_oxygen_connectivity();
    
    std::uniform_int_distribution<int> gen(0, 1);
    std::uniform_int_distribution<int> bond_gen(0, 2*n_mols-1);
    
    //fill in the oxygen positions
    double* position_array = initialise_oxygens();
    
    //now we need to make a randomly generated configuration that satisfies the ice rules with zero dipole mometn
    int outer_iter = 1;
    while (true) {
        for (i=0; i < n_mols; ++i) {
            coord_number[i] = 0;
        }
        //loop over all oxygen nearest neighbour pairs
        for (i=0; i < 2*n_mols; ++i) {
            first_oxygen = ox_connect[i][0];
            second_oxygen = ox_connect[i][1];
            //find unit separation vector
            for (dim=0; dim < 3; ++dim) {
                r_oo[dim] = position_array[9*first_oxygen+dim] - position_array[9*second_oxygen+dim];
                //periodic boundary conditions for if we're going outside the box
                if (r_oo[dim] <= -0.5 * boxl[dim]) {
                    r_oo[dim] += boxl[dim];
                }
                else if (r_oo[dim] > 0.5 * boxl[dim]) {
                    r_oo[dim] -= boxl[dim];
                }
                //and normalise the resulting vector
                r_oo[dim] /= rnn;
            }
            
            //pick an oxygen to attach the hydrogen to at random
            if (gen(rng)) {
                for (dim=0; dim < 3; ++dim) {
                    hyd_pos[3*i+dim] = position_array[9*first_oxygen+dim] - roheq * r_oo[dim];
                }
                hyd_connect[i] = first_oxygen;
                coord_number[first_oxygen] += 1;
            }
            else {
                for (dim=0; dim < 3; ++dim) {
                    hyd_pos[3*i+dim] = position_array[9*second_oxygen+dim] + roheq * r_oo[dim];
                }
                hyd_connect[i] = second_oxygen;
                coord_number[second_oxygen] += 1;
            }
        }
        int inner_iter = 1;
        while (true) {
            //check if each oxygen has a coordination number of two
            is_coord = true;
            for (i=0; i < n_mols; ++i) {
                if (coord_number[i] != 2) {
                    is_coord = false;
                }
            }
            if (is_coord) {
                //we are done
                if (outfile.is_open()) {
                    outfile << "configuration found which satisfies the ice rules in " << inner_iter << " iterations \n";
                    outfile.flush();
                }
                else {
                    std::cerr << "couldn't write, outfile appears to be closed \n";
                }
                break;
            }
            inner_iter += 1;
            
            //if not, randomly pick an O-O nearest neighbour pair to swap hydrogens
            bond_swap = bond_gen(rng);
            first_oxygen = ox_connect[bond_swap][0];
            second_oxygen = ox_connect[bond_swap][1];

            //current |c_i - c_j| for this pair
            old_diff = abs(coord_number[first_oxygen] - coord_number[second_oxygen]);
            if (hyd_connect[bond_swap] == first_oxygen) {
                //|c_i - c_j| if an attempted swap is made
                new_diff = abs(coord_number[first_oxygen] - coord_number[second_oxygen] - 2);
                if ((new_diff < old_diff) || ((new_diff==old_diff) && (gen(rng)==0))) {
                    //then the swap is accepted
                    //find the O-O unit separation vector for this pair
                    for (dim=0; dim < 3; ++dim) {
                        r_oo[dim] = position_array[9*first_oxygen+dim] - position_array[9*second_oxygen+dim];
                        //periodic boundary conditions for if we're going outside the box
                        if (r_oo[dim] <= -0.5 * boxl[dim]) {
                            r_oo[dim] += boxl[dim];
                        }
                        else if (r_oo[dim] > 0.5 * boxl[dim]) {
                            r_oo[dim] -= boxl[dim];
                        }
                        //and normalise the resulting vector
                        r_oo[dim] /= rnn;
                    }
                    for (dim=0; dim < 3; ++dim) {
                        hyd_pos[3*bond_swap+dim] = position_array[9*second_oxygen+dim] + roheq * r_oo[dim];
                    }
                    hyd_connect[bond_swap] = second_oxygen;
                    coord_number[second_oxygen] += 1;
                    coord_number[first_oxygen] -= 1;
                }
                //otherwise no change occurs
            }
            //the hydrogen is currently bonded to the second oxygen
            else {
                //|c_i - c_j| if an attempted swap is made
                new_diff = abs(coord_number[first_oxygen] - coord_number[second_oxygen] + 2);
                if ((new_diff < old_diff) || ((new_diff==old_diff) && (gen(rng)==0))) {
                    //then the swap is accepted
                    //find the O-O unit separation vector for this pair
                    for (dim=0; dim < 3; ++dim) {
                        r_oo[dim] = position_array[9*first_oxygen+dim] - position_array[9*second_oxygen+dim];
                        //periodic boundary conditions for if we're going outside the box
                        if (r_oo[dim] <= -0.5 * boxl[dim]) {
                            r_oo[dim] += boxl[dim];
                        }
                        else if (r_oo[dim] > 0.5 * boxl[dim]) {
                            r_oo[dim] -= boxl[dim];
                        }
                        //and normalise the resulting vector
                        r_oo[dim] /= rnn;
                    }
                    for (dim=0; dim < 3; ++dim) {
                        hyd_pos[3*bond_swap+dim] = position_array[9*first_oxygen+dim] - roheq * r_oo[dim];
                    }
                    hyd_connect[bond_swap] = first_oxygen;
                    coord_number[first_oxygen] += 1;
                    coord_number[second_oxygen] -= 1;
                }
                //otherwise no change occurs
            }
        }
        
        //keep track of how many bonds each oxygen currently has
        std::vector<double> ox_bonds(n_mols);
        for (i=0; i < n_mols; ++i) {
            ox_bonds[i] = 0;
        }
        
        //loop over all hydrogen atoms
        for (i=0; i < 2*n_mols; ++i) {
            //find the index of the oxygen bonded to it
            bond_oxygen = hyd_connect[i];
            //increment the number of bonds this oxygen currently has by one
            ox_bonds[bond_oxygen] += 1;
            pos_index = 9 * bond_oxygen + 3 * ox_bonds[bond_oxygen];
            for (dim=0; dim < 3; ++dim) {
                //store the position at the correct point in the array
                //so that its structure is such that after each oxygen come the two hydrogens which are bonded to it
                position_array[pos_index+dim] = hyd_pos[3*i+dim];
            }
        }
        
        dipole_moment = get_dipole_moment(position_array);
        if (outfile.is_open()) {
            outfile << "dipole moment is " << dipole_moment[0] << " " << dipole_moment[1] << " " << dipole_moment[2] << "\n";
            outfile.flush();
        }
        else {
            std::cerr << "couldn't write, outfile appears to be closed \n";
        }
        is_dipole_zero = true;
        for (dim=0; dim < 3; ++dim) {
            if (abs(dipole_moment[dim]) >= 1e-10) {
                is_dipole_zero = false;
            }
        }
        if (is_dipole_zero) {
            if (outfile.is_open()) {
                outfile << "configuration found that has zero dipole moment in " << outer_iter << " outer iterations \n";
                outfile.flush();
            }
            else {
                std::cerr << "couldn't write, outfile appears to be closed \n";
            }
            break;
        }
        outer_iter += 1;
    }
    
    //set initial positions in all other replicas to those in the first replica
    for (bead=1; bead < n_beads; ++bead) {
        for (i=0; i < n_dof; ++i) {
            position_array[bindex(bead, i)] = position_array[i];
        }
    }
    
    //return the initial positions
    return position_array;
}

double* initialise_positions(std::string cell_type, std::mt19937& rng, std::ofstream& outfile) {
    
    double* position_array;
    
    if (cell_type == "fcc") {
        position_array = initialise_fcc();
    }
    
    else {
        position_array = initialise_orth(rng, outfile);
    }
    
    return position_array;
}

