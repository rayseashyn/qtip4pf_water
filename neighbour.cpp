//function to construct a neighbour list so that short-range intermolecular interactions are only
//considered between atom pairs sufficiently close to each other
//and the list updated after atoms have moved far enough since the last update

#include <vector>
#include <iostream>
#include "neighbour.h"
#include "qtip4pfparams.h"
#include "indexing.h"

//generate neighbour list for lennard-jones interactions between oxygen atoms for each of the replicas
void neighbour_list(const double position_array[], int point[], int list[], double disp_array[], int bead) {
    //indexing parameters
    int i, j, k, dim;
    
    //squared distance
    double r_ij_sq;
    
    //array to hold displacement vectors between atoms i and j in Cartesian and cell units
    std::array<double, 3> r_ij;
    
    //set the neighbour list for this replica to zeros
    for (i=0; i < list_size; ++i) {
        list[lindex(bead, i)] = 0;
    }
    
    //and reset the displacement vectors for each replica back to zero, as its neighbour list is now being updated
    for (i=0; i < n_mols; ++i) {
        for (dim=0; dim < 3; ++dim) {
            disp_array[dindex(bead, 3*i+dim)] = 0.0;
        }
    }
    
    //construct the verlet neighbour list for this replica
    k = 0;
    for (i=0; i < n_mols; ++i) {
        point[pindex(bead, i)] = k;
        //loop over all oxygen atoms with index j>i
        for (j=i+1; j < n_mols; ++j) {
            for (dim=0; dim < 3; ++dim) {
                //only calculate separation vector between oxygen atoms
                r_ij[dim] = position_array[bindex(bead, 9*i+dim)] - position_array[bindex(bead, 9*j+dim)];
            }
            
            for (dim=0; dim < 3; ++dim) {
                r_ij[dim] -= boxl[dim] * std::round(r_ij[dim] / boxl[dim]);
            }
            
            //compute squared distance
            r_ij_sq = r_ij[0]*r_ij[0] + r_ij[1]*r_ij[1] + r_ij[2]*r_ij[2];
            
            //check if we should store these pairs in the neighbour list for this bead
            if (r_ij_sq < r_skin_sq) {
                list[lindex(bead, k)] = j;
                k += 1;
                if (k >= list_size) {
                    std::cerr << "list size exceeded!\n";
                }
            }
        }
    }
}
