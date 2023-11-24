//a function to convert arrays from the "bead representation" to the "normal mode representation"
//for calculations using path integral molecular dynamics
//utilises the fast fourier transform from the fftw3 library for efficiency

#include <fftw3.h>
#include <vector>
#include "fourier.h"
#include "qtip4pfparams.h"

//compute fourier transform of input array using fftw
//transforms to the normal mode representation if is_forwards=1 and to the bead representation if is_forwards=0
void ftransform (double input_array[], const int length, const bool is_forwards) {
    //input array is either q or p, being n_beads x n_dof arrays
    //plan is for the forwards or backwards fourier transform, depending on which is needed
    //1D array to do FT on
    static std::vector<double> qj(n_beads);
    static fftw_plan p_forwards = fftw_plan_r2r_1d(n_beads, qj.data(), qj.data(), FFTW_R2HC, FFTW_MEASURE);
    static fftw_plan p_backwards = fftw_plan_r2r_1d(n_beads, qj.data(), qj.data(), FFTW_HC2R, FFTW_MEASURE);
    
    //indexing parameters
    int bead, i;
    double scale = sqrt(1.0 / n_beads);
    double twoscale = sqrt(2.0 / n_beads);
    double root2 = sqrt(2.0);
    //account for the fact that the arrays to be fourier transformed could be different sizes
    //the first index is always n_beads, and the second is the dimension of the array
    //eg for force_array, position_array and momentum_array this will be 9*n_mols
    //but for disp_array this will be 3*n_mols as only oxygen atoms are included in the neighbour list
    //int length = input_array[0].size();
    
    for (i=0; i < length; ++i) {
        for (bead=0; bead < n_beads; ++bead) {
            qj[bead] = input_array[bead*length+i];
        }
        if (is_forwards) {
            //convert array to normal mode representation
            fftw_execute(p_forwards);
            
            input_array[i] = scale * qj[0];
            //all the degenerate normal modes need a factor of sqrt(2) along with the Fourier transform
            //this is required for thermostatting - to ensure that the kinetic energy is the same in both
            //the bead and normal mode representations, for instance
            for (bead=1; bead < n_beads; ++bead) {
                input_array[bead*length+i] = twoscale * qj[bead];
            }
            //if there are an even number of beads, then there is an extra non-degenerate normal mode
            //(the highest frequency one) - this appears at index n_beads/2 in the Fourier transformed array
            //so we get rid of the factor of sqrt(2) we added earlier in this case
            if (n_beads%2==0) {
                input_array[n_beads*length/2+i] /= root2;
            }
        }
        else {
            //remove the sqrt(2) scaling of the degenerate normal modes done in the forwards transform, to ensure that
            //we transform back to the bead representation correctly
            for (bead=1; bead < n_beads; ++bead) {
                qj[bead] /= root2;
            }
            if (n_beads%2==0) {
                qj[n_beads/2] *= root2;
            }
            //now convert array to bead representation
            fftw_execute(p_backwards);
            for (bead=0; bead < n_beads; ++bead) {
                input_array[bead*length+i] = scale * qj[bead];
            }
        }
    }
}
