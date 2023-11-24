//defining functions for quick, intuitive way to find the correct place
//in an array for the different arrays used in the simulation
#ifndef INDEXING
#define INDEXING

inline int bindex(int bead, int i) {
    return bead*n_dof + i;
};

inline int kindex(int k, int i) {
    return n_dof*k + i;
};

inline int lindex(int bead, int i) {
    return list_size*bead + i;
};

inline int dindex(int bead, int i) {
    return 3*n_mols*bead + i;
};

inline int pindex(int bead, int i) {
    return n_mols*bead + i;
};

#endif
