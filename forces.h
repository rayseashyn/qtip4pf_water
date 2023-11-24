#ifndef FORCES
#define FORCES

double lj_oxygen(const std::array<double, 3>& r_ij, double rijsq, double force_array[], int i, int j, int bead=0);

double intra_morse_bond(const std::array<double, 3>& r_oh, double roh, double force_array[], int i, int bond_step, int bead=0);

double intra_harm_angle(const std::array<double, 3>& r_oh1, const std::array<double, 3>& r_oh2, double roh1, double roh2, double force_array[], int i, int bead=0);

double rwald(const std::array<double, 3>& r_ij, double rijsq, double force_array[], const double charge_array[], int i, int j, double ewald_kappa, int bead=0);

double kwald(const double position_array[], double force_array[], const double charge_array[], int k_max, double ewald_kappa, double k_cut_sq, int bead=0);

double get_spline(const std::array<double, 3> r_ij, const double rij, double force_array[], const int i, const int j, const boost::math::cubic_b_spline<double>& spline, int bead=0);

double get_forces_nospline(const double position_array[], double force_array[], const double charge_array[], const int point[], const int list[], int bead=0);

double get_forces_wspline(const double position_array[], double force_array[], const double charge_array[], const int point[], const int list[], const boost::math::cubic_b_spline<double>& oo_spline, const boost::math::cubic_b_spline<double>& oh_spline, const boost::math::cubic_b_spline<double>& hh_spline, const double oo_spline_cutsq, const double oh_spline_cutsq, const double hh_spline_cutsq, int bead=0);

double force_eval_sc_nospline(const double position_array[], double force_array[], const double charge_array[], const double mass_array[], const int point[], const int list[], double& sc_pot, double fphys[]=nullptr);

double force_eval_sc_wspline(const double position_array[], double force_array[], const double charge_array[], const double mass_array[], const int point[], const int list[], const boost::math::cubic_b_spline<double>& oo_spline, const boost::math::cubic_b_spline<double>& oh_spline, const boost::math::cubic_b_spline<double>& hh_spline, const double oo_spline_cutsq, const double oh_spline_cutsq, const double hh_spline_cutsq, double& sc_pot, double fphys[]=nullptr);

#endif
