#ifndef UPDATE
#define UPDATE

double update_trajectory_class_wspline(double position_array[], double momentum_array[], const double mass_array[], double force_array[], double charge_array[], int point[], int list[], double disp_array[], const boost::math::cubic_b_spline<double>& oo_spline, const boost::math::cubic_b_spline<double>& oh_spline, const boost::math::cubic_b_spline<double>& hh_spline, const double oo_spline_cutsq, const double oh_spline_cutsq, const double hh_spline_cutsq);
double update_trajectory_class_nospline(double position_array[], double momentum_array[], const double mass_array[], double force_array[], double charge_array[], int point[], int list[], double disp_array[]);
double update_sc_wspline(double position_array[], double momentum_array[], const double mass_array[], double force_array[], double fphys[], double charge_array[], int point[], int list[], double disp_array[], std::mt19937& rng, bool is_nvt, double& heat, double& sc_pot, const boost::math::cubic_b_spline<double>& oo_spline, const boost::math::cubic_b_spline<double>& oh_spline, const boost::math::cubic_b_spline<double>& hh_spline, const double oo_spline_cutsq, const double oh_spline_cutsq, const double hh_spline_cutsq);
double update_sc_nospline(double position_array[], double momentum_array[], const double mass_array[], double force_array[], double fphys[], double charge_array[], int point[], int list[], double disp_array[], std::mt19937& rng, bool is_nvt, double& heat, double& sc_pot);

#endif
