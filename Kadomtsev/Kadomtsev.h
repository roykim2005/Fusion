#ifndef KADOMTSEV_H
#define KADOMTSEV_H

#include <vector>
#include <string>

/**
 * @class Kadomtsev
 * @brief Computes pre/post-sawtooth crash safety factor q(r) profiles using Kadomtsev helical flux reconnection model.
 */
class Kadomtsev {
private:
    double qc; // Central safety factor q_c
    double qa; // Edge safety factor (q(r) at r=a)
    double pc; // Central kinetic pressure (Pa)
    double epsa; //  Inverse aspect ratio a/R0 
    double B0; // Toroidal magnetic field (T)
    int Nr; // radial grid resolution
    double mu; // Current profile exponent
    double chi; // Thermal diffusivity [m^2/s]
    double n_e; // Electron density [m^-3]

    std::vector<double> r; // Normalized radial coordinate r = r/a in [0, 1]
    std::vector<double> q_pre; // Pre-crash safety factor profile q_pre(r)
    std::vector<double> p; // Pre-crash plasma pressure profile p(r)
    std::vector<double> B_theta; // Poloidal magnetic field profile B_theta(r) [Tesla]
    std::vector<double> dpsi_star_dr; // Radial derivative of helical flux d(psi_star)/dr
    std::vector<double> psi_star; // Helical magnetic flux profile psi_star(r)
    std::vector<double> q_post; // Post-crash relaxed safety factor profile q_post(r)
    std::vector<double> T_pre; // Pre-crash temperature profile
    std::vector<double> T_post; // Post-crash mixed temperature profile
    std::vector<double> eta_pre; // Pre-crash Spitzer resistivity
    std::vector<double> eta_post; // Post-crash Spitzer resistivity

    double r1; // Resonant q=1 surface radius (r1 / a)
    double r0; // Boundary radius (r0 / a)
    double Tc; // Central temperature (measured in keV)
    double current_time{0.0}; // for time evolution


    void validate_inputs() const;
    void allocate_memory();
    void compute_profiles();
    void integrate_helical_flux();
    void run_post_solve_sanity_checks() const;

    double interpolate(const std::vector<double>& x_arr, const std::vector<double>& y_arr, double x_val) const;
    double find_root_bisection(double (*func)(double, const Kadomtsev*), double low, double high, double tol = 1e-6);
    double find_r2_outer(double r1_inner) const;
    void compute_temperature_and_resistivity();

public:
    Kadomtsev(double qa_val, double pc_val, double epsa_val, double B0_val, double qc_val, int Nr_val);

    void solve();
    void print_summary() const;
    void set_current_time(double t) { current_time = t; }
    void write_data(const std::string& filename) const;

    double eval_q_minus_one(double radial_r) const;
    double eval_psi_diff(double radial_r) const;

    // time dependence
    void reset_post_to_pre();
    void step_ramp_phase(double dt);

    double get_q0() const { return q_pre[0]; }
    double get_T0() const { return T_pre[0]; }
};

#endif // KADOMTSEV_H
