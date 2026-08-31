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
    double qc; ///< Central safety factor q_c
    double qa; ///< Edge safety factor (q(r) at r=a)
    double pc; ///< Central kinetic pressure (Pa)
    double epsa; ///<  Inverse aspect ratio a/R0 
    double B0; ///< Toroidal magnetic field (T)
    int Nr; ///< radial grid resolution
    double mu; ///< Current profile exponent

    std::vector<double> r; ///< Normalized radial coordinate r = r/a in [0, 1]
    std::vector<double> q_pre; ///< Pre-crash safety factor profile q_pre(r)
    std::vector<double> p; ///< Pre-crash plasma pressure profile p(r)
    std::vector<double> B_theta; ///< Poloidal magnetic field profile B_theta(r) [Tesla]
    std::vector<double> dpsi_star_dr; ///< Radial derivative of helical flux d(psi_star)/dr
    std::vector<double> psi_star; ///< Helical magnetic flux profile psi_star(r)
    std::vector<double> q_post; ///< Post-crash relaxed safety factor profile q_post(r)


    double r1; ///< Resonant q=1 surface radius (r1 / a)
    double r0; ///< Boundary radius (r0 / a)



    void validate_inputs() const;
    void allocate_memory();
    void compute_profiles();
    void integrate_helical_flux();
    void run_post_solve_sanity_checks() const;

    double interpolate(const std::vector<double>& x_arr, const std::vector<double>& y_arr, double x_val) const;
    double find_root_bisection(double (*func)(double, const Kadomtsev*), double low, double high, double tol = 1e-6);

public:
    Kadomtsev(double qa_val, double pc_val, double epsa_val, double B0_val, double qc_val, int Nr_val);

    void solve();
    void print_summary() const;
    void write_data(const std::string& filename) const;

    double eval_q_minus_one(double radial_r) const;
    double eval_psi_diff(double radial_r) const;
};

#endif // KADOMTSEV_H