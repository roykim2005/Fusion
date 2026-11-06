#include "Kadomtsev.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <stdexcept>
#include <algorithm>
#include <utility>
#include <numbers>


// finds r1 such that q(r1)=1
static double q_root_wrapper(double radial_r, const Kadomtsev* model) {
    return model->eval_q_minus_one(radial_r);
}
// finds r0 such that psi(r0)=psi(0)
static double psi_root_wrapper(double radial_r, const Kadomtsev* model) {
    return model->eval_psi_diff(radial_r);
}
// imports parameters from Kadomtsev.json
Kadomtsev::Kadomtsev(double qa_val, double qc_val, int Nr_val, double kappa_val, double theta_val, double lam_val, double Tedge_hat_val, double alpha_chi_val)
    : qc(qc_val), qa(qa_val), Nr(Nr_val), kappa(kappa_val), theta(theta_val), lam(lam_val), Tedge_hat(Tedge_hat_val), alpha_chi(alpha_chi_val) {
    // sanity check
    validate_inputs();
    mu = qa / qc;
    // resizes arrays std::vector for Nr points for r, pre q, post q, B0, p, psi.
    allocate_memory();
    compute_profiles();
}

void Kadomtsev::validate_inputs() const {
    if (Nr < 100) {
        throw std::invalid_argument("Nr must be at least 100.");
    }
    if (qa <= qc) {
        throw std::invalid_argument("Edge safety factor qa must be greater than central parameter qc.");
    }
    if (qc >= 1.0) {
        throw std::domain_error("qc is already at equilibrium.");
    }
}

void Kadomtsev::allocate_memory() {
    r.resize(Nr);
    q_pre.resize(Nr);
    p.resize(Nr);
    B_theta.resize(Nr);
    dpsi_star_dr.resize(Nr);
    psi_star.resize(Nr);
    q_post.resize(Nr);
    T_pre.resize(Nr);
    T_post.resize(Nr);
    eta_pre.resize(Nr);
    eta_post.resize(Nr);
}

// 1D interpolation function
double interpolate_1d(const std::vector<double>& x_vals, 
                       const std::vector<double>& y_vals, 
                       double x_target) 
{
    if (x_vals.empty()) return 0.0;
    if (x_target <= x_vals.front()) return y_vals.front();
    if (x_target >= x_vals.back()) return y_vals.back();

    auto it = std::lower_bound(x_vals.begin(), x_vals.end(), x_target);
    size_t idx = std::distance(x_vals.begin(), it);

    if (idx == 0) return y_vals[0];

    double x0 = x_vals[idx - 1], x1 = x_vals[idx];
    double y0 = y_vals[idx - 1], y1 = y_vals[idx];

    return y0 + (x_target - x0) * (y1 - y0) / (x1 - x0);
}

// interpolates 
double Kadomtsev::interpolate(const std::vector<double>& x_arr, const std::vector<double>& y_arr, double x_val) const {
    if (x_val <= x_arr.front()) return y_arr.front();
    if (x_val >= x_arr.back()) return y_arr.back();

    for (size_t i = 0; i < x_arr.size() - 1; ++i) {
        if (x_val >= x_arr[i] && x_val <= x_arr[i + 1]) {
            double t = (x_val - x_arr[i]) / (x_arr[i + 1] - x_arr[i]);
            return y_arr[i] + t * (y_arr[i + 1] - y_arr[i]);
        }
    }
    return y_arr.back();
}
// finds r1 for q-1 = 0 and r0 for psi(r0)-psi(0) = 0
double Kadomtsev::find_root_bisection(double (*func)(double, const Kadomtsev*), double low, double high, double tol) {
    double f_low = func(low, this);
    double f_high = func(high, this);

    if (f_low * f_high > 0.0) {
        throw std::runtime_error("Root not in interval [" + 
                                 std::to_string(low) + ", " + std::to_string(high) + "].");
    }

    double mid = low;
    while ((high - low) / 2.0 > tol) {
        mid = low + (high - low) / 2.0;
        double f_mid = func(mid, this);

        if (std::abs(f_mid) < 1e-12) break;

        if (f_mid * f_low < 0.0) {
            high = mid;
            f_high = f_mid;
        } else {
            low = mid;
            f_low = f_mid;
        }
    }
    return mid;
}

double Kadomtsev::eval_q_minus_one(double radial_r) const {
    return interpolate(r, q_pre, radial_r) - 1.0;
}

double Kadomtsev::eval_psi_diff(double radial_r) const {
    return interpolate(r, psi_star, radial_r) - psi_star[0];
}

// temperature model
double Kadomtsev::find_r2_outer(double r1_inner) const {
    double target_psi = interpolate(r, psi_star, r1_inner);
    
    // Bisection search between r1 (q=1 surface) and r0 (mixing boundary)
    double low = r1;
    double high = r0;
    double mid = low;
    
    while ((high - low) / 2.0 > 1e-6) {
        mid = low + (high - low) / 2.0;
        double f_mid = interpolate(r, psi_star, mid) - target_psi;
        
        if (std::abs(f_mid) < 1e-10) break;

        double f_low = interpolate(r, psi_star, low) - target_psi;
        if (f_mid * f_low < 0.0) {
            high = mid;
        } else {
            low = mid;
        }
    }
    return mid;
}

void Kadomtsev::compute_temperature_and_resistivity() {
    std::vector<double> r_post_list;
    std::vector<double> T_mix_list;

    // reconnection surface counts
    int N_reconnect = 300; 

    // region r1 to 0
    for (int k = 0; k < N_reconnect; ++k) {
        double frac = (static_cast<double>(k) + 0.5) / N_reconnect;
        double r1_minus = r1 * frac;
        
        // finds r2_plus such that psi*(r1_minus) == psi*(r2_plus)
        double r2_plus = find_r2_outer(r1_minus);

        double r_post = std::sqrt(std::max(0.0, r2_plus * r2_plus - r1_minus * r1_minus));

        // Calculate derivative ratio dpsi/dr
        double dpsi_dr1 = std::abs(interpolate_1d(r, dpsi_star_dr, r1_minus));
        double dpsi_dr2 = std::abs(interpolate_1d(r, dpsi_star_dr, r2_plus));
        double dr1_dr2 = (dpsi_dr1 > 1e-3) ? (dpsi_dr2 / dpsi_dr1) : 1.0;

        double T1 = interpolate_1d(r, T_pre, r1_minus);
        double T2 = interpolate_1d(r, T_pre, r2_plus);

        // temp relaxation
        double denom = r1_minus * dr1_dr2 + r2_plus;
        double T_mix = (denom > 1e-9) ? ((T1 * r1_minus * dr1_dr2 + T2 * r2_plus) / denom) : T1; // use average, find how q is modified after crash

        r_post_list.push_back(r_post);
        T_mix_list.push_back(T_mix);
    }

    // for plot python plot generation, may or may not need this. currently it looks like its mapping things backwards, this should sort from low to high r_post
    std::vector<std::pair<double, double>> pairs(r_post_list.size());
    for (size_t i = 0; i < pairs.size(); ++i) {
        pairs[i] = {r_post_list[i], T_mix_list[i]};
    }
    std::sort(pairs.begin(), pairs.end());

    for (size_t i = 0; i < pairs.size(); ++i) {
        r_post_list[i] = pairs[i].first;
        T_mix_list[i]  = pairs[i].second;
    }

    // this calculates the core temp and volume averages the temp.
    double thermal_energy_integral = 0.0;
    double volume_integral = 0.0;
    for (int i = 0; i < Nr && r[i] <= r0; ++i) {
        double dr = (i > 0) ? (r[i] - r[i-1]) : r[0];
        thermal_energy_integral += T_pre[i] * r[i] * dr;
        volume_integral += r[i] * dr;
    }
    double T_core_avg = thermal_energy_integral / volume_integral;


    // spatial grid
    int idx_r0 = 0;
    while (idx_r0 < Nr - 1 && r[idx_r0] < r0) {
        idx_r0++;
    }


    for (int i = 0; i < Nr; ++i) {
        if (r[i] <= r0) { // it was doing T_post = T_core_avg before which just averages all T_post, it didn't call r[i] <= r1 because that never happened
            T_post[i] = interpolate_1d(r_post_list, T_mix_list, r[i]);
        } else {
            T_post[i] = T_pre[i]; // outer
        }
        T_post[idx_r0] = T_pre[idx_r0];

        // new resistivity
        eta_post[i] = std::pow(std::max(Tedge_hat, T_post[i]), -1.5);
    }
    
}

// sets up a 1D grid and calculates pre crash q, pre crash p, B0, d psi/dr = (1-q(r))B0, if T = 2 keV
void Kadomtsev::compute_profiles() {
    double dr = (1.0 - 1e-5) / (Nr - 1);

    mu = qa / qc;
    for (int i = 0; i < Nr; ++i) {
        r[i] = 1e-5 + i * dr;

        q_pre[i] = qc + (qa - qc) * r[i] * r[i]; // q profile
        double base_profile = std::max(0.01, 1.0 - r[i] * r[i]);
        
        T_pre[i] = std::pow(base_profile, mu); // temp, Tc is now 1
        eta_pre[i] = std::pow(std::max(Tedge_hat, T_pre[i]), -1.5); // resistivity
        B_theta[i] = r[i] / q_pre[i]; // B field
        dpsi_star_dr[i] = (1.0 - q_pre[i]) * B_theta[i];
    }

    f_chi = (1.0 + alpha_chi) / (std::pow(2.0, 1.0 + alpha_chi) - 1.0);
    chi_half.resize(Nr - 1);
    for (int i = 0; i < Nr - 1; ++i) {
        double r_h = r[i] + 0.5 * dr;
        chi_half[i] = f_chi * std::pow(1.0 + r_h * r_h, alpha_chi);
    }
}


// integrates d psi/dr
void Kadomtsev::integrate_helical_flux() {
    psi_star[0] = 0.0;
    for (int i = 1; i < Nr; ++i) {
        double dr = r[i] - r[i - 1];
        psi_star[i] = psi_star[i - 1] + 0.5 * (dpsi_star_dr[i - 1] + dpsi_star_dr[i]) * dr;
    }
}

void Kadomtsev::solve() {
    for (int i = 0; i < Nr; ++i)
        dpsi_star_dr[i] = (1.0 - q_pre[i]) * B_theta[i]; // now it recomputes the d psi*/dr before 
    integrate_helical_flux();

    r1 = find_root_bisection(q_root_wrapper, r[0], r[Nr - 1]);
    r0 = find_root_bisection(psi_root_wrapper, r1 + 1e-3, 0.99);

    double q_r0 = interpolate(r, q_pre, r0);
    for (int i = 0; i < Nr; ++i) {
        // internal and external, q=1 inside and nothing happens outside
        if (r[i] <= r0) {
            q_post[i] = 1.0 + (q_r0 - 1.0) * std::pow(r[i] / r0, 2);
        } else {
            q_post[i] = q_pre[i];
        }
    }
    compute_temperature_and_resistivity();
    run_post_solve_sanity_checks();

    std::cout << "Crash at t_hat = " << std::fixed << std::setprecision(6) << current_time << std::endl;
}

void Kadomtsev::run_post_solve_sanity_checks() const {
    // resonance surface and boundary radius surface sanity check
    if (!(r1 > 0.0 && r1 < r0 && r0 < 1.0)) {
        throw std::runtime_error("r0 and r1 must satisfy 0 < r1 < r0 < 1.");
    }
    // pre crash sanity check
    if (q_pre[0] >= 1.0) {
        throw std::runtime_error("Pre-crash q(0) must be less than 1.0.");
    }

    // post crash sanity check
    if (std::abs(q_post[0] - 1.0) > 1e-3) {
        throw std::runtime_error("Post-crash q(0) did not reset to 1.0.");
    }

    // continuity check
    double q_pre_r0 = interpolate(r, q_pre, r0);
    double q_post_r0 = interpolate(r, q_post, r0);
    if (std::abs(q_pre_r0 - q_post_r0) > 1e-3) {
        throw std::runtime_error("Discontinuity in q_post at r0.");
    }
}

// post crash, resets the conditions. maybe there's a better model that accounts for previous conditions 
void Kadomtsev::reset_post_to_pre() {
    for (int i = 0; i < Nr; ++i) {
        T_pre[i] = T_post[i];
        q_pre[i] = q_post[i];
        eta_pre[i] = eta_post[i];
        B_theta[i] = r[i] / q_pre[i];
    }
    be_steps_left = 4;
}


// tridiagonal matrix solver for big matrices from time dependence codes
static std::vector<double> solve_tridiagonal(const std::vector<double>& A, const std::vector<double>& B, const std::vector<double>& C, std::vector<double> D) {
    int N = D.size();
    std::vector<double> c_prime(N, 0.0);
    std::vector<double> x(N, 0.0);

    c_prime[0] = C[0] / B[0];
    D[0] = D[0] / B[0];

    for (int i = 1; i < N; ++i) {
        double m = 1.0 / (B[i] - A[i] * c_prime[i - 1]);
        c_prime[i] = C[i] * m;
        D[i] = (D[i] - A[i] * D[i - 1]) * m;
    }

    x[N - 1] = D[N - 1];
    for (int i = N - 2; i >= 0; --i) {
        x[i] = D[i] - c_prime[i] * x[i + 1];
    }

    return x;
}


// transport equation solver


void Kadomtsev::step_ramp_phase(double dt) { // dt is dt hat
    const int n_iter = 2;
    double dr = r[1] - r[0];
    double w = (be_steps_left > 0) ? 1.0 : 0.5; // 0.5 is CN, 1 is backward euler


    double alpha = kappa * dt / (dr * dr);

    auto make_eta = [&](const std::vector<double>& T, std::vector<double>& eta) {
        for (int i = 0; i < Nr; ++i) {
            double Tf = std::max(Tedge_hat, T[i]);
            eta[i] = 1.0 / (Tf * std::sqrt(Tf));
        }
    };
        auto make_j = [&](const std::vector<double>& u, std::vector<double>& j) {
        for (int i = 1; i < Nr - 1; ++i) j[i] = (u[i + 1] - u[i - 1]) / (2.0 * dr * r[i]);
        j[0] = j[1];
        j[Nr - 1] = (u[Nr - 1] - u[Nr - 2]) / (dr * r[Nr - 1]);
    };

    // initial timestep
    std::vector<double> T_old = T_pre, u_old(Nr), eta_old(Nr), j_old(Nr, 0.0), S_old(Nr);
    for (int i = 0; i < Nr; ++i) u_old[i] = r[i] * B_theta[i];
    make_eta(T_old, eta_old);
    make_j(u_old, j_old);
    for (int i = 0; i < Nr; ++i) S_old[i] = eta_old[i] * j_old[i] * j_old[i];

    // next step, one before is the first guess
    std::vector<double> T_new = T_old, u_new = u_old, eta_new = eta_old, j_new = j_old, S_new = S_old;
    std::vector<double> eta_w(Nr), S_w(Nr);

    for (int it = 0; it < n_iter; ++it) {
        for (int i = 0; i < Nr; ++i) {
            eta_w[i] = (1.0 - w) * eta_old[i] + w * eta_new[i];
            S_w[i]   = (1.0 - w) * S_old[i]   + w * S_new[i];
        }
        
        // temp
        std::vector<double> A(Nr, 0.0), B(Nr, 1.0), C(Nr, 0.0), D(T_old);
        B[0] = 1.0; C[0] = -1.0; D[0] = 0.0;
        for (int i = 1; i < Nr - 1; ++i) {
            double r_m = r[i] - 0.5 * dr, r_p = r[i] + 0.5 * dr;
            double chi_m = chi_half[i - 1], chi_p = chi_half[i];
            A[i] = -w * alpha * chi_m * r_m / r[i];
            C[i] = -w * alpha * chi_p * r_p / r[i];
            B[i] = 1.0 - A[i] - C[i];
            double lap = (chi_p * r_p / r[i]) * (T_old[i + 1] - T_old[i])
                       - (chi_m * r_m / r[i]) * (T_old[i] - T_old[i - 1]);
            D[i] = T_old[i] + (1.0 - w) * alpha * lap + dt * theta * S_w[i];
        }
        D[Nr - 1] = Tedge_hat;
        T_new = solve_tridiagonal(A, B, C, D);
        for (int i = 0; i < Nr; ++i) T_new[i] = std::max(Tedge_hat, T_new[i]);

        // let u = rB_theta, jhat = 1/r du/dr
        std::vector<double> Au(Nr, 0.0), Bu(Nr, 1.0), Cu(Nr, 0.0), Du(u_old);
        for (int i = 1; i < Nr - 1; ++i) {
            double r_m = r[i] - 0.5 * dr, r_p = r[i] + 0.5 * dr;
            double eta_m = 0.5 * (eta_w[i] + eta_w[i - 1]);
            double eta_p = 0.5 * (eta_w[i] + eta_w[i + 1]);
            double beta = lam * dt * r[i] / (dr * dr);
            Au[i] = -w * beta * eta_m / r_m;
            Cu[i] = -w * beta * eta_p / r_p;
            Bu[i] = 1.0 - Au[i] - Cu[i];
            Du[i] = u_old[i] + (1.0 - w) * beta * (eta_p * (u_old[i + 1] - u_old[i]) / r_p - eta_m * (u_old[i] - u_old[i - 1]) / r_m);
        }
        Du[0] = 0.0; // u(0)=0.
        u_new = solve_tridiagonal(Au, Bu, Cu, Du);

        // new coeff
        make_eta(T_new, eta_new);
        make_j(u_new, j_new);
        for (int i = 0; i < Nr; ++i) S_new[i] = eta_new[i] * j_new[i] * j_new[i];
    }

    T_pre = T_new;
    for (int i = 1; i < Nr; ++i) {
        B_theta[i] = u_new[i] / r[i];
        q_pre[i] = r[i] / B_theta[i];
    }
    q_pre[0] = q_pre[1] - (q_pre[2] - q_pre[1]) * (r[1]*r[1] - r[0]*r[0]) / (r[2]*r[2] - r[1]*r[1]);
    B_theta[0] = r[0] / q_pre[0];
    for (int i = 0; i < Nr; ++i) {
        double Tf = T_pre[i];
        eta_pre[i] = 1.0 / (Tf * std::sqrt(Tf));
    }
    if (be_steps_left > 0) --be_steps_left;
}

void Kadomtsev::print_summary() const {
    double q_r0_pre = interpolate(r, q_pre, r0);
    double q_r0_post = interpolate(r, q_post, r0);

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "============================================================== \n"; 

    std::cout << "Resonant Surface (r1): " << r1 << "\n";
    std::cout << "Boundary Radius (r0): " << r0 << "\n";
    std::cout << "q(0) - pre-crash: " << q_pre[0] << ", post-crash: " << q_post[0] << "\n";
    std::cout << "q(r1/4) - pre-crash: " << interpolate(r, q_pre, r1/4.0) << ", post-crash: " << interpolate(r, q_post, r1/4.0) << "\n";
    std::cout << "q(r1/2) - pre-crash: " << interpolate(r, q_pre, r1/2.0) << ", post-crash: " << interpolate(r, q_post, r1/2.0) << "\n";
    std::cout << "q(3r1/4) - pre-crash: " << interpolate(r, q_pre, 3*r1/4.0) << ", post-crash: " << interpolate(r, q_post, 3*r1/4.0) << "\n";
    std::cout << "q(r1) - pre-crash: " << interpolate(r, q_pre, r1) << ", post-crash: " << interpolate(r, q_post, r1) << "\n";
    std::cout << "q(r0) - pre-crash: " << q_r0_pre << ", post-crash: " << q_r0_post << "\n";
    std::cout << "============================================================== \n"; 

    std::cout << "T_hat(0) - pre-crash: " << T_pre[0] << ", post-crash: " << T_post[0] << " \n";
    std::cout << "T_hat(r1/2) - pre-crash: " << interpolate(r, T_pre, r1/2.0) << ", post-crash: " << interpolate(r, T_post, r1/2.0) << "\n";
    // T(r1) post is the temp after merging
    std::cout << "T_hat(r1) - pre-crash: " << interpolate(r, T_pre, r1) << ", post-crash: " << interpolate(r, T_post, r1) << "\n";
    std::cout << "T_hat(r0) - pre-crash: " << interpolate(r, T_pre, r0) << ", post-crash: " << interpolate(r, T_post, r0) << "\n";
    std::cout << std::scientific << std::setprecision(3);
    std::cout << "eta_hat(0) - pre-crash: " << eta_pre[0] << ", post-crash: " << eta_post[0] << "\n";
    std::cout << "eta_hat(r1) - pre-crash: " << interpolate(r, eta_pre, r1) << ", post-crash: " << interpolate(r, eta_post, r1) << "\n";
}



void Kadomtsev::write_data(const std::string& filename) const {
    std::ofstream file(filename);

    file << "# r/a\tq_pre\tq_post\tT_pre/Tc\tT_post/Tc\teta_pre/eta_c\teta_post/eta_c\tpsi_star\n";

    file << std::scientific << std::setprecision(6);

    for (int i = 0; i < Nr; ++i) {
        file << r[i] << "\t" 
             << q_pre[i] << "\t" 
             << q_post[i] << "\t" 
             << T_pre[i] << "\t" 
             << T_post[i] << "\t" 
             << eta_pre[i] << "\t" 
             << eta_post[i] << "\t" 
             << psi_star[i] << "\n";
    }

    file.close();
    std::cout << "Exported to: " << filename << std::endl;
}


// suggestion: first timestep, i should see something reasonable and nothing crazy. Boundary condition for triangular matrix. make sure i have bc right. 
// near axis, do finite difference. USe crank nicholson and you end up with tridiagonal. 

// plot T vs r and q vs r, differnt stages of crash

// timestep too small, not enough for time evolutions. NORMALIZE EVERYTHING! for timesteps, want D to be 1. 
