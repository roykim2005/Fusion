#include "Kadomtsev.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <stdexcept>
#include <algorithm>
#include <utility>

// finds r1 such that q(r1)=1
static double q_root_wrapper(double radial_r, const Kadomtsev* model) {
    return model->eval_q_minus_one(radial_r);
}
// finds r0 such that psi(r0)=psi(0)
static double psi_root_wrapper(double radial_r, const Kadomtsev* model) {
    return model->eval_psi_diff(radial_r);
}
// imports parameters from Kadomtsev.json
Kadomtsev::Kadomtsev(double qa_val, double pc_val, double epsa_val, double B0_val, double qc_val, int Nr_val)
    : qa(qa_val), pc(pc_val), epsa(epsa_val), B0(B0_val), qc(qc_val), Nr(Nr_val) {
    // sanity check
    validate_inputs();

    mu = qa / qc;
    // resizes arrays std::vector for Nr points for r, pre q, post q, B0, p, psi.
    allocate_memory();
}

void Kadomtsev::validate_inputs() const {
    if (Nr < 100) {
        throw std::invalid_argument("Nr must be at least 100.");
    }
    if (B0 <= 0.0) {
        throw std::invalid_argument("Toroidal magnetic field B0 must be positive.");
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
    int N_reconnect = 200; 

    // region r1 to 0
    for (int k = 0; k < N_reconnect; ++k) {
        double r1_minus = r1 * (1.0 - static_cast<double>(k) / N_reconnect);
        
        // finds r2_plus such that psi*(r1_minus) == psi*(r2_plus)
        double r2_plus = find_r2_outer(r1_minus);

        double r_post = std::sqrt(r2_plus * r2_plus - r1_minus * r1_minus);

        // Calculate derivative ratio dpsi/dr
        double dpsi_dr1 = std::abs(interpolate_1d(r, dpsi_star_dr, r1_minus));
        double dpsi_dr2 = std::abs(interpolate_1d(r, dpsi_star_dr, r2_plus));
        double dr1_dr2 = (dpsi_dr1 > 1e-9) ? (dpsi_dr2 / dpsi_dr1) : 1.0;

        double T1 = interpolate_1d(r, T_pre, r1_minus);
        double T2 = interpolate_1d(r, T_pre, r2_plus);

        // temp relaxation
        double T_mix = (T1 * r1_minus * dr1_dr2 + T2 * r2_plus) / (r1_minus * dr1_dr2 + r2_plus);

        r_post_list.push_back(r_post);
        T_mix_list.push_back(T_mix);
    }

    // for plot python plot generation, may or may not need this. currently it looks like its mapping things backwards
    std::vector<std::pair<double, double>> pairs(r_post_list.size());
    for (size_t i = 0; i < pairs.size(); ++i) {
        pairs[i] = {r_post_list[i], T_mix_list[i]};
    }
    std::sort(pairs.begin(), pairs.end());

    for (size_t i = 0; i < pairs.size(); ++i) {
        r_post_list[i] = pairs[i].first;
        T_mix_list[i]  = pairs[i].second;
    }

    for (int i = 0; i < Nr; ++i) {
        if (r[i] <= r0) {
            T_post[i] = interpolate_1d(r_post_list, T_mix_list, r[i]);
        } else {
            T_post[i] = T_pre[i]; // outer
        }

        // new resistivity
        eta_post[i] = 1.0e-6 * std::pow(T_post[i], -1.5);
    }
}


// sets up a 1D grid and calculates pre crash q, pre crash p, B0, d psi/dr = (1-q(r))B0, if T = 2 keV
void Kadomtsev::compute_profiles() {
    double dr = (1.0 - 1e-5) / (Nr - 1);
    
    // used Tc = 2 keV
    if (Tc <= 0.0) Tc = 2.0;

    for (int i = 0; i < Nr; ++i) {
        r[i] = 1e-5 + i * dr;
        
        q_pre[i] = qa * (r[i] * r[i]) / (1.0 - std::pow(1.0 - r[i] * r[i], mu));
        p[i] = pc * std::pow(1.0 - r[i] * r[i], mu);
        
        // Temperature profile
        T_pre[i] = Tc * std::pow(1.0 - r[i] * r[i], mu);
        
        // resistivity (spitzer, eta ~ T^(-3/2)), used 1e-6 ohm meter for eta_0
        eta_pre[i] = 1e-6 * std::pow(T_pre[i], -1.5);

        // B theta
        B_theta[i] = (r[i] * B0 * epsa) / q_pre[i];
        // helical flux
        dpsi_star_dr[i] = (1.0 - q_pre[i]) * B_theta[i];
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
    compute_profiles();
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

    std::cout << "T(0) - pre-crash: " << T_pre[0] << " keV, post-crash: " << T_post[0] << " keV\n";
    std::cout << "T(r1/2) - pre-crash: " << interpolate(r, T_pre, r1/2.0) << " keV, post-crash: " << interpolate(r, T_post, r1/2.0) << " keV\n";
    // T(r1) post is the temp after merging
    std::cout << "T(r1) - pre-crash: " << interpolate(r, T_pre, r1) << " keV, post-crash: " << interpolate(r, T_post, r1) << " keV\n";
    std::cout << "T(r0) - pre-crash: " << interpolate(r, T_pre, r0) << " keV, post-crash: " << interpolate(r, T_post, r0) << " keV\n\n";
    std::cout << std::scientific << std::setprecision(3);
    std::cout << "eta(0) - pre-crash: " << eta_pre[0] << " Ohm m, post-crash: " << eta_post[0] << " Ohm m\n";
    std::cout << "eta(r1) - pre-crash: " << interpolate(r, eta_pre, r1) << " Ohm m, post-crash: " << interpolate(r, eta_post, r1) << " Ohm m\n";
}



void Kadomtsev::write_data(const std::string& filename) const {
    std::ofstream file(filename);

    file << "# r/a\tq_pre\tq_post\tT_pre\tT_post\teta_pre\teta_post\tpsi_star\n";

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
