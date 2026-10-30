#include "Kadomtsev.h"
#include "json.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

using json = nlohmann::json;

int main(int argc, char* argv[]) {
    std::string config_filename = "Kadomtsev.json";
    if (argc > 1) {
        config_filename = argv[1];
    }

    try {
        json config;

        std::ifstream file(config_filename);
        if (!file.is_open()) {
            throw std::runtime_error("Could not open JSON configuration file: " + config_filename);
        }

        file >> config;
        double Z = config.value("Z", 2.0);
        double chi0 = config.value("chi0", 4.0);
        double qa = config["qa"].get<double>();
        double qc = config["qc"].get<double>();
        double epsa = config["epsa"].get<double>();
        double B0 = config["B0"].get<double>();
        double a = config["a"].get<double>(); 
        int Nr = config["Nr"].get<int>();

        double chi = chi0 / (a * a);
        double dt = config.value("dt", 1e-6); // s
        double t_final = config.value("t_final", 0.100); // s
        double q_crit = config.value("q_crit", 0.98);
        double n_e = config.value("n_e", 1.0e20); // m^-3
        double tau_R = config.value("tau_R", 0.1); // s
        double Tc = 2.0; // keV

        const double mu0 = 4.0e-7 * 3.14159265;
        const double e_keV = 1.60217663e-16; // J per kev
        const double m_e = 9.1093837e-31; // kg
        const double e_C = 1.602176634e-19; // C
        const double c_light = 2.99792458e8;
        const double pi = 3.14159265;

        double N_Z = 0.58 + 0.75 / (0.76 + Z);
        double alpha1 = 1.0 + 0.36 / Z;
        double ft_bar = 0.8 * (1.46 * std::sqrt(epsa) - (5.0 / 7.0) * 0.46 * std::pow(epsa, 1.5));
        double h_bar = 1.0 / (1.0 - alpha1 * ft_bar);
        double Bth_a = epsa * B0 / qa;
        double n20 = n_e / 1.0e20;

        auto lnLambda = [&](double T_keV) { return 15.2 - std::log(std::sqrt(n20) / T_keV); };
        double eta_coef = Z * N_Z * h_bar / (1.96 * 6.0 * std::sqrt(2.0) * std::pow(pi, 1.5)) * std::sqrt(m_e) * e_C * e_C * std::pow(c_light, 4) * mu0 * mu0;

        for (int it = 0; it < 100; ++it) {
            double T_J = std::pow(eta_coef * lnLambda(Tc) * Bth_a * Bth_a / (mu0 * mu0 * n_e * chi0), 0.4);
            double Tnew = T_J / e_keV;
            bool done = std::abs(Tnew - Tc) < 1e-10;
            Tc = Tnew;
            if (done) break;
        }
        double eta_c = eta_coef * lnLambda(Tc) * std::pow(Tc * e_keV, -1.5);

        double j_c = epsa * B0 / (mu0 * a); // density scaling


        double kappa = (2.0 / 3.0) * chi * t_final;
        double theta = (2.0 / 3.0) * t_final * eta_c * j_c * j_c / (n_e * e_keV * Tc);
        double tauR_hat = tau_R / t_final;
        double dt_hat = dt / t_final;
        double Tedge_hat = 0.05 / Tc;

        Kadomtsev model(qa, qc, Nr, kappa, theta, tauR_hat, Tedge_hat);

        std::ofstream time_history("time_history.dat");
        if (!time_history.is_open()) {
            throw std::runtime_error("Could not create output file: time_history.dat");
        }
        std::cout << "Scales: Tc = " << Tc << " keV, eta_c = " << eta_c << " ohm m\n";
        

        time_history << "# t_hat\tT_center/Tc\tq_center\n";

        int crash_count = 0;

        // time evolution, t_hat from 0 to 1
        long n_steps = std::llround(1.0 / dt_hat);
        for (long n = 0; n < n_steps; ++n) {
            double t_hat = n * dt_hat;
            double t_phys = t_hat * t_final; // seconds
            model.set_current_time(t_hat);

            if (model.get_q0() < q_crit) {
                crash_count++;
                std::cout << "[Sawtooth crash #" << crash_count << "] at t_hat = " << t_hat << " (q0 = " << model.get_q0() << ")\n";
                model.solve();
                model.print_summary();
                if (crash_count == 2) model.write_data("kadomtsev.dat");
                model.reset_post_to_pre();
            } else {
                model.step_ramp_phase(dt_hat);
            }

            time_history << t_hat << "\t" << model.get_T0() << "\t" << model.get_q0() << "\n";
        }

        time_history.close();

        std::cout << "\nTotal crashes recorded: " << crash_count << std::endl;
    }
    catch (const std::exception& e) {
        std::cerr << "Program Kadomtsev - Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
