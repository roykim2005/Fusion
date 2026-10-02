#include "Kadomtsev.h"
#include "json.hpp"
#include <iostream>
#include <fstream>
#include <string>

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

        double qa = config["qa"].get<double>();
        double pc = config["pc"].get<double>();
        double epsa = config["epsa"].get<double>();
        double B0 = config["B0"].get<double>();
        double qc = config["qc"].get<double>();
        int Nr = config["Nr"].get<int>();

        double dt = config.value("dt", 1e-6);
        double t_final = config.value("t_final", 0.100); // in seconds
        double q_crit  = config.value("q_crit", 0.98); // instead of 1.0 


        Kadomtsev model(qa, pc, epsa, B0, qc, Nr);

        std::ofstream time_history("time_history.dat");
        if (!time_history.is_open()) {
            throw std::runtime_error("Could not create output file: time_history.dat");
        }
        time_history << "# t[s]\tT_center[keV]\tq_center\n";

        int crash_count = 0;

        // time evolution
        for (double t = 0.0; t < t_final; t += dt) {
            model.set_current_time(t);

            time_history << t << "\t" << model.get_T0() << "\t" << model.get_q0() << "\n";
            if (model.get_q0() < q_crit) {
                crash_count++;
                std::cout << "[Sawtooth crash #" << crash_count << "] at t = " 
                          << t * 1000.0 << " ms (q0 = " << model.get_q0() << ")\n";
                model.solve();
                model.print_summary();
                model.reset_post_to_pre();
            } else {
                model.step_ramp_phase(dt);
            }
            time_history << t << "\t" << model.get_T0() << "\t" << model.get_q0() << "\n";
        }

        time_history.close();

        std::cout << "\nTotal crashes recorded: " << crash_count << std::endl;
        
        model.write_data("kadomtsev.dat");

    } 
    catch (const std::exception& e) {
        std::cerr << "Program Kadomtsev - Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
