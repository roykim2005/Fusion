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
            throw std::runtime_error("Could not open JSON file: " + config_filename);
        }

        file >> config;

        double qa   = config["qa"].get<double>();
        double pc   = config["pc"].get<double>();
        double epsa = config["epsa"].get<double>();
        double B0   = config["B0"].get<double>();
        double qc   = config["qc"].get<double>();
        int Nr      = config["Nr"].get<int>();

        Kadomtsev model(qa, pc, epsa, B0, qc, Nr);
        model.solve();
        model.print_summary();
        model.write_data("kadomtsev.dat");

    } 
    catch (const std::exception& e) {
        std::cerr << "Program Kadomtsev - Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
