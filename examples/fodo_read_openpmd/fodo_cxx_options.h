#ifndef FODO_CXX_OPTIONS_H_
#define FODO_CXX_OPTIONS_H_
#include <string>

struct Fodo_cxx_options {
    Fodo_cxx_options(int argc, char** argv);
    int gridx;
    int gridy;
    int gridz;
    int macroparticles;
    double real_particles;
    int turns;
    std::string input_openpmd_filename;
};
#endif /* FODO_CXX_OPTIONS_H_ */
