#include "fodo_cxx_options.h"
#include "synergia/utils/command_line_arg.h"

// this file was automatically generated by the command
//     synergia fodo_cxx_options.py --create-cxx-options-source
// DO NOT EDIT

Fodo_cxx_options::Fodo_cxx_options(int argc, char **argv) :
    gridx(32),
    gridy(32),
    gridz(128),
    macroparticles(1048576),
    real_particles(2940000000000.0),
    turns(10)
{
    for (int i = 1; i < argc; ++i) {
        Command_line_arg arg(argv[i]);
        if (arg.is_equal_pair()) {
            if (arg.get_lhs() == "gridx") {
                gridx = arg.extract_value<int >();
            } else if (arg.get_lhs() == "gridy") {
                gridy = arg.extract_value<int >();
            } else if (arg.get_lhs() == "gridz") {
                gridz = arg.extract_value<int >();
            } else if (arg.get_lhs() == "macroparticles") {
                macroparticles = arg.extract_value<int >();
            } else if (arg.get_lhs() == "real_particles") {
                real_particles = arg.extract_value<double >();
            } else if (arg.get_lhs() == "turns") {
                turns = arg.extract_value<int >();
            } else if (arg.get_lhs() == "synergia_executable") {
                // ignore
            } else if (arg.get_lhs() == "run") {
                // ignore
            } else if (arg.get_lhs() == "submit") {
                // ignore
            } else if (arg.get_lhs() == "overwrite") {
                // ignore
            } else if (arg.get_lhs() == "numproc") {
                // ignore
            } else if (arg.get_lhs() == "procspernode") {
                // ignore
            } else if (arg.get_lhs() == "jobdir") {
                // ignore
            } else {
                throw std::runtime_error("Unknown argument " + arg.get_lhs());
            }
        } else {
            throw std::runtime_error("Bad argument " + std::string(argv[i]));
        }
    }
}