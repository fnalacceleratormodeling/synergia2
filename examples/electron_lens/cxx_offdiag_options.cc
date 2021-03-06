#include "cxx_offdiag_options.h"
#include "synergia/utils/command_line_arg.h"

// this file was automatically generated by the command
//     synergia cxx_offdiag_options.py --create-cxx-options-source
// DO NOT EDIT

Cxx_offdiag_options::Cxx_offdiag_options(int argc, char **argv) :
    lattice_file("el_Xoffdiag.lsx"),
    rf_voltage(0.261977646285),
    elensradius(0.0041525),
    elensenergy(0.01),
    elenslength(1.0),
    elenscurrent(0.0),
    elenslongrms(0.5),
    elensdivide(2),
    elensadaptive(false),
    harmon(50),
    transport("maps"),
    bpms(true),
    stepper("splitoperator"),
    space_charge_3dh(1),
    space_charge_2dh(0),
    space_charge_rec(0),
    solver("3doh"),
    spc_comm_size(32),
    gridx(64),
    gridy(64),
    gridz(128),
    magiccomp(0),
    sccomp(1.0),
    steps_per_quad(1),
    num_steps_else(100),
    num_steps(72),
    map_order(1),
    xrms(0.0041525),
    yrms(0.0041525),
    zrms(0.5),
    dpoprms(0.001),
    k2l(0.0),
    beamparams(false),
    transhemit(1.09607333934e-06),
    transvemit(1.04757819876e-06),
    hbeta(18.3185523495),
    vbeta(18.0230840973),
    halpha(1.83173893189),
    valpha(-1.85824588289),
    flatbucket(false),
    emitx(1.00055865e-06),
    emity(1.00055865e-06),
    apertures(false),
    aperture_sigma(4.0),
    num_turns(1000),
    maxturns(3000),
    num_bunches(1),
    macroparticles(1000000),
    real_particles(200000000000.0),
    load_bunch(0),
    save_bunch(0),
    seed(13),
    turn_particles(1),
    turn_period(25),
    turn_track(10000),
    spc_tuneshift(1),
    cell_particles(false),
    xoffset(0.0),
    yoffset(0.0),
    zoffset(0.0),
    errsize(0.0),
    errelement(12),
    checkpointperiod(50),
    concurrentio(8),
    zparticle(false),
    set_nux(0.0),
    set_nuy(0.0),
    verbosity(1)
{
    for (int i = 1; i < argc; ++i) {
        Command_line_arg arg(argv[i]);
        if (arg.is_equal_pair()) {
            if (arg.get_lhs() == "lattice_file") {
                lattice_file = arg.extract_value<std::string >();
            } else if (arg.get_lhs() == "rf_voltage") {
                rf_voltage = arg.extract_value<double >();
            } else if (arg.get_lhs() == "elensradius") {
                elensradius = arg.extract_value<double >();
            } else if (arg.get_lhs() == "elensenergy") {
                elensenergy = arg.extract_value<double >();
            } else if (arg.get_lhs() == "elenslength") {
                elenslength = arg.extract_value<double >();
            } else if (arg.get_lhs() == "elenscurrent") {
                elenscurrent = arg.extract_value<double >();
            } else if (arg.get_lhs() == "elenslongrms") {
                elenslongrms = arg.extract_value<double >();
            } else if (arg.get_lhs() == "elensdivide") {
                elensdivide = arg.extract_value<int >();
            } else if (arg.get_lhs() == "elensadaptive") {
                elensadaptive = arg.extract_value<bool >();
            } else if (arg.get_lhs() == "harmon") {
                harmon = arg.extract_value<int >();
            } else if (arg.get_lhs() == "transport") {
                transport = arg.extract_value<std::string >();
            } else if (arg.get_lhs() == "bpms") {
                bpms = arg.extract_value<bool >();
            } else if (arg.get_lhs() == "stepper") {
                stepper = arg.extract_value<std::string >();
            } else if (arg.get_lhs() == "space_charge_3dh") {
                space_charge_3dh = arg.extract_value<int >();
            } else if (arg.get_lhs() == "space_charge_2dh") {
                space_charge_2dh = arg.extract_value<int >();
            } else if (arg.get_lhs() == "space_charge_rec") {
                space_charge_rec = arg.extract_value<int >();
            } else if (arg.get_lhs() == "solver") {
                solver = arg.extract_value<std::string >();
            } else if (arg.get_lhs() == "spc_comm_size") {
                spc_comm_size = arg.extract_value<int >();
            } else if (arg.get_lhs() == "gridx") {
                gridx = arg.extract_value<int >();
            } else if (arg.get_lhs() == "gridy") {
                gridy = arg.extract_value<int >();
            } else if (arg.get_lhs() == "gridz") {
                gridz = arg.extract_value<int >();
            } else if (arg.get_lhs() == "magiccomp") {
                magiccomp = arg.extract_value<int >();
            } else if (arg.get_lhs() == "sccomp") {
                sccomp = arg.extract_value<double >();
            } else if (arg.get_lhs() == "steps_per_quad") {
                steps_per_quad = arg.extract_value<int >();
            } else if (arg.get_lhs() == "num_steps_else") {
                num_steps_else = arg.extract_value<int >();
            } else if (arg.get_lhs() == "num_steps") {
                num_steps = arg.extract_value<int >();
            } else if (arg.get_lhs() == "map_order") {
                map_order = arg.extract_value<int >();
            } else if (arg.get_lhs() == "xrms") {
                xrms = arg.extract_value<double >();
            } else if (arg.get_lhs() == "yrms") {
                yrms = arg.extract_value<double >();
            } else if (arg.get_lhs() == "zrms") {
                zrms = arg.extract_value<double >();
            } else if (arg.get_lhs() == "dpoprms") {
                dpoprms = arg.extract_value<double >();
            } else if (arg.get_lhs() == "k2l") {
                k2l = arg.extract_value<double >();
            } else if (arg.get_lhs() == "beamparams") {
                beamparams = arg.extract_value<bool >();
            } else if (arg.get_lhs() == "transhemit") {
                transhemit = arg.extract_value<double >();
            } else if (arg.get_lhs() == "transvemit") {
                transvemit = arg.extract_value<double >();
            } else if (arg.get_lhs() == "hbeta") {
                hbeta = arg.extract_value<double >();
            } else if (arg.get_lhs() == "vbeta") {
                vbeta = arg.extract_value<double >();
            } else if (arg.get_lhs() == "halpha") {
                halpha = arg.extract_value<double >();
            } else if (arg.get_lhs() == "valpha") {
                valpha = arg.extract_value<double >();
            } else if (arg.get_lhs() == "flatbucket") {
                flatbucket = arg.extract_value<bool >();
            } else if (arg.get_lhs() == "emitx") {
                emitx = arg.extract_value<double >();
            } else if (arg.get_lhs() == "emity") {
                emity = arg.extract_value<double >();
            } else if (arg.get_lhs() == "apertures") {
                apertures = arg.extract_value<bool >();
            } else if (arg.get_lhs() == "aperture_sigma") {
                aperture_sigma = arg.extract_value<double >();
            } else if (arg.get_lhs() == "num_turns") {
                num_turns = arg.extract_value<int >();
            } else if (arg.get_lhs() == "maxturns") {
                maxturns = arg.extract_value<int >();
            } else if (arg.get_lhs() == "num_bunches") {
                num_bunches = arg.extract_value<int >();
            } else if (arg.get_lhs() == "macroparticles") {
                macroparticles = arg.extract_value<int >();
            } else if (arg.get_lhs() == "real_particles") {
                real_particles = arg.extract_value<double >();
            } else if (arg.get_lhs() == "load_bunch") {
                load_bunch = arg.extract_value<int >();
            } else if (arg.get_lhs() == "save_bunch") {
                save_bunch = arg.extract_value<int >();
            } else if (arg.get_lhs() == "seed") {
                seed = arg.extract_value<int >();
            } else if (arg.get_lhs() == "turn_particles") {
                turn_particles = arg.extract_value<int >();
            } else if (arg.get_lhs() == "turn_period") {
                turn_period = arg.extract_value<int >();
            } else if (arg.get_lhs() == "turn_track") {
                turn_track = arg.extract_value<int >();
            } else if (arg.get_lhs() == "spc_tuneshift") {
                spc_tuneshift = arg.extract_value<int >();
            } else if (arg.get_lhs() == "cell_particles") {
                cell_particles = arg.extract_value<bool >();
            } else if (arg.get_lhs() == "xoffset") {
                xoffset = arg.extract_value<double >();
            } else if (arg.get_lhs() == "yoffset") {
                yoffset = arg.extract_value<double >();
            } else if (arg.get_lhs() == "zoffset") {
                zoffset = arg.extract_value<double >();
            } else if (arg.get_lhs() == "errsize") {
                errsize = arg.extract_value<double >();
            } else if (arg.get_lhs() == "errelement") {
                errelement = arg.extract_value<int >();
            } else if (arg.get_lhs() == "checkpointperiod") {
                checkpointperiod = arg.extract_value<int >();
            } else if (arg.get_lhs() == "concurrentio") {
                concurrentio = arg.extract_value<int >();
            } else if (arg.get_lhs() == "zparticle") {
                zparticle = arg.extract_value<bool >();
            } else if (arg.get_lhs() == "set_nux") {
                set_nux = arg.extract_value<double >();
            } else if (arg.get_lhs() == "set_nuy") {
                set_nuy = arg.extract_value<double >();
            } else if (arg.get_lhs() == "verbosity") {
                verbosity = arg.extract_value<int >();
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
