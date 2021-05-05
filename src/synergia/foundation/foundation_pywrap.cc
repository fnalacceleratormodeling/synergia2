

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "four_momentum.h"
#include "reference_particle.h"
#include "distribution.h"
#include "pcg_distribution.h"
#include "math_constants.h"
#include "physical_constants.h"
#include "trigon.h"


namespace py = pybind11;
using namespace py::literals;

template<unsigned int P, unsigned int D>
size_t 
canonical_to_index(std::array<size_t, array_length(P)> const& ind)
{
    static auto map = fill_index_to_canonical<P, D>();

    arr_t<size_t, array_length(P)> val;
    for(int i=0; i<val.size(); ++i) val[i] = ind[i];

    return map[val];
}

template<unsigned int P, unsigned int D>
void add_py_trigon_index_to_canonical(pybind11::module& m)
{
    std::stringstream ss;
    ss << "Trigon_index_to_canonical_o" << P;

    m.def(ss.str().c_str(), [](size_t idx) {
        static auto inds = indices<P, D>();
        std::array<size_t, array_length(P)> ret;
        for(int i=0; i<ret.size(); ++i) ret[i] = inds[idx][i];
        return ret;
    });
}

void add_py_trigon_canonical_to_index(pybind11::module& m)
{
    const char* name = "Trigon_canonical_to_index";

    m.def(name, [](std::array<size_t, array_length(1)> ind) {
        std::sort(ind.begin(), ind.end());
        return canonical_to_index<1, 6>(ind);
    });

    m.def(name, [](std::array<size_t, array_length(2)> ind) {
        std::sort(ind.begin(), ind.end());
        return canonical_to_index<2, 6>(ind);
    });

    m.def(name, [](std::array<size_t, array_length(3)> ind) {
        std::sort(ind.begin(), ind.end());
        return canonical_to_index<3, 6>(ind);
    });

    m.def(name, [](std::array<size_t, array_length(4)> ind) {
        std::sort(ind.begin(), ind.end());
        return canonical_to_index<4, 6>(ind);
    });

    m.def(name, [](std::array<size_t, array_length(5)> ind) {
        std::sort(ind.begin(), ind.end());
        return canonical_to_index<5, 6>(ind);
    });

    m.def(name, [](std::array<size_t, array_length(6)> ind) {
        std::sort(ind.begin(), ind.end());
        return canonical_to_index<6, 6>(ind);
    });

    m.def(name, [](std::array<size_t, array_length(7)> ind) {
        std::sort(ind.begin(), ind.end());
        return canonical_to_index<7, 6>(ind);
    });
}


template<class T, unsigned int P, unsigned int D>
void add_py_trigon_class(pybind11::module& m, const char* name)
{
    using trigon_t = Trigon<T, P, D>;

    py::class_<trigon_t>(m, name)
        .def( "count", 
                [](trigon_t& self) { return trigon_t::count; } )
        .def( "dim", 
                [](trigon_t& self) { return trigon_t::dim; } )
        .def( "power", 
                [](trigon_t& self) { return trigon_t::power(); } )

        .def( py::init<>(), "Construct" )
        .def( py::init<T>(), "Construct with a val", "val"_a )
        .def( "value", 
                (T const& (trigon_t::*)() const)&trigon_t::value)
        .def( "set", 
                (void (trigon_t::*)(T))&trigon_t::set, 
                "val"_a )
        .def( "set", 
                (void (trigon_t::*)(T, size_t))&trigon_t::set, 
                "val"_a, "index"_a )

        .def( "get_term",
                (T (trigon_t::*)(size_t))&trigon_t::get_term,
                "idx"_a )

        .def( "get_term",
                [](trigon_t& self, std::array<size_t, array_length(P)> const& ind) { return self.get_term(canonical_to_index<P, D>(ind)); },
                "indices"_a )
            

        .def( "get_term",
                (T (trigon_t::*)(unsigned int, size_t))&trigon_t::get_term,
                "power"_a, "idx"_a )


        .def( "set_term",
                (void (trigon_t::*)(size_t, T const&))&trigon_t::set_term,
                "idx"_a, "val"_a )

        .def( "set_term",
                (void (trigon_t::*)(unsigned int, size_t, T const&))&trigon_t::set_term,
                "power"_a, "idx"_a, "val"_a )

        .def( "print_coeffs",
                [](trigon_t& self) { std::cout << self; } )


        ;
}

PYBIND11_MODULE(foundation, m)
{
    auto pc = m.def_submodule("pconstants");
    pc.attr("pdg_year")        = pconstants::pdg_year;
    pc.attr("mp")              = pconstants::mp;
    pc.attr("proton_mass")     = pconstants::proton_mass;
    pc.attr("me")              = pconstants::me;
    pc.attr("electron_mass")   = pconstants::electron_mass;
    pc.attr("mmu")             = pconstants::mmu;
    pc.attr("muon_mass")       = pconstants::muon_mass;
    pc.attr("e")               = pconstants::e;
    pc.attr("c")               = pconstants::c;
    pc.attr("mu0")             = pconstants::mu0;
    pc.attr("epsilon0")        = pconstants::epsilon0;
    pc.attr("re")              = pconstants::re;
    pc.attr("rp")              = pconstants::rp;
    pc.attr("rmu")             = pconstants::rmu;
    pc.attr("proton_charge")   = pconstants::proton_charge;
    pc.attr("antiproton_charge") = pconstants::antiproton_charge;
    pc.attr("electron_charge") = pconstants::electron_charge;
    pc.attr("positron_charge") = pconstants::positron_charge;
    pc.attr("muon_charge")     = pconstants::muon_charge;
    pc.attr("antimuon_charge") = pconstants::antimuon_charge;
 
    // Four_momentum
    py::class_<Four_momentum>(m, "Four_momentum")
        .def( py::init<double>(),
                "Construct a Four_momentum in the rest frame."
                "\n\tparam: mass in GeV/c^2",
                "mass"_a )

        .def( py::init<double, double>(),
                "Construct a Four_momentum with the given total energy."
                "\n\tparam: mass in GeV/c^2"
                "\n\tparam: total_energy in GeV",
                "mass"_a, "total_energy"_a )

        .def( "set_total_energy",
                &Four_momentum::set_total_energy,
                "Set the total energy (in GeV)",
                "total_energy"_a )


        .def( "set_kinetic_energy",
                &Four_momentum::set_kinetic_energy,
                "Set the kinetic energy (in GeV)",
                "kinetic_energy"_a )

        .def( "set_momentum",
                &Four_momentum::set_momentum,
                "Set the momentum (in GeV/c)",
                "momentum"_a )

        .def( "set_gamma",
                &Four_momentum::set_gamma,
                "Set the relativistic gamma factor (unitless)",
                "gamma"_a )

        .def( "set_beta",
                &Four_momentum::set_beta,
                "Set the relativistic beta factor (unitless)",
                "beta"_a )

        .def( "get_mass",
                &Four_momentum::get_mass,
                "Get the mass in GeV/c^2" )

        .def( "get_total_energy",
                &Four_momentum::get_total_energy,
                "Get the total energy in GeV" )

        .def( "get_kinetic_energy",
                &Four_momentum::get_kinetic_energy,
                "Get the kinetic energy in GeV" )

        .def( "get_momentum",
                &Four_momentum::get_momentum,
                "Get the momentum in GeV/c" )

        .def( "get_gamma",
                &Four_momentum::get_gamma,
                "Get the relativistic gamma factor" )

        .def( "get_beta",
                &Four_momentum::get_beta,
                "Get the relativistic beta factor" )

        .def( "equal",
                &Four_momentum::equal,
                "Check equality to the given tolerance"
                "\n\tparam: tolerance factional accuracy for beta and gamma",
                "four_momentum"_a, "tolerance"_a )
        ;

    // Reference_particle
    py::class_<Reference_particle>(m, "Reference_particle")
        .def( py::init<int, double, double>(),
                "Construct a Reference_particle with a given mass and total energy."
                "\n\tparam: charge in units of e"
                "\n\tparam: mass in GeV/c^2"
                "\n\tparam: total_energy in GeV in the lab frame",
                "charge"_a, "mass"_a, "total_energy"_a )

        .def( py::init<int, Four_momentum const&>(),
                "Construct a Reference_particle with a given four momentum"
                "\n\tparam: charge in units of e"
                "\n\tparam: four_momentum in the lab frame",
                "charge"_a, "four_momentum"_a )

        .def( py::init<int, Four_momentum const&, std::array<double, 6> const&>(),
                "Construct a Reference_particle with a given four momentum and state"
                " in the reference frame."
                "\n\tparam: charge in units of e"
                "\n\tparam: four_momentum in the lab frame"
                "\n\tparam: state is a six-dimensional state vector",
                "charge"_a, "four_momentum"_a, "state"_a )

        .def( "set_four_momentum",
                &Reference_particle::set_four_momentum,
                "Set the four momentum."
                "four_momentum"_a )

        .def( "set_state", 
                (void (Reference_particle::*)(std::array<double, 6> const&))&Reference_particle::set_state,
                //py::overload_cast<std::array<double, 6> const&>(&Reference_particle::set_state),
                "Set the state vector in the reference frame.",
                "state"_a )

        .def( "set_state",
                (void (Reference_particle::*)(double, double, double, double, double, double))&Reference_particle::set_state,
                //py::overload_cast<double, double, double, double, double, double>(&Reference_particle::set_state),
                "Set the state vector in the reference frame.",
                "x"_a, "xp"_a, "y"_a, "yp"_a, "cdt"_a, "dpop"_a )

        .def( "set_total_energy",
                &Reference_particle::set_total_energy,
                "Set the total energy."
                "\n\tparam: total_energy in GeV in the lab frame",
                "total_energy"_a )

        .def( "increment_trajectory", 
                &Reference_particle::increment_trajectory,
                "Increment the trajectory length."
                "\n\tparam: length in m",
                "length"_a )

        .def( "start_repetition", 
                &Reference_particle::start_repetition,
                "Start a new repetition." )

        .def( "set_trajectory", 
                &Reference_particle::set_trajectory,
                "Manually set the trajectory parameters."
                "\n\tparam: repetition starting at 0"
                "\n\tparam: repetition_length in m"
                "\n\tparam: s in m",
                "repetition"_a, "repetition_length"_a, "s"_a )

        .def( "get_charge", 
                &Reference_particle::get_charge,
                "Return the Reference_particle charge in units of e." )

        .def( "get_mass", 
                &Reference_particle::get_mass,
                "Return the Reference_particle mass in units of GeV/c." )

        .def( "get_four_momentum",
                &Reference_particle::get_four_momentum,
                "Get the four momentum in the lab frame." )

        .def( "get_state",
                &Reference_particle::get_state,
                "Get the six-dimensional state vector the reference frame." )

        .def( "get_beta",
                &Reference_particle::get_beta,
                "Get the relativistic beta in the lab frame." )

        .def( "get_gamma",
                &Reference_particle::get_gamma,
                "Get the relativistic gamma in the lab frame." )

        .def( "get_momentum",
                &Reference_particle::get_momentum,
                "Get the momentum in GeV/c in the lab frame." )

        .def( "get_total_energy",
                &Reference_particle::get_total_energy,
                "Get the total energy in GeV in the lab frame." )

        .def( "get_s", 
                &Reference_particle::get_s,
                "Get the total path length in m of the reference particle trajectory." )

        .def( "get_s_n", 
                &Reference_particle::get_s_n,
                "Get the distance traveled in m since the beginning of the current repetition." )

        .def( "get_repetition", 
                &Reference_particle::get_repetition,
                "Get the number of repetition." )

        .def( "get_repetition_length", 
                &Reference_particle::get_repetition_length,
                "Get the repetition length in m." )

        .def( "equal",
                &Reference_particle::equal,
                "Check equality to the given tolerance."
                "\n\tparam: reference_particle another Reference_particle"
                "\n\tparam: tolerance fractional accuracy",
                "reference_particle"_a, "tolerance"_a )
        ;

    py::class_<Distribution>(m, "Distribution")
        ;

    py::class_<Random_distribution, Distribution>(m, "Random_distribution")
        .def( py::init<unsigned long int, Commxx const&>(),
                "seed"_a, "comm"_a )
        ;

    py::class_<PCG_random_distribution, Distribution>(m, "PCG_random_distribution")
        .def( py::init<unsigned long int, Commxx const&>(),
                "seed"_a, "comm"_a )
        ;

    // Trigon related classes and functions
    add_py_trigon_index_to_canonical<1, 6>(m);
    add_py_trigon_index_to_canonical<2, 6>(m);
    add_py_trigon_index_to_canonical<3, 6>(m);
    add_py_trigon_index_to_canonical<4, 6>(m);
    add_py_trigon_index_to_canonical<5, 6>(m);
    add_py_trigon_index_to_canonical<6, 6>(m);
    add_py_trigon_index_to_canonical<7, 6>(m);

    add_py_trigon_canonical_to_index(m);

    add_py_trigon_class<double, 1, 6>(m, "Trigon_o1");
    add_py_trigon_class<double, 2, 6>(m, "Trigon_o2");
    add_py_trigon_class<double, 3, 6>(m, "Trigon_o3");
    add_py_trigon_class<double, 4, 6>(m, "Trigon_o4");
    add_py_trigon_class<double, 5, 6>(m, "Trigon_o5");
    add_py_trigon_class<double, 6, 6>(m, "Trigon_o6");
    add_py_trigon_class<double, 7, 6>(m, "Trigon_o7");


#if 0
    class_<Distribution, boost::noncopyable > ("Distribution", no_init);

    class_<Dummy1 >("mconstants",no_init)
        .def_readonly("pi", mconstants::pi)
        ;

    class_<Dummy2 >("pconstants",no_init)
        .def_readonly("pdg_year", pconstants::pdg_year)
        .def_readonly("mp", pconstants::mp)
        .def_readonly("proton_mass", pconstants::proton_mass)
        .def_readonly("me", pconstants::me)
        .def_readonly("electron_mass", pconstants::electron_mass)
        .def_readonly("mmu", pconstants::mmu)
        .def_readonly("muon_mass", pconstants::muon_mass)
        .def_readonly("e", pconstants::e)
        .def_readonly("c", pconstants::c)
        .def_readonly("mu0", pconstants::mu0)
        .def_readonly("epsilon0", pconstants::epsilon0)
        .def_readonly("re", pconstants::re)
        .def_readonly("rp", pconstants::rp)
        .def_readonly("rmu", pconstants::rmu)
        .def_readonly("proton_charge", pconstants::proton_charge)
        .def_readonly("antiproton_charge", pconstants::antiproton_charge)
        .def_readonly("electron_charge", pconstants::electron_charge)
        .def_readonly("positron_charge", pconstants::positron_charge)
        .def_readonly("muon_charge", pconstants::muon_charge)
        .def_readonly("antimuon_charge", pconstants::antimuon_charge)
        ;

    //    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(get_default_seed_overloads, get_default_seed, 0, 1);
    scope
    Random_distribution_scope =
        class_<Random_distribution, bases<Distribution > > ("Random_distribution",
                init<unsigned long int, Commxx const &>())
            .def(init<unsigned long int, Commxx const &, Random_distribution::Generator>())
            .def("get_default_seed",&Random_distribution::get_default_seed)
            .staticmethod("get_default_seed")
            .def("set_seed",&Random_distribution::set_seed)
            .def("get_original_seed",&Random_distribution::get_original_seed)
//            .def("get_original_seed",(void(*)(const char*))0,get_default_seed_overloads())
            .def("get",&Random_distribution::get)
            .def("fill_uniform",random_fill_uniform_ref)
            .def("fill_unit_gaussian",random_fill_unit_gaussian_ref)
            .def("fill_unit_disk",random_fill_unit_disk_ref)
            ;

        enum_<Random_distribution::Generator >("Generator")
            .value("ranlxd2", Random_distribution::ranlxd2)
            .value("mt19937", Random_distribution::mt19937)
            .export_values();

#endif
}
