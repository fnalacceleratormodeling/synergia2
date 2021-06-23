
#include "synergia/foundation/trigon.h"
#include "synergia/bunch/bunch.h"
#include "synergia/simulation/lattice_simulator.h"
#include "synergia/libFF/ff_element.h"
#include "synergia/utils/lsexpr.h"
#include "synergia/utils/simple_timer.h"


void run()
{
    Logger screen(0, LoggerV::DEBUG);

    auto lsexpr = read_lsexpr_file("sis18-6.lsx");
    Lattice lattice(lsexpr);

    // tune the lattice
    lattice.set_all_string_attribute("extractor_type", "libff");
    Lattice_simulator::tune_circular_lattice(lattice);

    // get the reference particle
    auto const& ref = lattice.get_reference_particle();

    screen << "reference momentum = " 
        << ref.get_momentum() << " GeV\n";


#if 0
    using Tri = Trigon<double, 3, 3>;
    Tri x(0.0, 0);
    Tri y(0.0, 1);
    Tri z(0.0, 2);

    auto f = exp(x)*exp(y);
    auto g = kt::qpow(f, 3);
    auto k = x + y + x*y + y*z + z*z;

    screen << "x = " << x << "\n";
    screen << "y = " << y << "\n";
    screen << "f = exp(x)*exp(y) = " << f << "\n";
    screen << "g = qpow(f,3) = " << g << "\n";
    screen << "k = x+y+xy+yz+zz = " << k << "\n";


    arr_t<double, 3> v{0.1, 0.2, 1.0};
    std::cout << "f(v) = " << f(v) << "\n";
    std::cout << "g(v) = " << g(v) << "\n";
    std::cout << "k(v) = " << k(v) << "\n";

    return;
#endif

#if 0
    using Tri = Trigon<double, 3, 2>;
    Tri x(0.0, 0);
    Tri y(0.0, 1);

    screen << "x = " << x << "\n";
    screen << "y = " << y << "\n";

    TMapping<Tri> a, b;
    a[0] = x*y*y + exp(x+y);
    a[1] = cos(y*x*x) / (x+2.0);

    screen << "a[0] = " << a[0] << "\n";
    screen << "a[1] = " << a[1] << "\n";

    Tri xx(a[0].value(), 0);
    Tri yy(a[1].value(), 1);

    b[0] = sin(xx) * cos(yy);
    b[1] = exp(xx*xx*xx) / (xx*yy);

    screen << "b[0] = " << b[0] << "\n";
    screen << "b[1] = " << b[1] << "\n";

    auto c = b(a, {xx.value(), yy.value()});

    screen << "c[0] = " << c[0] << "\n";
    screen << "c[1] = " << c[1] << "\n";

    Tri q = a[0];
    Tri v = a[1];

    Tri w = sin(q)*cos(v);
    Tri z = exp(q*q*q) / (q*v);

    screen << "w = " << w << "\n";
    screen << "z = " << z << "\n";

    return;
#endif

    auto nf = Lattice_simulator::calculate_normal_form<2>(lattice);
    //auto ma = nf.stationaryActions(0.1, 0.1, 0.1);
    //std::cout << "ma = " << ma[0] << ", " << ma[1] << ", " << ma[2] << "\n";

    auto nfd = nf.cnvDataToNormalForm({0.1, 0.15, 0.2, 0.25, 0.05, 0.01});
    auto hfd = nf.cnvDataFromNormalForm(nfd);

    std::cout << "normal form = \n";
    for(int i=0; i<6; ++i) std::cout << nfd[i] << "\n";
    std::cout << "\n";

    std::cout << "human form = \n";
    for(int i=0; i<6; ++i) std::cout << hfd[i] << "\n";
    std::cout << "\n";

    return;
 

#if 0
    // tunes and cdt
    auto res = Lattice_simulator::calculate_tune_and_cdt(lattice, 0);

    std::cout << "tune_h = " << res[0] << "\n";
    std::cout << "tune_v = " << res[1] << "\n";
    std::cout << "c_delta_t = " << res[2] << "\n\n";

    // chromaticities
    double dpp = 0.01;
    auto chroms = Lattice_simulator::get_chromaticities(lattice, dpp);

    std::cout << "slip_factor = " << chroms.slip_factor << "\n";
    std::cout << "slip_factor_p = " << chroms.slip_factor_prime << "\n";
    std::cout << "momentum_compaction = " << chroms.momentum_compaction << "\n";
    std::cout << "h_chrom = " << chroms.horizontal_chromaticity << "\n";
    std::cout << "h_chrom_p = " << chroms.horizontal_chromaticity_prime << "\n";
    std::cout << "v_chrom = " << chroms.vertical_chromaticity << "\n";
    std::cout << "v_chrom_p = " << chroms.vertical_chromaticity_prime << "\n\n";
#endif

    // closed orbit
    auto probe = Lattice_simulator::calculate_closed_orbit(lattice, 0.01);

    std::cout << "closed orbit probe: ";
    for(int i=0; i<6; ++i) std::cout << probe[i] << ", ";
    std::cout << "\n";

    // trigon bunch
    using trigon_t = Trigon<double, 7, 6>;

    Commxx comm;
    bunch_t<trigon_t> tb(ref, comm.size(), comm);

    auto parts = tb.get_host_particles();

    // init value
    for(int i=0; i<6; ++i) 
        parts(0, i).set(probe[i], i);

    // check in
    tb.checkin_particles();

    // propagate trigon
    simple_timer_start("trigon propagate");

    for(auto & ele : lattice.get_elements())
    {
        if (ele.get_type() == element_type::rfcavity)
        {
            ele.set_double_attribute("volt", 0.0);
        }
        else if (ele.get_type() == element_type::quadrupole)
        {
            //ele.set_double_attribute("k1", 0.0);
            //ele.set_double_attribute("k1s", 0.0);
        }

        FF_element::apply(ele, tb);
    }

    simple_timer_stop("trigon propagate");
    std::cout << "\ntrigon propagated\n";

    tb.checkout_particles();
    auto jac = tb.get_jacobian(0);

    screen << std::setprecision(16);

    for(int i=0; i<6; ++i) {
        for(int j=0; j<6; ++j) {
            screen << jac(i, j) << ", ";
        }
        screen << "\n";
    }
    screen << "\n";

    simple_timer_print(screen);
    return;
}

int main(int argc, char ** argv)
{
    MPI_Init(&argc, &argv);
    Kokkos::initialize(argc, argv);

    run();

    Kokkos::finalize();
    MPI_Finalize();
    return 0;
}

