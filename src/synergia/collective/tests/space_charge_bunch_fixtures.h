#ifndef SPACE_CHARGE_BUNCH_FIXTURES_H_
#define SPACE_CHARGE_BUNCH_FIXTURES_H_

const int charge = pconstants::proton_charge;
const double mass = pconstants::mp;
const double real_num = 1.7e11;
const int total_num = 10000;
const double total_energy = 125.0;
struct Ellipsoidal_bunch_fixture
{
    Ellipsoidal_bunch_fixture() :
        four_momentum(mass, total_energy), reference_particle(charge,
                four_momentum), comm_sptr(new Commxx), bunch(reference_particle,
                total_num, real_num, comm_sptr), seed(718281828), distribution(seed, *comm_sptr),
                grid_shape(3)
    {
        BOOST_TEST_MESSAGE("setup ellipsoidal bunch fixture");
        MArray2d covariances(boost::extents[6][6]);
        MArray1d means(boost::extents[6]);
        for (int i = 0; i < 6; ++i) {
            means[i] = 0.0;
            for (int j = i; j < 6; ++j) {
                covariances[i][j] = 0.0;
            }
        }
        stdx = 1.1e-3;
        stdy = 2.3e-3;
        stdz = 3.5e-3;
        covariances[0][0] = stdx * stdx;
        covariances[2][2] = stdy * stdy;
        covariances[4][4] = stdz * stdz;
        covariances[1][1] = covariances[3][3] = covariances[5][5] = 0.00001;
        populate_6d(distribution, bunch, means, covariances);
        grid_shape[0] = 16;
        grid_shape[1] = 24;
        grid_shape[2] = 32;
    }

    ~Ellipsoidal_bunch_fixture()
    {
        BOOST_TEST_MESSAGE("tear down ellipsoidal bunch fixture");
    }

    Four_momentum four_momentum;
    Reference_particle reference_particle;
    Commxx_sptr comm_sptr;
    Bunch bunch;
    unsigned long int seed;
    Random_distribution distribution;
    double stdx, stdy, stdz;
    std::vector<int > grid_shape;
};

struct Spherical_bunch_fixture
{
    Spherical_bunch_fixture() :
        four_momentum(mass, total_energy), reference_particle(charge,
                four_momentum), comm_sptr(new Commxx), bunch(reference_particle,
                total_num, real_num, comm_sptr), seed(718281828), distribution(seed, *comm_sptr),
                grid_shape(3)
    {
        BOOST_TEST_MESSAGE("setup Spherical bunch fixture");
        MArray2d covariances(boost::extents[6][6]);
        MArray1d means(boost::extents[6]);
        for (int i = 0; i < 6; ++i) {
            means[i] = 0.0;
            for (int j = i; j < 6; ++j) {
                covariances[i][j] = 0.0;
            }
        }
        sigma = 1.3e-3;
        covariances[0][0] = sigma * sigma;
        covariances[2][2] = sigma * sigma;
        covariances[4][4] = sigma * sigma;
        covariances[1][1] = covariances[3][3] = covariances[5][5] = 0.00001;
        populate_6d(distribution, bunch, means, covariances);
        grid_shape[0] = 16;
        grid_shape[1] = 24;
        grid_shape[2] = 32;
    }

    ~Spherical_bunch_fixture()
    {
        BOOST_TEST_MESSAGE("tear down Spherical bunch fixture");
    }

    Four_momentum four_momentum;
    Reference_particle reference_particle;
    Commxx_sptr comm_sptr;
    Bunch bunch;
    unsigned long int seed;
    Random_distribution distribution;
    double sigma;
    std::vector<int > grid_shape;
};

struct Spherical_bunch_fixture_2d
{
    Spherical_bunch_fixture_2d() :
                    four_momentum(mass, total_energy),
                    reference_particle(charge, four_momentum),
                    comm_sptr(new Commxx),
                    bunch(reference_particle, total_num, real_num, comm_sptr),
                    seed(718281828),
                    distribution(seed, *comm_sptr),
                    grid_shape(3)
    {
        BOOST_TEST_MESSAGE("setup Spherical bunch fixture");
        MArray2d covariances(boost::extents[6][6]);
        MArray1d means(boost::extents[6]);
        for (int i = 0; i < 6; ++i) {
            means[i] = 0.0;
            for (int j = i; j < 6; ++j) {
                covariances[i][j] = 0.0;
            }
        }
        sigma = 1.3e-3;
        sigmaz = 1.3e-1;
        covariances[0][0] = sigma * sigma;
        covariances[2][2] = sigma * sigma;
        covariances[4][4] = sigmaz * sigmaz;
        covariances[1][1] = covariances[3][3] = covariances[5][5] = 0.00001;
        populate_6d(distribution, bunch, means, covariances);
        grid_shape[0] = 16;
        grid_shape[1] = 24;
        grid_shape[2] = 32;
    }

    ~Spherical_bunch_fixture_2d()
    {
        BOOST_TEST_MESSAGE("tear down Spherical bunch fixture");
    }

    Four_momentum four_momentum;
    Reference_particle reference_particle;
    Commxx_sptr comm_sptr;
    Bunch bunch;
    unsigned long int seed;
    Random_distribution distribution;
    double sigma, sigmaz;
    std::vector<int > grid_shape;
};

struct Cylindrical_bunch_fixture
{
    Cylindrical_bunch_fixture() :
        four_momentum(mass, total_energy), reference_particle(charge,
                four_momentum), comm_sptr(new Commxx), bunch(reference_particle,
                total_num, real_num, comm_sptr), seed(718281828), distribution(seed, *comm_sptr),
                grid_shape(3)
    {
        BOOST_TEST_MESSAGE("setup Cylindrical bunch fixture");
        MArray2d covariances(boost::extents[6][6]);
        MArray1d means(boost::extents[6]);
        for (int i = 0; i < 6; ++i) {
            means[i] = 0.0;
            for (int j = i; j < 6; ++j) {
                covariances[i][j] = 0.0;
            }
        }
        sigma = 1.3e-3;
        covariances[0][0] = sigma * sigma;
        covariances[2][2] = sigma * sigma;
        covariances[4][4] = sigma * sigma * 1000;
        covariances[1][1] = covariances[3][3] = covariances[5][5] = 0.00001;
        populate_6d(distribution, bunch, means, covariances);
        grid_shape[0] = 16;
        grid_shape[1] = 16;
        grid_shape[2] = 16;
    }

    ~Cylindrical_bunch_fixture()
    {
        BOOST_TEST_MESSAGE("tear down Cylindrical bunch fixture");
    }

    Four_momentum four_momentum;
    Reference_particle reference_particle;
    Commxx_sptr comm_sptr;
    Bunch bunch;
    unsigned long int seed;
    Random_distribution distribution;
    double sigma;
    std::vector<int > grid_shape;
};

const int fine_num_particles = 20000;
struct Cylindrical_bunch_fixture_fine
{
    Cylindrical_bunch_fixture_fine() :
        four_momentum(mass, total_energy), reference_particle(charge,
                four_momentum), comm_sptr(new Commxx), bunch(reference_particle,
                fine_num_particles, real_num,
                            comm_sptr),
                    seed(718281828),
                    distribution(seed, *comm_sptr),
                grid_shape(3)
    {
        BOOST_TEST_MESSAGE("setup Cylindrical bunch fixture");
        MArray2d covariances(boost::extents[6][6]);
        MArray1d means(boost::extents[6]);
        for (int i = 0; i < 6; ++i) {
            means[i] = 0.0;
            for (int j = i; j < 6; ++j) {
                covariances[i][j] = 0.0;
            }
        }
        sigma = 1.3e-3;
        covariances[0][0] = sigma * sigma;
        covariances[2][2] = sigma * sigma;
        covariances[4][4] = sigma * sigma;
        covariances[1][1] = covariances[3][3] = covariances[5][5] = 0.00001;
        populate_6d(distribution, bunch, means, covariances);
        grid_shape[0] = 64;
        grid_shape[1] = 64;
        grid_shape[2] = 16;
    }

    ~Cylindrical_bunch_fixture_fine()
    {
        BOOST_TEST_MESSAGE("tear down Cylindrical bunch fixture");
    }

    Four_momentum four_momentum;
    Reference_particle reference_particle;
    Commxx_sptr comm_sptr;
    Bunch bunch;
    unsigned long int seed;
    Random_distribution distribution;
    double sigma;
    std::vector<int > grid_shape;
};


// ----------------------------------------------------------------------------

const int rod_num_particles = 1 + 2501*8+2; // 1 + (odd number) * 8 + 2 to expand the domein on both sides
const double rod_betagamma = 11.0/60.0;
const double rod_gamma = 61.0/60.0;
const double rod_beta = 11.0/61.0;
const double rod_real_num = 5.0e9;
const double rod_length = 0.1;
const double rod_radius = 1.0e-6; //radius of rod particles
//const double rod_radius = 1.0e-3; //radius of rod particles
const double rod_probe = 1.0e-3;  // radius of the probe particle
struct Rod_bunch_fixture
{
    Rod_bunch_fixture() :
        four_momentum(mass, mass*rod_gamma), reference_particle(charge,
                four_momentum), comm_sptr(new Commxx), bunch(reference_particle,
                rod_num_particles, rod_real_num,
                            comm_sptr, rod_length),
                grid_shape(3)
    {
        BOOST_TEST_MESSAGE("setup Rod bunch fixture");
        bunch.set_sort_period(-1);
        MArray2d_ref local_particles(bunch.get_local_particles());
        // a ring of 8 particles around each longitudinal location
        int num_longitudinal = (rod_num_particles-3)/8;
        double dz = rod_length/(num_longitudinal-1);

        double r2o2 = std::sqrt(2.0)/2.0;

        double z = -rod_length/2.0;
        for (int i=3; i<rod_num_particles; i+=8, z+=dz) {
            local_particles[i][Bunch::x] = rod_radius;
            local_particles[i][Bunch::y] = 0.0;

            local_particles[i+1][Bunch::x] = rod_radius*r2o2;
            local_particles[i+1][Bunch::y] = 2*rod_radius*r2o2;

            local_particles[i+2][Bunch::x] = 0.0;
            local_particles[i+2][Bunch::y] = 2*rod_radius;

            local_particles[i+3][Bunch::x] = -rod_radius*r2o2;
            local_particles[i+3][Bunch::y] =  2*rod_radius*r2o2;

            local_particles[i+4][Bunch::x] = -rod_radius;
            local_particles[i+4][Bunch::y] = 0.0;

            local_particles[i+5][Bunch::x] = -rod_radius*r2o2;
            local_particles[i+5][Bunch::y] = -2*rod_radius*r2o2;

            local_particles[i+6][Bunch::x] = 0.0;
            local_particles[i+6][Bunch::y] = -2*rod_radius;

            local_particles[i+7][Bunch::x] = rod_radius*r2o2;
            local_particles[i+7][Bunch::y] = -2*rod_radius*r2o2;

            for (int j=i; j<i+8; ++j) {
                local_particles[j][Bunch::cdt] = z/rod_beta;
                local_particles[j][Bunch::xp] = 0.0;
                local_particles[j][Bunch::yp] = 0.0;
                local_particles[j][Bunch::dpop] = 0.0;
                local_particles[j][Bunch::id] = j;
            }
        }
        local_particles[0][Bunch::x] = rod_probe;
        local_particles[0][Bunch::y] = 0.0;
        local_particles[0][Bunch::xp] = 0.0;
        local_particles[0][Bunch::yp] = 0.0;
        local_particles[0][Bunch::cdt] = 0.0;
        local_particles[0][Bunch::dpop] = 0.0;
        local_particles[0][Bunch::id] = 0.0;

        local_particles[1][Bunch::x] = rod_probe*1.2;
        local_particles[1][Bunch::y] = rod_probe*1.2;
        local_particles[1][Bunch::xp] = 0.0;
        local_particles[1][Bunch::yp] = 0.0;
        local_particles[1][Bunch::cdt] = 0.0;
        local_particles[1][Bunch::dpop] = 0.0;
        local_particles[1][Bunch::id] = 1.0;

        local_particles[2][Bunch::x] = -rod_probe*1.2;
        local_particles[2][Bunch::y] = -rod_probe*1.2;
        local_particles[2][Bunch::xp] = 0.0;
        local_particles[2][Bunch::yp] = 0.0;
        local_particles[2][Bunch::cdt] = 0.0;
        local_particles[2][Bunch::dpop] = 0.0;
        local_particles[2][Bunch::id] = 2.0;

//        grid_shape[0] = 64;
//        grid_shape[1] = 64;
//        grid_shape[2] = 128;
        grid_shape[0] = 4;
        grid_shape[1] = 8;
        grid_shape[2] = 12;
    }

    ~Rod_bunch_fixture()
    {
        BOOST_TEST_MESSAGE("tear down Rod bunch fixture");
    }

    Four_momentum four_momentum;
    Reference_particle reference_particle;
    Commxx_sptr comm_sptr;
    Bunch bunch;
    unsigned long int seed;
    std::vector<int > grid_shape;
};












// ----------------------------------------------------------------------------

const double domain_min = -2.0;
const double domain_max = 2.0;
const double domain_offset = 0.0;
const int toy_grid_shape[] = { 4, 6, 8 };
const int toy_total_num = 1;
const double toy_real_num = 1.0e20;

struct Toy_bunch_fixture
{
    Toy_bunch_fixture() :
                four_momentum(mass, total_energy),
                reference_particle(pconstants::proton_charge, four_momentum),
                comm_sptr(new Commxx),
                bunch(reference_particle, total_num, toy_real_num, comm_sptr),
                physical_size(3),
                physical_offset(3),
                cell_size(3),
                grid_shape(3),
                expected(
                        boost::extents[toy_grid_shape[2]][toy_grid_shape[1]][toy_grid_shape[0]])

    {
        for (int i = 0; i < 3; ++i) {
            grid_shape[i] = toy_grid_shape[i];
            physical_offset[i] = domain_offset;
            physical_size[i] = domain_max - domain_min;
            cell_size[i] = (domain_max - domain_min) / grid_shape[i];
        }
        for (unsigned int i = 0; i < expected.shape()[0]; ++i) {
            for (unsigned int j = 0; j < expected.shape()[1]; ++j) {
                for (unsigned int k = 0; k < expected.shape()[2]; ++k) {
                    expected[i][j][k] = 0.0;
                }
            }
        }
        density_norm = (toy_real_num / total_num) * pconstants::e
                / (cell_size[0] * cell_size[1] * cell_size[2]);
    }

    ~Toy_bunch_fixture()
    {
    }

    Four_momentum four_momentum;
    Reference_particle reference_particle;
    Commxx_sptr comm_sptr;
    Bunch bunch;
    double density_norm;

    std::vector<double > physical_size, physical_offset, cell_size;
    std::vector<int > grid_shape;
    Rectangular_grid_sptr rho_grid_sptr;
    MArray3d expected;
};

const int toy_grid_shape_xyz[] = { 4, 6, 8 };
struct Toy_bunch_fixture_xyz
{
    Toy_bunch_fixture_xyz() :
                four_momentum(mass, total_energy),
                reference_particle(pconstants::proton_charge, four_momentum),
                comm_sptr(new Commxx),
                bunch(reference_particle, total_num, toy_real_num, comm_sptr),
                physical_size(3),
                physical_offset(3),
                cell_size(3),
                grid_shape(3),
                expected(
                        boost::extents[toy_grid_shape_xyz[0]][toy_grid_shape_xyz[1]][toy_grid_shape_xyz[2]])

    {
        for (int i = 0; i < 3; ++i) {
            grid_shape[i] = toy_grid_shape_xyz[i];
            physical_offset[i] = domain_offset;
            physical_size[i] = domain_max - domain_min;
            cell_size[i] = (domain_max - domain_min) / grid_shape[i];
        }
        for (unsigned int i = 0; i < expected.shape()[0]; ++i) {
            for (unsigned int j = 0; j < expected.shape()[1]; ++j) {
                for (unsigned int k = 0; k < expected.shape()[2]; ++k) {
                    expected[i][j][k] = 0.0;
                }
            }
        }
        density_norm = (toy_real_num / total_num) * pconstants::e
                / (cell_size[0] * cell_size[1] * cell_size[2]);
    }

    ~Toy_bunch_fixture_xyz()
    {
    }

    Four_momentum four_momentum;
    Reference_particle reference_particle;
    Commxx_sptr comm_sptr;
    Bunch bunch;
    double density_norm;

    std::vector<double > physical_size, physical_offset, cell_size;
    std::vector<int > grid_shape;
    Rectangular_grid_sptr rho_grid_sptr;
    MArray3d expected;
};


struct Toy_bunch_fixture_2d
{
    Toy_bunch_fixture_2d() :
                four_momentum(mass, total_energy),
                reference_particle(pconstants::proton_charge, four_momentum),
                comm_sptr(new Commxx),
                bunch(reference_particle, total_num, toy_real_num, comm_sptr),
                physical_size(3),
                physical_offset(3),
                cell_size(3),
                grid_shape(3),
                expected_2dc(boost::extents[2*toy_grid_shape[0]]
                        [2*toy_grid_shape[1]]),
                expected_1d(boost::extents[toy_grid_shape[2]])
    {
        for (int i = 0; i < 3; ++i) {
            grid_shape[i] = toy_grid_shape[i];
            physical_offset[i] = domain_offset;
            physical_size[i] = domain_max - domain_min;
            cell_size[i] = (domain_max - domain_min) / grid_shape[i];
        }
        for (unsigned int i = 0; i < expected_2dc.shape()[0]; ++i) {
            for (unsigned int j = 0; j < expected_2dc.shape()[1]; ++j) {
                expected_2dc[i][j] = 0.0;
            }
        }
        for (unsigned int k = 0; k < expected_1d.shape()[0]; ++k) {
            expected_1d[k] = 0.0;
        }
        density_norm_2d = (toy_real_num / total_num) * pconstants::e
                / (cell_size[0] * cell_size[1]);
        density_norm_1d = 1.0;
    }

    ~Toy_bunch_fixture_2d()
    {
    }

    Four_momentum four_momentum;
    Reference_particle reference_particle;
    Commxx_sptr comm_sptr;
    Bunch bunch;
    double density_norm_2d, density_norm_1d;

    std::vector<double > physical_size, physical_offset, cell_size;
    std::vector<int > grid_shape;
    Rectangular_grid_sptr rho_grid_sptr;
    MArray2dc expected_2dc;
    MArray1d expected_1d;
};

#endif /* SPACE_CHARGE_BUNCH_FIXTURES_H_ */
