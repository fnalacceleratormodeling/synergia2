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
                four_momentum), comm(MPI_COMM_WORLD), bunch(reference_particle,
                total_num, real_num, comm), distribution(0, comm),
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
        covariances[1][1] = covariances[3][3] = covariances[5][5] = 1.0;
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
    Commxx comm;
    Bunch bunch;
    Random_distribution distribution;
    double stdx, stdy, stdz;
    std::vector<int > grid_shape;
};

struct Spherical_bunch_fixture
{
    Spherical_bunch_fixture() :
        four_momentum(mass, total_energy), reference_particle(charge,
                four_momentum), comm(MPI_COMM_WORLD), bunch(reference_particle,
                total_num, real_num, comm), distribution(0, comm),
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
        covariances[1][1] = covariances[3][3] = covariances[5][5] = 1.0;
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
    Commxx comm;
    Bunch bunch;
    Random_distribution distribution;
    double sigma;
    std::vector<int > grid_shape;
};

struct Cylindrical_bunch_fixture
{
    Cylindrical_bunch_fixture() :
        four_momentum(mass, total_energy), reference_particle(charge,
                four_momentum), comm(MPI_COMM_WORLD), bunch(reference_particle,
                total_num, real_num, comm), distribution(0, comm),
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
        covariances[1][1] = covariances[3][3] = covariances[5][5] = 1.0;
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
    Commxx comm;
    Bunch bunch;
    Random_distribution distribution;
    double sigma;
    std::vector<int > grid_shape;
};

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
                comm(MPI_COMM_WORLD),
                bunch(reference_particle, total_num, toy_real_num, comm),
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
    Commxx comm;
    Bunch bunch;
    double density_norm;

    std::vector<double > physical_size, physical_offset, cell_size;
    std::vector<int > grid_shape;
    Rectangular_grid_sptr rho_grid_sptr;
    MArray3d expected;
};

#endif /* SPACE_CHARGE_BUNCH_FIXTURES_H_ */
