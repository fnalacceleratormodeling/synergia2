#ifndef RECTANGULAR_GRID_DOMAIN_FIXTURE_H_
#define RECTANGULAR_GRID_DOMAIN_FIXTURE_H_

const double domain_min = -1.0;
const double domain_max = 1.0;
const double domain_offset = 5.0;
const int grid_size0 = 4;
const int grid_size1 = 5;
const int grid_size2 = 3;
const bool is_periodic = false;

struct Rectangular_grid_domain_fixture
{
    Rectangular_grid_domain_fixture() :
        physical_size(3), physical_offset(3), grid_shape(3)
    {
        for (int i = 0; i < 3; ++i) {
            physical_offset[i] = domain_offset;
            physical_size[i] = domain_max - domain_min;
        }
        grid_shape[0] = grid_size0;
        grid_shape[1] = grid_size1;
        grid_shape[2] = grid_size2;
        rectangular_grid_domain_sptr = Rectangular_grid_domain_sptr(
                new Rectangular_grid_domain(physical_size, physical_offset,
                        grid_shape, is_periodic));
    }

    ~Rectangular_grid_domain_fixture()
    {
    }

    std::vector<double > physical_size, physical_offset;
    std::vector<int > grid_shape;
    Rectangular_grid_domain_sptr rectangular_grid_domain_sptr;
};

#endif /* RECTANGULAR_GRID_DOMAIN_FIXTURE_H_ */
