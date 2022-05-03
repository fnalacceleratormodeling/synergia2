#ifndef DISTRIBUTED_RECTANGULAR_GRID_H_
#define DISTRIBUTED_RECTANGULAR_GRID_H_
#include "synergia/collective/rectangular_grid_domain.h"
#include "synergia/utils/multi_array_typedefs.h"
#include "synergia/utils/commxx.h"

class Distributed_rectangular_grid
{
private:
    Rectangular_grid_domain_sptr domain_sptr;
    boost::shared_ptr<Raw_MArray3d > grid_points_sptr;
    boost::shared_ptr<Raw_MArray2dc > grid_points_2dc_sptr;
    boost::shared_ptr<Raw_MArray1d > grid_points_1d_sptr;
    int lower, upper;
    int lower_guard, upper_guard;
    std::vector<int> uppers, lengths;
    double normalization;
    Commxx_sptr comm_sptr;
    void
    construct_hockney(int lower, int upper, std::vector<int > const & array_shape);
    void
    construct_rectangular(int lower, int upper, std::vector<int > const & array_shape);
    void
    calculate_uppers_lengths();
public:
    Distributed_rectangular_grid(std::vector<double > const & physical_size,
            std::vector<double > const & physical_offset,
            std::vector<int > const & grid_shape, bool periodic, int lower,
            int upper, Commxx_sptr comm_sptr, std::string const solver="hockney");
    Distributed_rectangular_grid(
            Rectangular_grid_domain_sptr rectangular_grid_domain_sptr,
            int lower, int upper, Commxx_sptr comm_sptr, std::string const solver="hockney");
    Distributed_rectangular_grid(
            Rectangular_grid_domain_sptr rectangular_grid_domain_sptr,
            int lower, int upper, std::vector<int > const & padded_shape,
            Commxx_sptr comm_sptr);
    Rectangular_grid_domain const&
    get_domain() const
    {
        return *domain_sptr;
    }
    Rectangular_grid_domain &
    get_domain()
    {
        return *domain_sptr;
    }
    Rectangular_grid_domain_sptr
    get_domain_sptr() const;
    Rectangular_grid_domain_sptr
    get_domain_sptr();
    int
    get_lower() const;
    int
    get_upper() const;
    int
    get_lower_guard() const;
    int
    get_upper_guard() const;
    std::vector<int > const&
    get_uppers();
    std::vector<int > const&
    get_lengths();
    MArray3d_ref const&
    get_grid_points() const;
    MArray3d_ref &
    get_grid_points();
    MArray2dc_ref const&
    get_grid_points_2dc() const;
    MArray2dc_ref &
    get_grid_points_2dc();
    MArray1d_ref const&
    get_grid_points_1d() const;
    MArray1d_ref &
    get_grid_points_1d();
    void
    set_normalization(double val);
    double
    get_normalization() const;
    Commxx const &
    get_comm() const;
    Commxx_sptr
    get_comm_sptr() const;
    void
    fill_guards();
};

typedef boost::shared_ptr<Distributed_rectangular_grid >
        Distributed_rectangular_grid_sptr;  // syndoc:include
#endif /* DISTRIBUTED_RECTANGULAR_GRID_H_ */
