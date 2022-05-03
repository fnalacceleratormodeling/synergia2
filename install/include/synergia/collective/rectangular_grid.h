#ifndef RECTANGULAR_GRID_H_
#define RECTANGULAR_GRID_H_

#include "synergia/utils/multi_array_typedefs.h"


class Rectangular_grid_1d
{

private:

    using kokkos_noinit = Kokkos::ViewAllocateWithoutInitializing;

    karray1d_row_dev grid_;
    karray1d_row_hst hgrid_;

    double normalization;

public:

    Rectangular_grid_1d(size_t shape_x, bool zero = true)
    : grid_( zero ? karray1d_row_dev("g1d", shape_x) 
                  : karray1d_row_dev(kokkos_noinit("g1d"), shape_x) )
    , hgrid_(Kokkos::create_mirror_view(grid_))
    , normalization(1.0)
    { }

    size_t shape(size_t dim) const
    { return grid_.extent(dim); }

    size_t span() const
    { return grid_.span(); }
    
    double* data(size_t x = 0) const
    { return &grid_(x); }

    //void set_zero()
    //{ points_.setZero(); }

    void set_normalization(double val)
    { normalization = val; }

    double get_normalization() const
    { return normalization; }
};


class Rectangular_grid_2dc
{

private:

    using kokkos_noinit = Kokkos::ViewAllocateWithoutInitializing;

    karray2dc_row_dev grid_;
    karray2dc_row_hst hgrid_;

    double normalization;

public:

    Rectangular_grid_2dc(size_t shape_x, size_t shape_y, bool zero = true)
    : grid_( zero ? karray2dc_row_dev("g2dc", shape_x, shape_y) 
                  : karray2dc_row_dev(kokkos_noinit("g2dc"), shape_x, shape_y) )
    , hgrid_(Kokkos::create_mirror_view(grid_))
    , normalization(1.0)
    { }

    size_t shape(size_t dim) const
    { return grid_.extent(dim); }

    size_t span() const
    { return grid_.span(); }
    
    Kokkos::complex<double>* data(size_t x = 0, size_t y = 0) const
    { return &grid_(x, y); }

    //void set_zero()
    //{ points_.setZero(); }

    void set_normalization(double val)
    { normalization = val; }

    double get_normalization() const
    { return normalization; }
};


#if 0
template<typename T = double>
class Rectangular_grid_eigen
{
public:

    typedef Eigen::Tensor<T, 3, Eigen::RowMajor> EArray3d;

private:

    EArray3d  points_;
    std::array<int, 3> shape_;

    double normalization;

public:

    Rectangular_grid_eigen(std::array<int, 3> const & grid_shape, bool zero = true)
        : points_(grid_shape[0], grid_shape[1], grid_shape[2])
        , shape_(grid_shape)
        , normalization(1.0)
    { if (zero) set_zero(); }

    EArray3d const & get_grid_points() const
    { return points_; }

    EArray3d & get_grid_points()
    { return points_; }

    std::array<int, 3> const & shape() const
    { return shape_; }

    typename EArray3d::Scalar const &
    grid(Eigen::Index x, Eigen::Index y, Eigen::Index z) const
    { return points_(x, y, z); }

    typename EArray3d::Scalar &
    grid(Eigen::Index x, Eigen::Index y, Eigen::Index z)
    { return points_(x, y, z); }

    typename EArray3d::Scalar const *
    data(Eigen::Index x = 0, Eigen::Index y = 0, Eigen::Index z = 0) const
    { return points_.data() + x * shape_[1] * shape_[2] + y * shape_[2] + z; }

    void set_zero()
    { points_.setZero(); }

    void set_normalization(double val)
    { normalization = val; }

    double get_normalization() const
    { return normalization; }
};
#endif

#endif /* RECTANGULAR_GRID_EIGEN_H_ */
