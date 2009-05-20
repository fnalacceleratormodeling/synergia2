#ifndef HAVE_FIELD_DOMAIN_H
#define HAVE_FIELD_DOMAIN_H

#include <vector>

class Field_domain
{
private:
    std::vector<double> physical_size;
    std::vector<double> physical_offset;
    std::vector<int> grid_shape;
    std::vector<bool> periodic;

    std::vector<double> left;
    std::vector<double> cell_size;

    void construct(const std::vector<double> &physical_size,
        const std::vector<double> &physical_offset,
        const std::vector<int> &grid_shape,
        const std::vector<bool> &periodic);

public:
    Field_domain();
    Field_domain(const std::vector<double> &physical_size,
        const std::vector<double> &physical_offset,
        const std::vector<int> &grid_shape,
        const std::vector<bool> &periodic);

    void set_params(const std::vector<double> &physical_size,
        const std::vector<double> &physical_offset,
        const std::vector<int> &grid_shape,
        const std::vector<bool> &periodic);

    std::vector<int> get_grid_shape() const;
    std::vector<double> get_physical_size() const;
    std::vector<double> get_cell_size() const;
    std::vector<bool> get_periodic() const;
    void get_leftmost_indices_offsets(double c1, double c2, double c3,
        std::vector<int> &indices, std::vector<double> &offsets) const;
};

class Cylindrical_field_domain
{
    private:
        double radius;
        double length,half_length;
        std::vector<int> grid_shape;
        bool periodic_z;

        std::vector<double> cell_size;

    public:
        Cylindrical_field_domain(double radius, double length,
                                 const std::vector<int> &grid_shape,
                                 bool periodic_z);

        // jfa: The name get_leftmost_indices_offsets is possibly misleading
        void get_leftmost_indices_offsets(double c0, double c1, double c2,
                                          std::vector<int> &indices,
                                          std::vector<double> &offsets) const;
        const std::vector<int> &get_grid_shape() const;
        const std::vector<double> &get_cell_size() const;
        const double get_length() const;
        const double get_radius() const;
        const bool get_periodic_z() const;
};

#endif // HAVE_FIELD_DOMAIN_H
