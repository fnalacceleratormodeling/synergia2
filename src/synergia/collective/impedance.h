#ifndef IMPEDANCE_H_
#define  IMPEDANCE_H_
#include "synergia/simulation/operator.h"
#include "synergia/simulation/step.h"
#include "synergia/bunch/bunch.h"

class Impedance : public Collective_operator
{
private:
    std::string pipe_symmetry;
    int z_grid;
    int nstored_turns;
    double orbit_length;
    double wake_factor;
    double bunch_spacing;
    std::string wake_file;
    std::vector<double> z_coord;
    std::vector<double>  x_wake;
    std::vector<double> y_wake;
    std::vector<double> z_wake;

public:

    Impedance(std::string const & wake_file, double const & orbit_length, double const & bunchsp, 
                int const  & zgrid, std::string const & pipe_symmetry, int const nstored_turns);


    int get_z_grid() const;
    double get_orbit_length() const;
    double get_wake_factor() const;
    double get_bunch_spacing() const;
    std::string get_pipe_symmetry() const;
    std::string get_wake_file_name() const;
    std::vector<double> get_z_coord() const;
    std::vector<double> get_x_wake() const;
    std::vector<double> get_y_wake() const;
    std::vector<double> get_z_wake() const;
    virtual int 
    get_nstored_turns() const;
    virtual void
    apply(Bunch & bunch, double time_step, Step & step);
    virtual
    ~Impedance();
};

typedef boost::shared_ptr<Impedance> Impedance_sptr;

#endif /* IMPEDANCE_H_ */

