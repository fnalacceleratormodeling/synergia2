#include <iostream>
#include <math.h>
#include <nlopt.hpp>

int count = 0;

double myfunc(const std::vector<double > &x, std::vector<double > &grad, void *my_func_data)
{
    ++count;
    if (!grad.empty()) {
        grad[0] = 0.0;
        grad[1] = 0.5 / sqrt(x[1]);
    }
    std::cout << "try " << count << ": x[0] = " << x[0] << ", x[1] = " << x[1]
            << ", f(x[0], x[1]) = " << sqrt(x[1]) << std::endl;
    return sqrt(x[1]);
}

typedef struct {
    double a, b;
} my_constraint_data;

double myconstraint(const std::vector<double> &x, std::vector<double> &grad, void *data)
{
    my_constraint_data *d = reinterpret_cast<my_constraint_data*>(data);
    double a = d->a, b = d->b;
    if (!grad.empty()) {
        grad[0] = 3 * a * (a*x[0] + b) * (a*x[0] + b);
        grad[1] = -1.0;
    }
    return ((a*x[0] + b) * (a*x[0] + b) * (a*x[0] + b) - x[1]);
}

int main() {
    //nlopt::opt opt(nlopt::LD_MMA, 2);
    nlopt::opt opt(nlopt::LN_COBYLA, 2);

    std::vector<double> lb(2);
    lb[0] = -HUGE_VAL; lb[1] = 0;
    opt.set_lower_bounds(lb);

    opt.set_min_objective(myfunc, NULL);

    my_constraint_data data[2] = { {2,0}, {-1,1} };
    opt.add_inequality_constraint(myconstraint, &data[0], 1e-8);
    opt.add_inequality_constraint(myconstraint, &data[1], 1e-8);

    opt.set_xtol_rel(1e-4);
    std::vector<double> x(2);
    x[0] = 1.234; x[1] = 5.678;  /* some initial guess */
    double minf;                 /* the minimum objective value, upon return */
    nlopt::result result = opt.optimize(x, minf);

    std::cout << "found minimum at f(" << x[0] << ", " << x[1] << ") = "
            << minf << std::endl;
    std::cout << "found minimum after " << count << " evaluation" << std::endl;

    return 0;
}
