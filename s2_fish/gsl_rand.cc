#include <iostream>

#include <boost/python.hpp>
#include <boost/python/numeric.hpp>
#include <numpy/arrayobject.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace boost::python;

numeric::array create_double_array(tuple &dims)
{
    int cdims[2];
    cdims[0] = extract<int>(dims[0]);
    cdims[1] = extract<int>(dims[1]);
    object obj(handle<>(PyArray_FromDims(2, cdims, PyArray_DOUBLE)));
    return extract<numeric::array>(obj);
}

class GSL_rand
{
private:
    gsl_rng * rng;
public:
    GSL_rand();
    numeric::array get_array(tuple &dims);
    ~GSL_rand();
};

GSL_rand::GSL_rand()
{
    gsl_rng_env_setup();
    rng = gsl_rng_alloc(gsl_rng_ranlxd2);
    gsl_rng_set(rng, 0);
    import_array();
}

numeric::array
GSL_rand::get_array(tuple &dims)
{
    numeric::array retval = create_double_array(dims);
    PyArrayObject* array_object = reinterpret_cast<PyArrayObject*>
                                  (retval.ptr());
    double *data = reinterpret_cast<double*>(array_object->data);
    int dim0 = extract<int>(dims[0]);
    int dim1 = extract<int>(dims[1]);
    double sigma = 2.0;
    //~ for(int i = 0; i<dim0; ++i) {
    //~ for(int j = 0; j<dim1; ++j) {
    //~ data[i*dim1+j] = gsl_ran_gaussian(rng,sigma);
    //~ }
    //~ }
    for (int i = 0; i < dim0*dim1; ++i) {
        data[i] = gsl_ran_ugaussian_ratio_method(rng);
    }
    return retval;
}

GSL_rand::~GSL_rand()
{
    gsl_rng_free(rng);
}

BOOST_PYTHON_MODULE(gsl_rand)
{
    numeric::array::set_module_and_type("numpy", "ndarray");
    class_<GSL_rand>("GSL_rand",
                     init<>())
    .def("get_array",
         &GSL_rand::get_array)
    ;
}

