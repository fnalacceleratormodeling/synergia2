#include <iostream>
#include <vector>
#include <iterator>
#include "array_2d.h"
#include "vector_helper.h"

void
output_vector(std::vector<int> v)
{
    std::cout << "[";
    std::copy(v.begin(),v.end(),std::ostream_iterator<int>(std::cout,","));
    std::cout << "]";
}

int
main()
{
    Array_2d<double> x1;
    x1.describe();
    x1.reshape(3,2);
    x1.describe();
    x1.set_all(0.0);
    x1.print("x1");
    std::cout << "adding 3.14 to x1\n";
    x1.add(3.14);
    x1.print("new x1");
    
    std::cout << "making x2 a copy of x1\n";
    Array_2d<double> x2 = x1;
    x2.print("x2");

    std::cout << "making x3 point to data from x1\n";
    Array_2d<double> x3(3,2,x1.get_data_ptr());
    x3.print("x3");
    
    std::cout << "scaling x1 by 20\n";
    x1.scale(20);
    x1.print("new x1");
    x2.print("new x2");
    x3.print("new x3");
    
    Array_2d<double> x5;
    double *mydata;
    std::vector<int> shape2(vector2(1,2));
    mydata = (double*)malloc(shape2[0]*shape2[1]*sizeof(double));
    x5.reshape(shape2[0],shape2[1],mydata);
    std::vector<int> index(2);
    for (index[0]=0; index[0]<shape2[0]; ++index[0]) {
        for (index[1]= 0; index[1]<shape2[1]; ++index[1]) {
            x5.at(index[0],index[1]) = 2.0*(index[0]+ index[1]*0.1);
        }
    }
    x5.print("x5");
    
    x5.freeze_shape();
    std::cout << "attempting to reshape frozen array... ";
    try {
        x5.reshape(3,4);
        std::cout << "no error!\n";
    } catch(std::out_of_range) {
        std::cout << "caught out_of_range error!\n";
    }

    std::cout << "attempting to set x5.at(-1,0)... ";
    try {
        x5.at(-1,0) = -999;
        std::cout << "no error!\n";
    } catch(std::out_of_range) {
        std::cout << "caught out_of_range error!\n";
    }
    
    std::cout << "x1.owns_data() = " << x1.owns_data() << std::endl;
    std::cout << "x2.owns_data() = " << x2.owns_data() << std::endl;
    std::cout << "x3.owns_data() = " << x3.owns_data() << std::endl;
    std::cout << "x5.owns_data() = " << x5.owns_data() << std::endl;
    
    std::vector<int> retrieved_shape(x5.get_shape());
    std::cout << "x5.get_shape() = ";
    output_vector(retrieved_shape);
    std::cout << std::endl;

    std::vector<int> retrieved_strides(x5.get_strides());
    std::cout << "x5.get_strides() = ";
    output_vector(retrieved_strides);
    std::cout << std::endl;
 
    std::cout << "x5.offset(1,1) = " << x5.offset(1,1) << std::endl; 
    
    std::cout << "x1.get_size() = " << x1.get_size() << std::endl;
    std::cout << "x2.get_size() = " << x2.get_size() << std::endl;
    std::cout << "x3.get_size() = " << x3.get_size() << std::endl;
    std::cout << "x5.get_size() = " << x5.get_size() << std::endl;
 
    std::cout << "x1.shape_is_frozen() = " << x1.shape_is_frozen() << std::endl;
    std::cout << "x2.shape_is_frozen() = " << x2.shape_is_frozen() << std::endl;
    std::cout << "x3.shape_is_frozen() = " << x3.shape_is_frozen() << std::endl;
    std::cout << "x5.shape_is_frozen() = " << x5.shape_is_frozen() << std::endl;
 

    return 0;
}
