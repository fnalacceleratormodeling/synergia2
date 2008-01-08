#include <iostream>
#include <vector>
#include <iterator>
#include "array_nd.h"
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
    std::vector<int> shape1(vector2(3,2));    
    std::cout << "shape 1 is ";
    output_vector(shape1);
    std::cout << std::endl;
    
    Array_nd<double> x1;
    x1.describe();
    x1.reshape(shape1);
    x1.describe();
    x1.set_all(0.0);
    x1.print("x1");
    std::cout << "adding 3.14 to x1\n";
    x1.add(3.14);
    x1.print("new x1");
    
    std::cout << "making x2 a copy of x1\n";
    Array_nd<double> x2 = x1;
    x2.print("x2");

    std::cout << "making x3 point to data from x1\n";
    Array_nd<double> x3(shape1,x1.get_data_ptr());
    x3.print("x3");
    
    std::cout << "scaling x1 by 20\n";
    x1.scale(20);
    x1.print("new x1");
    x2.print("new x2");
    x3.print("new x3");
    
    std::vector<int> shape2(vector3(2,3,4)),index(3);
    Array_nd<double> x4(shape2);
    for (index[0]=0; index[0]<shape2[0]; ++index[0]) {
        for (index[1]= 0; index[1]<shape2[1]; ++index[1]) {
            for (index[2]= 0; index[2]<shape2[2]; ++index[2]) {
                x4(index) = index[0]+ index[1]*0.1 + index[2]*0.01;
            }
        }
    }
    x4.print("x4");
    std::cout << "(not) resizing x4\n";
    x4.reshape(shape2);
    x4.print("x4");
    
    Array_nd<double> x5;
    double *mydata;
    mydata = (double*)malloc(shape2[0]*shape2[1]*shape2[2]*sizeof(double));
    x5.reshape(shape2,mydata);
    for (index[0]=0; index[0]<shape2[0]; ++index[0]) {
        for (index[1]= 0; index[1]<shape2[1]; ++index[1]) {
            for (index[2]= 0; index[2]<shape2[2]; ++index[2]) {
                x5.at(index) = 2.0*(index[0]+ index[1]*0.1 + index[2]*0.01);
            }
        }
    }
    x5.print("x5");
    
    x5.freeze_shape();
    std::cout << "attempting to reshape frozen array... ";
    try {
        x5.reshape(shape1);
        std::cout << "no error!\n";
    } catch(std::out_of_range) {
        std::cout << "caught out_of_range error!\n";
    }

    std::cout << "attempting to set x5.at(-1,0,0)... ";
    try {
        index = vector3(-1,0,0);
        x5.at(index) = -999;
        std::cout << "no error!\n";
    } catch(std::out_of_range) {
        std::cout << "caught out_of_range error!\n";
    }
    
    std::cout << "x5.bounds_check(-1,0,0) = "
        << x5.bounds_check(index) << std::endl;
    
    std::cout << "attempting to set x5.at(1,2,4)... ";
    try {
        index = vector3(1,2,4);
        x5.at(index) = -999;
        std::cout << "no error!\n";
    } catch(std::out_of_range) {
        std::cout << "caught out_of_range error!\n";
    }
    
    double dummy;
    std::cout << "attempting to get x5.at(0,-1,0)... ";
    try {
        index = vector3(0,-1,0);
        dummy = x5.at(index);
        std::cout << "no error!\n";
    } catch(std::out_of_range) {
        std::cout << "caught out_of_range error!\n";
    }
    
    std::cout << "attempting to get x5.at(1,99,3)... ";
    try {
        index = vector3(1,99,3);
        dummy = x5.at(index);
        std::cout << "no error!\n";
    } catch(std::out_of_range) {
        std::cout << "caught out_of_range error!\n";
    }
    
    std::cout << "x1.owns_data() = " << x1.owns_data() << std::endl;
    std::cout << "x2.owns_data() = " << x2.owns_data() << std::endl;
    std::cout << "x3.owns_data() = " << x3.owns_data() << std::endl;
    std::cout << "x4.owns_data() = " << x4.owns_data() << std::endl;
    std::cout << "x5.owns_data() = " << x5.owns_data() << std::endl;
    
    std::vector<int> retrieved_shape(x5.get_shape());
    std::cout << "x5.get_shape() = ";
    output_vector(retrieved_shape);
    std::cout << std::endl;
 
    index = vector3(1,2,3);
    std::cout << "x5.offset(1,2,3) = " << x5.offset(index) << std::endl; 
    
    std::cout << "x1.get_size() = " << x1.get_size() << std::endl;
    std::cout << "x2.get_size() = " << x2.get_size() << std::endl;
    std::cout << "x3.get_size() = " << x3.get_size() << std::endl;
    std::cout << "x4.get_size() = " << x4.get_size() << std::endl;
    std::cout << "x5.get_size() = " << x5.get_size() << std::endl;
 
    std::cout << "x1.shape_is_frozen() = " << x1.shape_is_frozen() << std::endl;
    std::cout << "x2.shape_is_frozen() = " << x2.shape_is_frozen() << std::endl;
    std::cout << "x3.shape_is_frozen() = " << x3.shape_is_frozen() << std::endl;
    std::cout << "x4.shape_is_frozen() = " << x4.shape_is_frozen() << std::endl;
    std::cout << "x5.shape_is_frozen() = " << x5.shape_is_frozen() << std::endl;
 

    return 0;
}
