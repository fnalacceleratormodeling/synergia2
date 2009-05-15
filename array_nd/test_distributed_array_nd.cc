#include <iostream>
#include <vector>
#include <iterator>
#include "distributed_array_nd.h"
#include "vector_helper.h"


int
main()
{
	Distributed_array_nd<double> x1(vector2(0,0),vector2(8,8),
    		vector2(0,0),vector2(4,4),vector2(false,true),1);
    Distributed_array_nd<double> x2(vector2(0,0),vector2(8,8),
    		vector2(4,4),vector2(8,8),vector2(false,true),1);
    x1.set_all(1.0);
    x2.set_all(2.0);
    x1.at(vector2(0,0)) = 123.0;
    x2.at(vector2(4,4)) = 456.0;
    x1.get_local_array_wguards().print("x1 with guards");
    x2.get_local_array_wguards().print("x2 with guards");
    x1.get_local_array().print("x1 no guards");
    x2.get_local_array().print("x2 no guards");
    return 0;
}
