#include <iostream>
#include "array_nd/array_nd.h"
#include "array_nd/array_1d.h"

int main()
{
    int imax=7;
    int jmax=5;
    Array_nd<double> a2n(vector2(imax,jmax));
    a2n.set_all(0.0);
    for(int i=0;i<imax;++i) {
        for(int j=0; j<jmax; ++j) {
            a2n(vector2(i,j)) = i*10+j;
        }
    }
    a2n.print("a2n");
    
    Array_nd<double> s1 = a2n.slice(vector2(Range(0,3),Range(1)));
    s1.print("s1");
    
    Array_nd<double> s1t = a2n.slice(vector2(Range(1),Range(0,3)));
    s1t.print("s1t");
    
    Array_nd<double> s2= a2n.slice(vector2(Range(0),Range(0,Range::end)));
    s2.print("s2");
    
    Array_1d<double> s22= a2n.slice(vector2(Range(0),Range(0,Range::end,2)));
    s22.describe();
    s22.print("s22");
    
    //~ Array_nd<double>::Iterator itx = s22.begin();
    //~ double * the_end = s22.end();
    //~ std::cout << "itx = " << *itx << ", end = " << the_end << std::endl;
    
    std::cout << "s22 iteration:\n";
    for(Array_nd<double>::Iterator it = s22.begin(); 
            it != s22.end();
            ++it) {
        std::cout << *it << std::endl;
    }
    
    std::cout << "s22 iteration:\n";
    for(Array_nd<double>::Iterator it = s22.begin(); 
            it != s22.end();
            ++it) {
        *it = -99.9;
    }
    s22.print("s22");
    
    //~ std::cout << "\na2n iteration:\n";
    //~ for(Array_nd<double>::Iterator it = a2n.begin(); 
            //~ *it != a2n.end();
            //~ ++it) {
        //~ std::cout << **it << std::endl;
    //~ }
    
    return 0;
}
