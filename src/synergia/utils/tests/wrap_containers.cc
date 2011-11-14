#include "wrap_containers.h"

void
take_vector_int(std::vector<int >  const & int_list)
{
    std::cout<<" list accepted as vector int in C "<<std::endl;
    for (std::vector<int>::const_iterator it=int_list.begin(); it !=int_list.end(); ++it){
        std::cout<<" it="<<*it<<std::endl;    
    }    
           
}
    
void
take_list_int(std::list<int > const & int_list) 
{
    std::cout<<" list accepted as list int in C "<<std::endl;
    for (std::list<int>::const_iterator it=int_list.begin(); it !=int_list.end(); ++it){
        std::cout<<" it="<<*it<<std::endl;    
    }    
           
}   
    
void
take_vector_dd(std::vector<double >  const & dd_list)
{
    std::cout<<" list accepted as vector double in C "<<std::endl;
    for (std::vector<double>::const_iterator it=dd_list.begin(); it !=dd_list.end(); ++it){
        std::cout<<" it="<<*it<<std::endl;    
    }     
   
}

    