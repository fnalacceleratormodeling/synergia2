#include "split_operator_stepper_choice.h"


//Split_operator_stepper_choice

template<class Archive>
    void
    Kicks::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP(collective_operators);
        ar & BOOST_SERIALIZATION_NVP(num_steps);
    }

template
void
Kicks::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Kicks::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Kicks::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Kicks::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);
BOOST_CLASS_EXPORT_IMPLEMENT(Kicks);



void
Split_operator_stepper_choice::construct_per_element_else()
{
    bool verbose=false;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    double t_turn, t_turn1;
    t_turn= MPI_Wtime();
    // check if all num_steps are >1
    for(List_choice_map::const_iterator lch_it=this->list_choice_map.begin(); lch_it!= this->list_choice_map.end();
                 ++lch_it){
            if (lch_it->second.num_steps < 1) {
                 std::cout<<"element name="<<  lch_it->first<<" num steps= "<<lch_it->second.num_steps<<std::endl;
                throw std::runtime_error(
                    "Split_operator_stepper_choice: all num_steps must be >= 1");
        }
    }

    if (this->num_steps_else<1) throw std::runtime_error(
                    "Split_operator_stepper_choice: num_steps_else must be >= 1");
    // find the length of the lattice without the chosen elements
    double length_else=0.;
    for (Lattice_elements::iterator it =
            get_lattice_simulator().get_lattice_sptr()->get_elements().begin(); it
            != get_lattice_simulator().get_lattice_sptr()->get_elements().end(); ++it){

        List_choice_map::const_iterator lch_it= this->list_choice_map.find((*it)->get_name());
        if (lch_it == this->list_choice_map.end()) length_else +=(*it)->get_length();

    }

    double step_length_else = length_else/this->num_steps_else;

    for (Lattice_elements::iterator it =
            get_lattice_simulator().get_lattice_sptr()->get_elements().begin(); it
            != get_lattice_simulator().get_lattice_sptr()->get_elements().end(); ++it){
        double length = (*it)->get_length();
        if (length == 0.0) {
            Independent_operator_sptr
                    ind_op(
                            new Independent_operator(
                                    "step",
                                    get_lattice_simulator().get_operation_extractor_map_sptr(),
                                    get_lattice_simulator().get_aperture_operation_extractor_map_sptr()));
            Lattice_element_slice_sptr slice(new Lattice_element_slice(*it));
            ind_op->append_slice(slice);
            Step_sptr step(new Step(0.0));
            step->append(ind_op, 1.0);
            get_steps().push_back(step);
            if ((rank==0) && verbose) std::cout<<" element: "<<(*it)->get_name()<<" zero length step"<<std::endl;
         }
         else{
            List_choice_map::const_iterator lch_it= this->list_choice_map.find((*it)->get_name());
            if  (lch_it != this->list_choice_map.end()){
               double step_length= length/lch_it->second.num_steps;
               double half_step_length = 0.5 * step_length;
               int steps_per_element=lch_it->second.num_steps;
                for (int i = 0; i <steps_per_element;  ++i) {
                    double left = i * step_length;
                    double middle = left + half_step_length;
                    double right = (i + 1) * step_length;

                    Step_sptr step(new Step(step_length));
                    //1st Half
                    Independent_operator_sptr
                        ind_op_first_half(
                                new Independent_operator(
                                        "first_half_chosen",
                                        get_lattice_simulator().get_operation_extractor_map_sptr(),
                                        get_lattice_simulator().get_aperture_operation_extractor_map_sptr()));
                    Lattice_element_slice_sptr slice_1st_half(
                            new Lattice_element_slice(*it, left, middle));
                    ind_op_first_half->append_slice(slice_1st_half);
                    step->append(ind_op_first_half, 0.5);


                     //Collective Effects
                    for (Collective_operators::const_iterator coll_op_it =
                        lch_it->second.collective_operators.begin(); coll_op_it
                        != lch_it->second.collective_operators.end();
                        ++coll_op_it) {
                        Collective_operator_sptr copied_collective_operator_sptr(
                                                                                    (*coll_op_it)->clone());
                    step->append(copied_collective_operator_sptr, 1.0);
                    }

                        //2nd Half
                    Independent_operator_sptr
                        ind_op_second_half(
                                new Independent_operator(
                                        "second_half_chosen",
                                        get_lattice_simulator().get_operation_extractor_map_sptr(),
                                        get_lattice_simulator().get_aperture_operation_extractor_map_sptr()));
                    Lattice_element_slice_sptr slice_2nd_half(
                        new Lattice_element_slice(*it, middle, right));
                    ind_op_second_half->append_slice(slice_2nd_half);
                    step->append(ind_op_second_half, 0.5);

                    get_steps().push_back(step);

                }
            if ((rank==0) && verbose) std::cout<<" element: "<<(*it)->get_name()<<" nr of collective kicks ="<<steps_per_element<<std::endl;
            }
            else{
                        int steps_per_element=int(ceil(length/step_length_else));
                        double step_length=length/steps_per_element;
                        double half_step_length=0.5*step_length;

                        for(int i = 0; i < steps_per_element; ++i){
                            double left = i * step_length;
                            double middle = left + half_step_length;
                            double right = (i + 1) * step_length;

                            Step_sptr step(new Step(step_length));
                          //1st Half
                            Independent_operator_sptr
                                ind_op_first_half(
                                    new Independent_operator(
                                            "first_half_else",
                                                get_lattice_simulator().get_operation_extractor_map_sptr(),
                                                    get_lattice_simulator().get_aperture_operation_extractor_map_sptr()));
                            Lattice_element_slice_sptr slice_1st_half(
                                    new Lattice_element_slice(*it, left, middle));
                            ind_op_first_half->append_slice(slice_1st_half);
                            step->append(ind_op_first_half, 0.5);

                            //Collective Effects
                            for (Collective_operators::const_iterator coll_op_it =
                                this->list_choice_map["else"].collective_operators.begin(); coll_op_it
                                != this->list_choice_map["else"].collective_operators.end(); ++coll_op_it) {
                                Collective_operator_sptr copied_collective_operator_sptr(
                                                                                            (*coll_op_it)->clone());
                                step->append(copied_collective_operator_sptr, 1.0);
                            }
                              //2nd Half
                            Independent_operator_sptr
                                ind_op_second_half(
                                    new Independent_operator(
                                            "second_half_else",
                                            get_lattice_simulator().get_operation_extractor_map_sptr(),
                                            get_lattice_simulator().get_aperture_operation_extractor_map_sptr()));
                            Lattice_element_slice_sptr slice_2nd_half(
                                    new Lattice_element_slice(*it, middle, right));

                            ind_op_second_half->append_slice(slice_2nd_half);
                            step->append(ind_op_second_half, 0.5);

                            get_steps().push_back(step);

                   }
            if ((rank==0) && verbose) std::cout<<" element: "<<(*it)->get_name()<<" nr of \"else\" collective kicks ="<<steps_per_element<<std::endl;
            }
         }

    }

    get_lattice_simulator().set_slices(extract_slices(get_steps()));

    t_turn1= MPI_Wtime();
     if (rank==0){
          std::cout<<" stepper_choice per element:"<<std::endl;
          std::cout<<"                total number of steps=  "<<this->get_steps().size()<<std::endl;
          std::cout<<" time: stepper choice done in "<<t_turn1-t_turn<<std::endl;
     }

}





void
Split_operator_stepper_choice::
make_stepper_else(Lattice_elements::iterator const& begin, Lattice_elements::iterator const& end,
    double const & length_between, double const & max_step_length)
{
    int num_steps=int(ceil(length_between/max_step_length));
    double step_length=length_between/num_steps;
    double  half_step_length=0.5*step_length;
    double left = 0.0;
    double offset_fudge = 0.0;
    Lattice_elements::iterator lattice_it =begin;
   // Lattice_elements::iterator lattice_end =end;
    for (int i = 0; i < num_steps; ++i) {
        Step_sptr step(new Step(step_length));
        step->append(
                Stepper::get_fixed_step("first_half_else", lattice_it, left, end,
                       half_step_length, offset_fudge, false), 0.5);



        for (Collective_operators::const_iterator coll_op_it =
                this->list_choice_map["else"].collective_operators.begin(); coll_op_it
                != this->list_choice_map["else"].collective_operators.end(); ++coll_op_it) {
            step->append(*coll_op_it, 1.0);
        }
        step->append(
                Stepper::get_fixed_step("second_half_else", lattice_it, left, end,
                        half_step_length, offset_fudge, false), 0.5);
        get_steps().push_back(step);
    }
    if (lattice_it != end) {
        throw(std::runtime_error(
                "internal error: split_operator_stepper_choice: make_stepper_else: did not make it to the end of the lattice\n"));
    }
}



void
Split_operator_stepper_choice::construct_split_else()
{
    bool verbose=false;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    double t_turn, t_turn1;
    t_turn= MPI_Wtime();
    // check if all num_steps are >1
    for(List_choice_map::const_iterator lch_it=this->list_choice_map.begin(); lch_it!= this->list_choice_map.end();
                 ++lch_it){
            if (lch_it->second.num_steps < 1) {
                 std::cout<<"element name="<<  lch_it->first<<" num steps= "<<lch_it->second.num_steps<<std::endl;
                throw std::runtime_error(
                    "Split_operator_stepper_choice: all num_steps must be >= 1");
        }
    }

    if (this->num_steps_else<1) throw std::runtime_error(
                    "Split_operator_stepper_choice: num_steps_else must be >= 1");

    // find the length of the lattice without the chosen elements
    double length_else=0.;
    for (Lattice_elements::iterator it =
            get_lattice_simulator().get_lattice_sptr()->get_elements().begin(); it
            != get_lattice_simulator().get_lattice_sptr()->get_elements().end(); ++it){

        List_choice_map::const_iterator lch_it= this->list_choice_map.find((*it)->get_name());
        if (lch_it == this->list_choice_map.end()) length_else +=(*it)->get_length();

    }


    double step_length_else = length_else/this->num_steps_else;
    if (rank==0) std::cout<<" length_else="<< length_else<<" step_length_else="<< step_length_else<<std::endl;

    double  length_between=0.; // length between chosen elements
    Lattice_elements::iterator begin=get_lattice_simulator().get_lattice_sptr()->get_elements().begin();
    Lattice_elements::iterator end=get_lattice_simulator().get_lattice_sptr()->get_elements().begin();
    for (Lattice_elements::iterator it =
            get_lattice_simulator().get_lattice_sptr()->get_elements().begin(); it
            != get_lattice_simulator().get_lattice_sptr()->get_elements().end(); ++it){

        double length = (*it)->get_length();
        List_choice_map::const_iterator found_it= this->list_choice_map.find((*it)->get_name());
        if  (found_it != this->list_choice_map.end()){ // this is a chosen element
            // first make the step for "else" if length_else (i.e. the length between chosen elements) is nonzero
            if (begin != end) {
                make_stepper_else(begin, end, length_between, step_length_else);
            }
            if (length == 0.0) {
                Independent_operator_sptr
                        ind_op(
                                new Independent_operator(
                                        "step_zero_chosen",
                                        get_lattice_simulator().get_operation_extractor_map_sptr(),
                                        get_lattice_simulator().get_aperture_operation_extractor_map_sptr()));
                Lattice_element_slice_sptr slice(new Lattice_element_slice(*it));
                ind_op->append_slice(slice);
                Step_sptr step(new Step(0.0));
                step->append(ind_op, 1.0);
                get_steps().push_back(step);
                if ((rank==0) && verbose) std::cout<<" element: "<<(*it)->get_name()<<" zero length step"<<std::endl;
            }
            else{
               int steps_per_element=found_it->second.num_steps;
               double step_length= length/steps_per_element;
               double half_step_length = 0.5 * step_length;

                for (int i = 0; i <steps_per_element;  ++i) {
                    double left = i * step_length;
                    double middle = left + half_step_length;
                    double right = (i + 1) * step_length;

                    Step_sptr step(new Step(step_length));
                    //1st Half
                    Independent_operator_sptr
                        ind_op_first_half(
                                new Independent_operator(
                                        "first_half_chosen",
                                        get_lattice_simulator().get_operation_extractor_map_sptr(),
                                        get_lattice_simulator().get_aperture_operation_extractor_map_sptr()));
                    Lattice_element_slice_sptr slice_1st_half(
                            new Lattice_element_slice(*it, left, middle));
                    ind_op_first_half->append_slice(slice_1st_half);
                    step->append(ind_op_first_half, 0.5);


                     //Collective Effects
                    for (Collective_operators::const_iterator coll_op_it =
                        found_it->second.collective_operators.begin(); coll_op_it
                        != found_it->second.collective_operators.end();
                        ++coll_op_it) {
                    step->append(*coll_op_it, 1.0);
                    }

                        //2nd Half
                    Independent_operator_sptr
                        ind_op_second_half(
                                new Independent_operator(
                                        "second_half_chosen",
                                        get_lattice_simulator().get_operation_extractor_map_sptr(),
                                        get_lattice_simulator().get_aperture_operation_extractor_map_sptr()));
                    Lattice_element_slice_sptr slice_2nd_half(
                        new Lattice_element_slice(*it, middle, right));
                    ind_op_second_half->append_slice(slice_2nd_half);
                    step->append(ind_op_second_half, 0.5);

                    get_steps().push_back(step);

                }
            if ((rank==0) && verbose) std::cout<<" element: "<<(*it)->get_name()<<" nr of collective kicks ="<<steps_per_element<<std::endl;
            }
        ++end;
        begin=end;
        length_between=0.;
        }
        else{
                 ++end;
                 length_between +=length;
        }
    }

    if (begin != end) {
                make_stepper_else(begin, end, length_between, step_length_else);
     }
     get_lattice_simulator().set_slices(extract_slices(get_steps()));


    t_turn1= MPI_Wtime();
     if (rank==0){
          std::cout<<" stepper_choice split_else:"<<std::endl;
          std::cout<<"                total number of steps=  "<<this->get_steps().size()<<std::endl;
          std::cout<<" time: stepper choice done in "<<t_turn1-t_turn<<std::endl;
     }

}



Split_operator_stepper_choice::Split_operator_stepper_choice(
        Lattice_sptr lattice_sptr, int map_order, List_choice_map const & list_choice_map, bool split_else):
Stepper(lattice_sptr, map_order),    list_choice_map(list_choice_map)
{

   try{
        (this->list_choice_map.find("else")!= this->list_choice_map.end()) ?
            this->num_steps_else=this->list_choice_map["else"].num_steps:
            throw std::runtime_error(
                    "Split_operator_stepper_choice: if there is no  keyword \"else\" in the list_choice_map"
                     " you should use the constructor which provides \"num_steps_else\"");

        split_else ? construct_split_else():  construct_per_element_else();
    }
    catch (std::exception const& e){
        std::cout<<e.what()<<std::endl;
        MPI_Abort(MPI_COMM_WORLD, 222);
    }
}

Split_operator_stepper_choice::Split_operator_stepper_choice(int num_steps_else,
        Lattice_sptr lattice_sptr, int map_order, List_choice_map const & list_choice_map,  bool split_else):
Stepper(lattice_sptr, map_order),    list_choice_map(list_choice_map), num_steps_else(num_steps_else)
{

    try{
        split_else ? construct_split_else():  construct_per_element_else();
    }
    catch (std::exception const& e){
        std::cout<<e.what()<<std::endl;
        MPI_Abort(MPI_COMM_WORLD, 222);
    }

}

Split_operator_stepper_choice::Split_operator_stepper_choice(
                    Lattice_simulator const& lattice_simulator, List_choice_map const & list_choice_map, bool split_else):
Stepper(lattice_simulator),    list_choice_map(list_choice_map)
{

   try{
        (this->list_choice_map.find("else")!= this->list_choice_map.end()) ?
            this->num_steps_else=this->list_choice_map["else"].num_steps:
            throw std::runtime_error(
                    "Split_operator_stepper_choice: if there is no  keyword \"else\" in the list_choice_map"
                     " you should use the constructor which provides \"num_steps_else\"");

        split_else ? construct_split_else():  construct_per_element_else();
    }
    catch (std::exception const& e){
        std::cout<<e.what()<<std::endl;
        MPI_Abort(MPI_COMM_WORLD, 222);
    }
}

Split_operator_stepper_choice::Split_operator_stepper_choice(int num_steps_else,
                    Lattice_simulator const& lattice_simulator, List_choice_map const & list_choice_map,  bool split_else):
Stepper(lattice_simulator),    list_choice_map(list_choice_map), num_steps_else(num_steps_else)
{

    try{
        split_else ? construct_split_else():  construct_per_element_else();
    }
    catch (std::exception const& e){
        std::cout<<e.what()<<std::endl;
        MPI_Abort(MPI_COMM_WORLD, 222);
    }

}

Split_operator_stepper_choice::Split_operator_stepper_choice()
{
}

template<class Archive>
    void
    Split_operator_stepper_choice::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Stepper);
        ar & BOOST_SERIALIZATION_NVP(list_choice_map );
        ar & BOOST_SERIALIZATION_NVP(num_steps_else);
    }

template
void
Split_operator_stepper_choice::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Split_operator_stepper_choice::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Split_operator_stepper_choice::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Split_operator_stepper_choice::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Split_operator_stepper_choice::~Split_operator_stepper_choice()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Split_operator_stepper_choice);
