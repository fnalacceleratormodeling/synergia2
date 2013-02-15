#include "synergia/collective/space_charge_rectangular.h"
#include "deposit.h"
#include "synergia/utils/simple_timer.h"
#include "synergia/foundation/math_constants.h"
using mconstants::pi;
#include "synergia/foundation/physical_constants.h"
using pconstants::epsilon0;
#include <fftw3.h>
#include <fftw3-mpi.h>


Space_charge_rectangular::Space_charge_rectangular(Commxx_sptr comm_f_sptr, std::vector<double > const & pipe_size, std::vector<int > const & grid_shape):
Collective_operator("space_charge_rectangular"), pipe_size(pipe_size), grid_shape(grid_shape),  comm_f_sptr(comm_f_sptr)
{

 try{
    this->domain_sptr = Rectangular_grid_domain_sptr(
                    new Rectangular_grid_domain(pipe_size, grid_shape , true));
    this->have_fftw_helper=false;
    construct_fftw_helper(comm_f_sptr);

 }
 catch (std::exception const& e){
        std::cout<<e.what()<<std::endl;
        MPI_Abort(MPI_COMM_WORLD, 111);
    }
}


Space_charge_rectangular::Space_charge_rectangular(std::vector<double > const & pipe_size, std::vector<int > const & grid_shape):
Collective_operator("space_charge_rectangular"),  pipe_size(pipe_size),  grid_shape(grid_shape)
{

 try{
    this->domain_sptr = Rectangular_grid_domain_sptr(
                    new Rectangular_grid_domain(pipe_size,  grid_shape , true));
    this->have_fftw_helper=false;

 }
 catch (std::exception const& e){
        std::cout<<e.what()<<std::endl;
        MPI_Abort(MPI_COMM_WORLD, 111);
    }
}
 
Space_charge_rectangular::Space_charge_rectangular()
{
}









template<class Archive>
    void
    Space_charge_rectangular::save(Archive & ar, const unsigned int version) const
    {
       
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Collective_operator);
        ar & BOOST_SERIALIZATION_NVP(comm_f_sptr)
                & BOOST_SERIALIZATION_NVP(grid_shape)
                & BOOST_SERIALIZATION_NVP(pipe_size) 
                & BOOST_SERIALIZATION_NVP(have_fftw_helper);
    }

template<class Archive>
    void
    Space_charge_rectangular::load(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Collective_operator);
        ar & BOOST_SERIALIZATION_NVP(comm_f_sptr)
                & BOOST_SERIALIZATION_NVP(grid_shape)
                & BOOST_SERIALIZATION_NVP(pipe_size)
                & BOOST_SERIALIZATION_NVP(have_fftw_helper); 
     
        domain_sptr = Rectangular_grid_domain_sptr(
                    new Rectangular_grid_domain(pipe_size, grid_shape , true));
        if (have_fftw_helper) { this->have_fftw_helper=false;          
                                construct_fftw_helper(comm_f_sptr);  
                               }
                                       
    }

template
void
Space_charge_rectangular::save<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version) const;
template
void
Space_charge_rectangular::save<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version) const;

template
void
Space_charge_rectangular::load<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);
template
void
Space_charge_rectangular::load<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);





BOOST_CLASS_EXPORT_IMPLEMENT(Space_charge_rectangular)













Space_charge_rectangular::~Space_charge_rectangular()
{
 //fftw_mpi_cleanup();
}

void
Space_charge_rectangular::construct_fftw_helper(Commxx_sptr comm_sptr)
{
    if (!have_fftw_helper){
        this->fftw_helper_sptr=  Fftw_rectangular_helper_sptr (  new   Fftw_rectangular_helper(grid_shape, comm_sptr));
        this->comm_f_sptr=comm_sptr;
        this->have_fftw_helper=true;
    }
    else {
        throw std::runtime_error(
                "Space_charge_rectangular::construct_fftw_helper:   already has fftw_helper ");
    }
}


void
Space_charge_rectangular::set_fftw_helper(Commxx_sptr comm_sptr)
{

 try{
    if (!have_fftw_helper){
            construct_fftw_helper(comm_sptr);
        }
        else {
            this->fftw_helper_sptr->reset_comm_f(comm_sptr);
            this->comm_f_sptr=comm_sptr;
        }
    }
 catch (std::exception const& e){
        std::cout<<e.what()<<std::endl;
        MPI_Abort(MPI_COMM_WORLD, 111);
    }
}

std::vector<double >
Space_charge_rectangular::get_pipe_size() const
{
return pipe_size;
}

std::vector<int >
Space_charge_rectangular::get_grid_shape() const
{
return grid_shape;
}

Rectangular_grid_domain_sptr
Space_charge_rectangular::get_domain_sptr() const
{
return domain_sptr;
}


Fftw_rectangular_helper_sptr
Space_charge_rectangular::get_fftw_helper_sptr() const
{
    return fftw_helper_sptr;
}



Rectangular_grid_sptr
Space_charge_rectangular::get_charge_density(Bunch const& bunch)
{
    double t;
    t = simple_timer_current();

    Rectangular_grid_sptr rho_sptr(new Rectangular_grid(domain_sptr));
    deposit_charge_rectangular_xyz(*rho_sptr, bunch);
    //t = simple_timer_show(t, "sc_get_charge_density: depozit_xyz");

    t = simple_timer_current();
    int error = MPI_Allreduce(MPI_IN_PLACE, (void*)  rho_sptr->get_grid_points().origin(),
                               rho_sptr->get_grid_points().num_elements(), MPI_DOUBLE, MPI_SUM, bunch.get_comm().get());


     t = simple_timer_show(t, "sc_get_charge_density: allmpireduce");

    if (error != MPI_SUCCESS) {
        throw std::runtime_error(
                "MPI error in Space_charge_rectangular: MPI_Allreduce in get_charge_density");
    }

//     double global_normalization;
//     double local_normalization(rho.get_normalization());
//     error = MPI_Allreduce(&local_normalization, &global_normalization,
//                                1, MPI_DOUBLE, MPI_SUM, comm.get());
//     if (error != MPI_SUCCESS) {
//         throw std::runtime_error(
//                 "MPI error in Space_charge_rectangular(MPI_Allreduce in get_phi: global_rho normalization)");
//     }
   //rho.set_normalization(global_normalization);

    return rho_sptr;
}

Distributed_rectangular_grid_sptr
Space_charge_rectangular::get_phi_local( Rectangular_grid & rho, Bunch const& bunch)
{

    if (!have_fftw_helper)  throw std::runtime_error(
                "Space_charge_rectangular::get_phi_local  space_charge does not have have_fftw_helper defined");

//    double t;
//    t = simple_timer_current();


    MArray3d_ref rho_ref(rho.get_grid_points());
    std::vector<int > shape(domain_sptr->get_grid_shape());

    ptrdiff_t local_nx=get_fftw_helper_sptr()->get_local_nx();
    ptrdiff_t local_x_start=get_fftw_helper_sptr()->get_local_x_start();



    int lower=local_x_start;
    int upper=local_x_start+local_nx;
    std::string solver("rectangular");
    std::vector<double > local_physical_size(domain_sptr->get_physical_size());
    std::vector<double > local_physical_offset(domain_sptr->get_physical_offset());
     // local_physical_size[0]=domain_sptr->get_cell_size()[0]*local_nx;
     //local_physical_offset[0]=domain_sptr->get_left()[0]+local_x_start*domain_sptr->get_cell_size()[0];
    std::vector<int > shape_local(shape);
    shape_local[0]=local_nx;
    Distributed_rectangular_grid_sptr phi_local(
        new Distributed_rectangular_grid(local_physical_size, local_physical_offset, shape_local,
            true, lower, upper, comm_f_sptr, solver)
            );


    MArray3d_ref phi_local_ref(phi_local->get_grid_points());
    MArray3d_ref rho_local_ref(rho_ref.origin()+local_x_start*shape[1]*shape[2],boost::extents[local_nx][shape[1]][shape[2]]);

    // t = simple_timer_current();
    get_fftw_helper_sptr()->transform(rho_local_ref, phi_local_ref);

   // t = simple_timer_show(t, "sc_get_phi_local: fftw_dst_direct");





    const int memory_fudge_factor = 1;
    fftw_complex *rho_nmp_local;
    rho_nmp_local= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) *local_nx*shape[1]*(shape[2]/2+1)*memory_fudge_factor);
    int dim[] = {shape[2]};
    fftw_plan plan=fftw_plan_many_dft_r2c(1, dim, local_nx*shape[1], phi_local_ref.origin(), NULL,
                                       1, shape[2], rho_nmp_local, NULL, 1,shape[2]/2+1,
                                        FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
   // t = simple_timer_show(t, "sc_get_phi_local: fftw1d_direct");


    MArray3dc_ref rho_nmp_ref(reinterpret_cast<std::complex<double>*>(rho_nmp_local), boost::extents[local_nx][shape[1]][shape[2]/2+1]);
    std::vector<double > pipe_size(rho.get_domain().get_physical_size());
    for (int n=0; n < local_nx; ++n){
         int nt=n+1+local_x_start;
         for (int m=0; m< shape[1]; ++m){
             int mt=m+1;
             for (int p=0; p< shape[2]/2+1; ++p){
                 double denominator=pi*pi*
                     (nt*nt/(pipe_size[0]*pipe_size[0])+mt*mt/(pipe_size[1]*pipe_size[1])+4.*p*p/(pipe_size[2]*pipe_size[2]));
                        rho_nmp_ref[n][m][p] /= denominator; // delta Phi =- rho
             }
         }
    }
    //t = simple_timer_show(t, "sc_get_phi_local: loop_phi_nmp");

    plan=fftw_plan_many_dft_c2r(1, dim, local_nx*shape[1], rho_nmp_local, NULL, 1,shape[2]/2+1,
                                        phi_local_ref.origin(), NULL, 1,shape[2],
                                        FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    fftw_free(rho_nmp_local);


    //t = simple_timer_show(t, "sc_get_phi_local: fftw1d_inverse");
    //t = simple_timer_current();
    get_fftw_helper_sptr()->inv_transform(phi_local_ref,phi_local_ref);
   // t = simple_timer_show(t, "sc_get_phi_local: fftw_dst_inverse");


    phi_local->set_normalization(1./(4.*shape[0]*shape[1]*shape[2]*epsilon0));


    return phi_local;
 }



void
Space_charge_rectangular::fill_guards_pplanes(Distributed_rectangular_grid & phi, int lower, int upper, int lengthx,
                          MArray2d & g_lower, MArray2d &g_upper)
{
    int mpi_compare;
    MPI_Comm_compare(phi.get_comm().get(), comm_f_sptr->get(), &mpi_compare) ;
    if  (mpi_compare != MPI_IDENT)    {
        throw std::runtime_error("space charge rectangular, phi comm is not the same as space_charge comm_f");
    }

    int lrank=comm_f_sptr->get_rank();
    std::vector<int > shape_phi(phi.get_domain().get_grid_shape());
    int message_size = shape_phi[1] * shape_phi[2];
    int shapex=upper-lower;

    void *send_buffer, *recv_buffer;
    MPI_Status status;

     // send to the right

    if ((upper < lengthx) &&  (upper >0)) {
        send_buffer=reinterpret_cast<void*>(phi.get_grid_points().origin()+(shapex-1)*message_size);
        MPI_Send(send_buffer, message_size, MPI_DOUBLE, lrank + 1, lrank, comm_f_sptr->get());
    }
    if (lower > 0) {
        recv_buffer=reinterpret_cast<void*>(g_lower.data());
        MPI_Recv(recv_buffer, message_size, MPI_DOUBLE, lrank - 1, lrank - 1,
                 comm_f_sptr->get(), &status);
    }
 // send to the left

    if (lower > 0) {
        send_buffer=reinterpret_cast<void*>(phi.get_grid_points().origin());
        MPI_Send(send_buffer, message_size, MPI_DOUBLE, lrank - 1, lrank,
                 comm_f_sptr->get());
    }
    if ((upper < lengthx) &&  (upper >0)){
        recv_buffer=reinterpret_cast<void*>(g_upper.data());
        MPI_Recv(recv_buffer, message_size, MPI_DOUBLE, lrank + 1, lrank + 1,
                 comm_f_sptr->get(), &status);
    }
}


Rectangular_grid_sptr
Space_charge_rectangular::get_En(Distributed_rectangular_grid &phi, Bunch const& bunch, int component)
{
    if ((component < 0) || (component > 2)) {
        std::stringstream message("");
        message << "calculate_E_n: invalid argument component=" << component
                << ". Argument be in range 0<=component<=2";
        throw std::invalid_argument(message.str());
    }
    double t;
    t = simple_timer_current();


    std::vector<int > shape_phi(phi.get_domain().get_grid_shape());
    MArray3d E_local(boost::extents[shape_phi[0]][shape_phi[1]][shape_phi[2]]);
    std::vector<double > hi(domain_sptr->get_cell_size());
    double h(hi[component]);


    int size=comm_f_sptr->get_size();
    int lrank=comm_f_sptr->get_rank();
    std::vector<int > shape(domain_sptr->get_grid_shape());

    int lower =  phi.get_lower();
    int upper =  phi.get_upper();
    int lengthx=shape[0];

    MArray2d guard_lower(boost::extents[shape[1]][shape[2]]);
    MArray2d guard_upper(boost::extents[shape[1]][shape[2]]);


    Rectangular_grid_domain_sptr domain_local_sptr(phi.get_domain_sptr());
    Rectangular_grid_sptr En_local(new Rectangular_grid(domain_local_sptr));
    MArray3d_ref En_local_a(En_local->get_grid_points());
    MArray3d_ref phi_a(phi.get_grid_points());
    boost::array<MArray3d::index, 3 > center, left, right;


    if (component==0) {
        fill_guards_pplanes(phi,  lower, upper, lengthx, guard_lower, guard_upper);
        for (int i = 0; i < upper-lower; ++i) {
                left[0] = i;
                center[0] = i;
                right[0] = i;
                for (int j = 0; j < shape[1]; ++j) {
                    left[1] = j;
                    center[1] =j;
                    right[1] = j;
                    for (int k = 0; k < shape[2]; ++k) {
                        left[2] = k;
                        center[2] = k;
                        right[2] = k;
                        right[component] = std::min(int(center[component] + 1), upper-lower-1);
                        left[component] = std::max(int(center[component] - 1),0);
                        double delta=  (right[component]-left[component])*h;
                        double phi_right, phi_left;
                        if ((center[0]==upper-lower-1) && (upper < shape[0])) {
                            phi_right=guard_upper[j][k];
                            delta=2*h;
                        }
                        else  phi_right=phi_a(right);
                        if ((center[0]==0) && (lower>0)){
                           phi_left=guard_lower[j][k];
                           delta=2*h;
                        }
                        else phi_left=phi_a(left);
                         //  $\vec{E} = - \grad \phi$
                        En_local_a(center)=-(phi_right - phi_left) / delta;
                     }
                 }
            }


    }
    else if (component==2) { // we have periodicity
        for (int i = 0; i < upper-lower; ++i) {
            left[0] = i;
            center[0] = i;
            right[0] = i;
            for (int j = 0; j < shape[1]; ++j) {
                left[1] = j;
                center[1] =j;
                right[1] = j;
                for (int k = 0; k < shape[2]; ++k) {
                    left[2] = k;
                    center[2] = k;
                    right[2] = k;
                    right[component] = (center[component] + 1) % shape[2];
                    left[component] = (center[component]-1 + shape[2])% shape[2];
                    double delta=2.*h;
                      //  $\vec{E} = - \grad \phi$
                    En_local_a(center)= -(phi_a(right) - phi_a(left)) / delta;
                 }
             }
        }
    }
    else {
        for (int i = 0; i < upper-lower; ++i) {
            left[0] = i;
            center[0] = i;
            right[0] = i;
            for (int j = 0; j < shape[1]; ++j) {
                left[1] = j;
                center[1] =j;
                right[1] = j;
                for (int k = 0; k < shape[2]; ++k) {
                    left[2] = k;
                    center[2] = k;
                    right[2] = k;
                    right[component] = std::min(int(center[component] + 1),shape[component] - 1);
                    left[component] = std::max(int(center[component] - 1),0);
                    double delta=(right[component]-left[component])*h;
                      //  $\vec{E} = - \grad \phi$
                    En_local_a(center)= -(phi_a(right) - phi_a(left)) / delta;
                 }
             }
        }
    }

  //  t = simple_timer_show(t, "get_En: calculate E local");

    Rectangular_grid_sptr En(new Rectangular_grid(domain_sptr));

     std::vector<int> uppers(phi.get_uppers());
     std::vector<int> receive_counts(phi.get_lengths()), receive_offsets(size);
     for (int i=0; i< size; ++i) {
         receive_offsets.at(i) = uppers.at(i)*shape[1]*shape[2]-receive_counts.at(i);
     }

    t = simple_timer_current();
    int error = MPI_Allgatherv(reinterpret_cast<void*>(En_local_a.origin()),
               receive_counts[lrank], MPI_DOUBLE,
               reinterpret_cast<void*>(En->get_grid_points().origin()),
                                       &receive_counts[0], &receive_offsets[0], MPI_DOUBLE, comm_f_sptr->get());

    if (error != MPI_SUCCESS) {
        throw std::runtime_error(
            "MPI error in Space_charge_rectangular(MPI_Allgatherv in get_En: En_local)");
    }

    t = simple_timer_show(t, "get_En:  gather En");
    En->set_normalization(phi.get_normalization()); // we should have here  \div $\vec{E}=rho/epsilon
    return En;

}

void
Space_charge_rectangular::apply_kick(Bunch & bunch, Rectangular_grid const& En, double  delta_t, int component)
{
    /// En is electric field in units of N/C
    double q = bunch.get_particle_charge() * pconstants::e; // [C]
    // delta_t_beam: [s] in beam frame
    double delta_t_beam = delta_t / bunch.get_reference_particle().get_gamma();
    // unit_conversion: [kg m/s] to [Gev/c]
    double unit_conversion = pconstants::c / (1.0e9 * pconstants::e);
    // scaled p = p/p_ref
    double p_scale = 1.0 / bunch.get_reference_particle().get_momentum();
    double factor = unit_conversion * q * delta_t_beam * En.get_normalization()
            * p_scale;

    int ps_component = 2 * component + 1;
    Rectangular_grid_domain & domain(*En.get_domain_sptr());
    MArray3d_ref grid_points(En.get_grid_points());
    for (int part = 0; part < bunch.get_local_num(); ++part) {
        double x = bunch.get_local_particles()[part][Bunch::x];
        double y = bunch.get_local_particles()[part][Bunch::y];
        double z = bunch.get_local_particles()[part][Bunch::z];
        double grid_val = interpolate_rectangular_xyz(x, y, z, domain,
                grid_points);
        bunch.get_local_particles()[part][ps_component] += factor * grid_val;
    }
}


void
Space_charge_rectangular::apply(Bunch & bunch, double time_step, Step & step, int verbosity, Logger & logger)
{
    double t,t1;
    t = simple_timer_current();
    t1 = simple_timer_current();

    bunch.convert_to_state(Bunch::fixed_t_bunch);
    t = simple_timer_show(t, "sc_apply: convert-to-state");
    Rectangular_grid_sptr rho(get_charge_density(bunch)); // [C/m^3]
    t = simple_timer_show(t, "sc_apply: get-rho");

    Distributed_rectangular_grid_sptr phi_local(get_phi_local(*rho, bunch)); // \nabla phi= -rho/epsilon0; [phi]=kg*m^2*C^{-1}*s^{-2}
    t = simple_timer_show(t, "sc_apply: get-phi_local");

    int max_component;
    max_component = 3;

    for (int component = 0; component < max_component; ++component) {
       Rectangular_grid_sptr  En(get_En(*phi_local, bunch, component)); // E=-/grad phi; [E]=kg*m/(C*s^2)=N/C
     //  t = simple_timer_show(t, "sc_apply: get_En");
       apply_kick(bunch, *En, time_step, component);
       // t = simple_timer_show(t, "sc_apply: apply_kick");
    }
     t = simple_timer_show(t, "sc_apply: 3x apply_kick and get En");
     t1 = simple_timer_show(t1, "sc_aplly total");
}

