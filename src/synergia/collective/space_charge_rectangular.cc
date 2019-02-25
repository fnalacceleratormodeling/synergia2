#include "synergia/collective/space_charge_rectangular.h"
#include "deposit.h"
#include "synergia/utils/simple_timer.h"
#include "synergia/foundation/math_constants.h"
using mconstants::pi;
#include "synergia/foundation/physical_constants.h"
using pconstants::epsilon0;
#include <fftw3.h>
#include <fftw3-mpi.h>
#include "interpolate_rectangular_xyz.h"


Space_charge_rectangular::Space_charge_rectangular(Commxx_sptr comm_f_sptr, std::vector<double > const & pipe_size, 
			std::vector<int > const & grid_shape, bool equally_spread):
Collective_operator("space_charge_rectangular"),
pipe_size(pipe_size),
grid_shape(grid_shape),
domain_sptr(),
comm_f_sptr(comm_f_sptr),
fftw_helper_sptr(),
have_fftw_helper(),
have_domain(false),
equally_spread(equally_spread),
diagnostics_list(),
have_diagnostics(false)//,fftw_helper_sptr(),domain_sptr()
{

 try{  
    this->have_fftw_helper=false;
    construct_fftw_helper(comm_f_sptr);
     if ((!comm_f_sptr->has_this_rank()) && (equally_spread)) throw std::runtime_error(
		  "Space_charge_rectangular:: equally_spread is incompatible with this choice of comm_f_sptr ");

 }
 catch (std::exception const& e){
        std::cout<<e.what()<<std::endl;
        MPI_Abort(MPI_COMM_WORLD, 111);
    }
}


Space_charge_rectangular::Space_charge_rectangular(std::vector<double > const & pipe_size, std::vector<int > const & grid_shape):
Collective_operator("space_charge_rectangular"),  pipe_size(pipe_size),  grid_shape(grid_shape), 
have_domain(false),have_diagnostics(false)//,comm_f_sptr(),fftw_helper_sptr(),domain_sptr()
{

 try{   
    this->have_fftw_helper=false;
    this->equally_spread=false; 

 }
 catch (std::exception const& e){
        std::cout<<e.what()<<std::endl;
        MPI_Abort(MPI_COMM_WORLD, 111);
    }
}
 
Space_charge_rectangular::Space_charge_rectangular()
{
}

Space_charge_rectangular *
Space_charge_rectangular::clone()
{
    return new Space_charge_rectangular(*this);
}


template<class Archive>
    void
    Space_charge_rectangular::save(Archive & ar, const unsigned int version) const
    {
       
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Collective_operator);
        ar & BOOST_SERIALIZATION_NVP(comm_f_sptr)
	    &  BOOST_SERIALIZATION_NVP(grid_shape)
	    &  BOOST_SERIALIZATION_NVP(pipe_size) 
	    &  BOOST_SERIALIZATION_NVP(have_fftw_helper)
	    &  BOOST_SERIALIZATION_NVP(equally_spread)
        &  BOOST_SERIALIZATION_NVP(have_diagnostics)
        &  BOOST_SERIALIZATION_NVP(diagnostics_list);
    }

template<class Archive>
    void
    Space_charge_rectangular::load(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Collective_operator);
        ar & BOOST_SERIALIZATION_NVP(comm_f_sptr)
        &  BOOST_SERIALIZATION_NVP(grid_shape)
        &  BOOST_SERIALIZATION_NVP(pipe_size) 
        &  BOOST_SERIALIZATION_NVP(have_fftw_helper)
        &  BOOST_SERIALIZATION_NVP(equally_spread)
        &  BOOST_SERIALIZATION_NVP(have_diagnostics)
        &  BOOST_SERIALIZATION_NVP(diagnostics_list);
        domain_sptr.reset();
        have_domain=false;       
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
	 if (comm_sptr->has_this_rank()){
	    this->fftw_helper_sptr=  Fftw_rectangular_helper_sptr (  new   Fftw_rectangular_helper(grid_shape, comm_sptr));
	    this->have_fftw_helper=true;
	 }
	  this->comm_f_sptr=comm_sptr;
      }
      else {
	  throw std::runtime_error(
		  "Space_charge_rectangular::construct_fftw_helper:   already has fftw_helper ");
      }
    
}


void
Space_charge_rectangular::set_fftw_helper(Commxx_sptr comm_sptr, bool equally_spread)
{

 try{
    if (!have_fftw_helper){
	construct_fftw_helper(comm_sptr);
    }
    else {
	if (comm_sptr->has_this_rank()) this->fftw_helper_sptr->reset_comm_f(comm_sptr);
	this->comm_f_sptr=comm_sptr;	    
   }
   this->equally_spread=equally_spread; 
   if ((!comm_sptr->has_this_rank()) && (equally_spread)) throw std::runtime_error(
		  "Space_charge_rectangular:: set fftw: equally_spread is incompatible with this choice of comm_sptr ");
 }
 catch (std::exception const& e){
        std::cout<<e.what()<<std::endl;
        MPI_Abort(MPI_COMM_WORLD, 111);
    }
}

bool  
Space_charge_rectangular::get_have_fftw_helper() const
{
  return have_fftw_helper;
}


Commxx_sptr 
Space_charge_rectangular::get_comm_sptr() const
{
  return comm_f_sptr;
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

void
Space_charge_rectangular::add_diagnostics(Diagnostics_space_charge_rectangular_sptr ddiagnostics_sptr)
{
  diagnostics_list.push_back(ddiagnostics_sptr);
  this->have_diagnostics=true;
}

void
Space_charge_rectangular::set_diagnostics_list(Diagnostics_space_charge_rectangulars diagnostics_list)
{
  this->diagnostics_list =diagnostics_list;
  this->have_diagnostics=true;
}



Rectangular_grid_domain_sptr
Space_charge_rectangular::get_domain_sptr() const
{
return domain_sptr;
}

void
Space_charge_rectangular::set_domain(Bunch const & bunch)
{
  if (!have_domain){
        double beta = bunch.get_reference_particle().get_beta();
        std::vector<double >  dsize=pipe_size;
        dsize[2] /= beta;    // size in z_lab frame, longitudinal cdt coordinate
         // A.M physical_offsets of the domain should be rescaled too, but in this case they are zero    
        this->domain_sptr = Rectangular_grid_domain_sptr(
                      new Rectangular_grid_domain(dsize,  grid_shape, true));
         have_domain=true;  
  }
}

bool 
Space_charge_rectangular::get_have_domain() const
{
  return have_domain;
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

    if (!have_domain) set_domain(bunch);   
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


    return rho_sptr;
}



Distributed_rectangular_grid_sptr
Space_charge_rectangular::get_phi_local( Rectangular_grid & rho, double const& gamma)
{

    if (!have_fftw_helper)  throw std::runtime_error(
                "Space_charge_rectangular::get_phi_local  space_charge does not have have_fftw_helper defined");

    if (!comm_f_sptr->has_this_rank()) throw 
             std::runtime_error("space charge rectangular, get_phi_local, comm_f_sptr has no rank");
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

    for (int n=0; n < local_nx; ++n){
	int nt=n+1+local_x_start;
	for (int m=0; m< shape[1]; ++m){
	    int mt=m+1;
	    for (int p=0; p< shape[2]/2+1; ++p){
		double denominator=pi*pi*
		    (nt*nt/(pipe_size[0]*pipe_size[0])+mt*mt/(pipe_size[1]*pipe_size[1])+4.*p*p/(pipe_size[2]*pipe_size[2]*gamma*gamma));
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
  
   if (!comm_f_sptr->has_this_rank()) throw 
             std::runtime_error("space charge rectangular, fill_guards_pplanes, comm_f_sptr has no rank");

    int mpi_compare;
    MPI_Comm_compare(phi.get_comm().get(), comm_f_sptr->get(), &mpi_compare) ;
    if  (mpi_compare != MPI_IDENT)    {
        throw std::runtime_error("space charge rectangular, fill_guards_pplanes, phi comm is not the same as space_charge comm_f");
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
Space_charge_rectangular::get_En(Distributed_rectangular_grid &phi, int component)
{
  
    if (!comm_f_sptr->has_this_rank()) throw 
             std::runtime_error("space charge rectangular, get_En, comm_f_sptr has no rank");

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

    if (equally_spread){ 
	int error = MPI_Allgatherv(reinterpret_cast<void*>(En_local_a.origin()),
		  receive_counts[lrank], MPI_DOUBLE,
		  reinterpret_cast<void*>(En->get_grid_points().origin()),
					  &receive_counts[0], &receive_offsets[0], MPI_DOUBLE, comm_f_sptr->get());

	if (error != MPI_SUCCESS) {
	    throw std::runtime_error(
		"MPI error in Space_charge_rectangular(MPI_Allgatherv in get_En: En_local)");
	}
	En->set_normalization(phi.get_normalization()); // we should have here  \div $\vec{E}=rho/epsilon   
    }
    else{
	int error = MPI_Gatherv(reinterpret_cast<void*>(En_local_a.origin()),
		  receive_counts[lrank], MPI_DOUBLE,
		  reinterpret_cast<void*>(En->get_grid_points().origin()),
					  &receive_counts[0], &receive_offsets[0], MPI_DOUBLE, 0,comm_f_sptr->get());

	if (error != MPI_SUCCESS) {
	    throw std::runtime_error(
		"MPI error in Space_charge_rectangular(MPI_Gatherv in get_En: En_local)");
	}
	
    }   
    
   //AM!  make sure the field is zero at the edge of the grid
// THIS toghether with zero charge distribution at the edge of the grid is essential for a conservative approximation
    MArray3d_ref grid_points(En->get_grid_points());
    for (int j=0; j<grid_points.shape()[1];++j){
        for (int k=0; k<grid_points.shape()[2];++k){
            grid_points[0][j][k]=0.;
            grid_points[grid_points.shape()[0]-1][j][k]=0.;
            
        }
    }    
    for (int i=0; i<grid_points.shape()[0];++i){
        for (int k=0; k<grid_points.shape()[2];++k){
            grid_points[i][0][k]=0.;         
            grid_points[i][grid_points.shape()[1]-1][k]=0.;
        }
    } 
    
    
    
    t = simple_timer_show(t, "get_En:  gather En");
     
    
    
    
    return En;

}

void
Space_charge_rectangular::do_diagnostics(Rectangular_grid const& En, int component, double time_step, Step & step, 
                                          Bunch & bunch)
{   
   if (have_diagnostics) {
      if ((component==0) || (component==1)){
         double step_beta=step.get_betas()[component];
         for (Diagnostics_space_charge_rectangulars::const_iterator d_it = diagnostics_list.begin();
            d_it != diagnostics_list.end(); ++d_it){
            if (bunch.is_bucket_index_assigned()){
                if ((*d_it)->get_bunch().get_bucket_index()==bunch.get_bucket_index()){                  
                    (*d_it)->update(bunch, En, component, time_step, step_beta); 
                    if (component==1) (*d_it)->write();
                }
            }
            else{
                    (*d_it)->update(bunch, En, component, time_step, step_beta); 
                    if (component==1) (*d_it)->write();
            }                           
         }
      }    
   } 
   
}  


void
Space_charge_rectangular::apply_kick(Bunch & bunch, Rectangular_grid const& En, double  delta_t, int component)
{
  
 //AM: kicks  in the z_lab frame 
 //Delta p_x&=& F_x \Delta t&=& - q \frac{1}{\gamma^2} \frac{\partial \Phi'}{\partial x} \Delta t=q \frac{1}{\beta \gamma^2} E_{grid~x} \Delta t\\
 //Delta E &= & q E_z \Delta s&=& q \frac{1}{\gamma^2 \beta} \frac{\partial \Phi'}{\partial ct} \beta c\Delta t=-q \frac{c}{\beta \gamma^2 }E_{grid~z} \Delta t\\
 // 1/beta factor in E_grid from charge deposition on (x,y,cdt) coordinates grid 

 
 

    bunch.convert_to_state(Bunch::fixed_z_lab);
    double q = bunch.get_particle_charge() * pconstants::e; // [C]
    double gamma=bunch.get_reference_particle().get_gamma();
    double beta=bunch.get_reference_particle().get_beta();
// unit_conversion: [kg m/s] to [Gev/c] 
    double unit_conversion = pconstants::c / (1.0e9 * pconstants::e);
// scaled p = p/p_ref
    double p_ref=bunch.get_reference_particle().get_momentum();
    double factor = unit_conversion * q * delta_t* En.get_normalization()/
            (p_ref*gamma*gamma*beta); // transverse kicks
   
      
    int ps_component = 2 * component + 1;
  
 
    Rectangular_grid_domain & domain(*En.get_domain_sptr());
    MArray3d_ref grid_points(En.get_grid_points());
   
    if (component==2){
       factor *= -p_ref; 
       double m = bunch.get_mass();
       for (int part = 0; part < bunch.get_local_num(); ++part) {
          double x = bunch.get_local_particles()[part][Bunch::x];
          double y = bunch.get_local_particles()[part][Bunch::y];
          double z = bunch.get_local_particles()[part][Bunch::z];   
          double grid_val =  interpolate_rectangular_xyz(x, y, z, domain,
                grid_points);
          double p=p_ref +bunch.get_local_particles()[part][Bunch::dpop] * p_ref;        
          double Eoc_i = std::sqrt(p * p + m * m);
          double Eoc_f=  Eoc_i + factor * grid_val;
          double delta_dpop=(std::sqrt(Eoc_f*Eoc_f-m*m)-std::sqrt(Eoc_i*Eoc_i-m*m))/p_ref;
          bunch.get_local_particles()[part][ps_component] += delta_dpop;
       }      
    }
    else{
        for (int part = 0; part < bunch.get_local_num(); ++part) {
          double x = bunch.get_local_particles()[part][Bunch::x];
          double y = bunch.get_local_particles()[part][Bunch::y];
          double z = bunch.get_local_particles()[part][Bunch::z];        
          double grid_val = interpolate_rectangular_xyz(x, y, z, domain,
                  grid_points);   
          bunch.get_local_particles()[part][ps_component] += factor * grid_val; 
        }
    }
   
}

std::vector<Rectangular_grid_sptr>
Space_charge_rectangular::get_Efield(Rectangular_grid & rho,Bunch const& bunch, int max_component, double const & gamma )
{	
   if (equally_spread) throw std::runtime_error
              	("Space_charge_rectangular get_Efield: don't call this function for true equally_spread ");
  
   std::vector<Rectangular_grid_sptr> Efield;
    if (comm_f_sptr->has_this_rank()){
       Distributed_rectangular_grid_sptr phi_local(get_phi_local(rho,gamma)); // \nabla phi= -rho/epsilon0; [phi]=kg*m^2*C^{-1}*s^{-2}            
       for (int component = 0; component < max_component; ++component) {
          Efield.push_back(get_En(*phi_local, component));
       }	
    } 
    else{
      for (int component = 0; component < max_component; ++component){ 
        Efield.push_back(Rectangular_grid_sptr(new Rectangular_grid(domain_sptr)));			
      }
      int mpi_compare;
      MPI_Comm_compare(comm_f_sptr->get_parent_sptr()->get(), bunch.get_comm_sptr()->get(), &mpi_compare);
      if ((mpi_compare != MPI_IDENT) && ( mpi_compare != MPI_CONGRUENT)){	       
                  throw std::runtime_error
                  ("Space_charge_rectangular get_Efield: comm_f_sptr parent and bunch.comm are not congruent");
      } 
    }  // comm_f_sptr->has_this_rank()
   
   
   
   

   // cast Efield from rank=0 of comm_spc to whole bunch.get_comm
   std::vector<int > shape(domain_sptr->get_grid_shape());
   int count=shape[0]*shape[1]*shape[2];
   for (int component = 0; component < max_component; ++component){
      int error=MPI_Bcast(Efield[component]->get_grid_points().origin(), count, 
				  MPI_DOUBLE, 0, bunch.get_comm().get());
	if (error != MPI_SUCCESS) {
		throw std::runtime_error(
		  "MPI error in Space_charge_rectangular, get_Efield: MPI_Bcast Efield failed)");
	}	  
	Efield[component]->set_normalization(1./(4.*shape[0]*shape[1]*shape[2]*epsilon0));
   }    
   return Efield;   
}  


void
Space_charge_rectangular::apply(Bunch & bunch, double time_step, Step & step, int verbosity, Logger & logger)
{ 
   // t = simple_timer_current();

    double t = 0.0;
    double t1 = 0.0;
    bunch.convert_to_state(Bunch::fixed_z_lab);
    set_domain(bunch);
    Rectangular_grid_sptr rho_sptr(get_charge_density(bunch)); // [C/m^3], [Q/DxDyD(ct)] in z_lab frame
    t = simple_timer_show(t, "sc_apply: get-rho");
    double gamma=bunch.get_reference_particle().get_gamma();
    int max_component(3);   
    if (equally_spread){     
        if (comm_f_sptr->get_parent_sptr().get()!= 0) {
	    int mpi_compare;
	    MPI_Comm_compare(comm_f_sptr->get_parent_sptr()->get(), bunch.get_comm_sptr()->get(), &mpi_compare);
            if ((mpi_compare != MPI_IDENT) && ( mpi_compare != MPI_CONGRUENT)){
               throw std::runtime_error
               ("Space_charge_rectangular apply, equally_spread=1: comm_f_sptr parent and bunch.comm are not congruent");
            } 	  
	    }	   
        Distributed_rectangular_grid_sptr phi_local(get_phi_local(*rho_sptr,gamma)); // Phi_local is Phi'(\gamma z)=\Phi(z), but the grid is (x,y,cdt)
	    for (int component = 0; component < max_component; ++component) {	  
	        Rectangular_grid_sptr  En(get_En(*phi_local, component)); // E=-/grad phi; [E]=kg*m/(C*s^2)=N/C	
            do_diagnostics(*En,component, time_step,step, bunch);
	        apply_kick(bunch, *En, time_step, component);	 
	   }
    }
    else{
        std::vector<Rectangular_grid_sptr>  Efield(get_Efield(*rho_sptr, bunch, max_component,gamma));
        for (int component = 0; component < max_component; ++component) { 
            do_diagnostics(*(Efield[component]),component, time_step,step, bunch);
            apply_kick(bunch, *(Efield[component]), time_step, component);  
        }
    }

    t = simple_timer_show(t, "sc_apply: 3x apply_kick and get En");
    t1 = simple_timer_show(t1, "sc_aplly total");
}

