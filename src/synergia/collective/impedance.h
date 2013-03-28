#ifndef IMPEDANCE_H_
#define  IMPEDANCE_H_
#include "wake_field.h"
#include "synergia/simulation/operator.h"
#include "synergia/simulation/step.h"
#include "synergia/bunch/bunch.h"
#include "synergia/utils/serialization.h"
#include "synergia/utils/commxx.h"

 struct Bunch_properties
 {
     double x_mean;
     double y_mean;
     double z_mean;
     double realnum;
     int bucket_index;
     template<class Archive>
         void
         serialize(Archive & ar, const unsigned int version);
 };


class Impedance : public Collective_operator
{
private:
    Wake_field_sptr wake_field_sptr;
   
    int z_grid;
    int nstored_turns;
    int num_buckets;
    double orbit_length;
    double wake_factor;
    double bunch_spacing;
    bool full_machine;
    std::vector<int > wn;
    
    /// a list which contains informations about prevoius turns, the last element is the last (the earliest) turn stored
    /// it is updated at every step where impedance kick is applied
    /// the element is a vector of size num_bunches. The vector elements correspund to different bunches
    std::list< std::vector<Bunch_properties> > stored_vbunches;
    
    //boost::shared_ptr<Raw_MArray1d >  xmom_sptr;
    boost::shared_ptr<MArray1d >  xmom_sptr;
    boost::shared_ptr<MArray1d >  ymom_sptr;
    boost::shared_ptr<MArray1d >  zdensity_sptr;
    boost::shared_ptr<MArray1int > bin_partition_sptr;
    

    boost::shared_ptr<MArray1d > xwake_leading_sptr;
    boost::shared_ptr<MArray1d > xwake_trailing_sptr;
    boost::shared_ptr<MArray1d > ywake_leading_sptr;
    boost::shared_ptr<MArray1d > ywake_trailing_sptr;
    boost::shared_ptr<MArray1d > zwake0_sptr;
    /// xwake_leading=sum (wake_field_sptr->xw_lead)*xmom
    /// ywake_leading=sum (wake_field_sptr->yw_lead)*ymom   
    /// xwake_trailingl =sum (wake_field_sptr->xw_trail)
    /// zwake0= sum (wake_field_sptr->z_wake)
    /// Dpx=wake_factor*(xwake_leading+ xwake_trailingl*x)*time_step/(gamma*beta)
    /// Dpy=wake_factor*(ywake_leading+ ywake_trailingl*y)*time_step/(gamma*beta)
    /// Dpz=wake_factor*zwake0*time_step/(gamma*beta);    	   

    
    double N_factor;
    double cell_size_z;
    double bunch_z_mean;
    int bunch_bucket;


    void construct(); 
    void store_bunches_data(Bunch_train & bunch_train);
    void calculate_moments_and_partitions(Bunch & bunch);
    void calculate_kicks(Commxx_sptr const & comm_sptr);
    void  apply_impedance_kick(Bunch & bunch, double wake_factor);
  

  public:
    Impedance();
    Impedance(std::string const & wake_file, std::string const & wake_type, int const  & zgrid,
		  double const & orbit_length, double const & bunchsp, int const nstored_turns,
			  bool full_machine=false,std::vector<int > wn=std::vector<int >());

    Impedance(std::string const & wake_file, std::string const & wake_type, int const  & zgrid,
		  double const & orbit_length, int const & num_buckets, int const nstored_turns,
			  bool full_machine=false, std::vector<int > wn=std::vector<int >());
			  
    Impedance(Impedance const& impedance);
    virtual Impedance *
    clone();   
    
    void set_z_grid(int const  & zgrid);
    int get_z_grid() const;
    Wake_field_sptr 
    get_wake_field_sptr() const; 
    double get_orbit_length() const;
    double get_wake_factor() const;
    double get_bunch_spacing() const;
    int get_num_buckets() const;
    std::vector<int >  get_train_wave() const;

    
    MArray1d_ref &  get_xmom();
    MArray1d_ref &  get_ymom();
    MArray1d_ref & get_zdensity();
    MArray1int_ref &  get_bin_partition();
    MArray1d_ref const &  get_xmom() const;
    MArray1d_ref const &  get_ymom() const;
    MArray1d_ref const & get_zdensity() const;
    MArray1int_ref const &  get_bin_partition() const;
    
    MArray1d_ref &  get_xwake_leading(); 
    MArray1d_ref const & get_xwake_leading() const;    
    MArray1d_ref &  get_xwake_trailing();
    MArray1d_ref const & get_xwake_trailing() const;    
    MArray1d_ref &  get_ywake_leading();
    MArray1d_ref const & get_ywake_leading() const; 
    MArray1d_ref &  get_ywake_trailing();
    MArray1d_ref const &  get_ywake_trailing() const;     
    MArray1d_ref &  get_zwake0();
    MArray1d_ref const &  get_zwake0() const;
    
    std::list< std::vector<Bunch_properties> > &
    get_stored_vbunches();
    
    virtual bool 
    is_full_machine() const;
    virtual int 
    get_nstored_turns() const;
    virtual
    void  apply(Bunch & bunch, double time_step, Step & step, int verbosity, Logger & logger);
    virtual void
    apply(Bunch_train & bunch_train, double time_step, Step & step, int verbosity,
            Train_diagnosticss const& per_operation_train_diagnosticss, Logger & logger);
	    
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version); 
    virtual
    ~Impedance();
};
BOOST_CLASS_EXPORT_KEY(Impedance);
typedef boost::shared_ptr<Impedance> Impedance_sptr;

#endif /* IMPEDANCE_H_ */



