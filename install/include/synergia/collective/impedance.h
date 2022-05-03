#ifndef IMPEDANCE_H_
#define IMPEDANCE_H_

#include "wake_field.h"

#include "synergia/simulation/operator.h"
#include "synergia/simulation/collective_operator_options.h"

class Impedance;

struct Impedance_options 
    : public CO_base_options<Impedance_options, Impedance>
{
    std::string wake_file;
    std::string wake_type;

    int z_grid;
    bool full_machine;

    int nstored_turns;
    int num_buckets;
    double orbit_length;
    double bunch_spacing;

    std::array<int, 3> wn;

    Impedance_options(
            std::string const& wake_file = "",
            std::string const& wake_type = "",
            int z_grid = 1000 )
        : wake_file(wake_file)
        , wake_type(wake_type)
        , z_grid(z_grid)
        , full_machine(false)
        , nstored_turns(15)
        , num_buckets(1)
        , orbit_length(1)
        , bunch_spacing(1)
    { }

    template<class Archive>
    void serialize(Archive & ar)
    { 
        ar(cereal::base_class<CO_base_options>(this));

        ar(wake_file);
        ar(wake_type);
        ar(z_grid);
        ar(full_machine);
        ar(nstored_turns);
        ar(num_buckets);
        ar(orbit_length);
        ar(bunch_spacing);
    }
};


CEREAL_REGISTER_TYPE(Impedance_options)

struct Bunch_props
{
    Bunch_props(int num_bunches, int max_turns)
        : num_bunches(num_bunches)
        , max_turns(max_turns)
        , registered_turns(0)
        , write_pos(0)
        , xmean("xmean", num_bunches*max_turns)
        , ymean("ymean", num_bunches*max_turns)
        , zmean("zmean", num_bunches*max_turns)
        , realnum("realnum", num_bunches*max_turns)
        , bucket_index("bucket", num_bunches*max_turns)
    { }

    // done writing a turn. increment the pointers
    void increment_registered_turns()
    {
        if (++registered_turns >= max_turns) registered_turns = max_turns;
        if (++write_pos >= max_turns) write_pos = 0;
    }

    // write position for current turn and the specified bunch
    KOKKOS_INLINE_FUNCTION
    int get_write_index(int bunch) const
    { return write_pos*num_bunches + bunch; }

    // turn: 0 is current, -1 is prev, -2 ...
    KOKKOS_INLINE_FUNCTION
    int get_read_index(int turn, int bunch) const
    { 
        // read_pos is one behind the write_pos
        int pos = write_pos + turn - 1;
        if (pos < 0) pos += max_turns;
        return pos*num_bunches + bunch; 
    }

    int num_bunches;
    int max_turns;

    int registered_turns;
    int write_pos;

    karray1d_dev xmean;
    karray1d_dev ymean;
    karray1d_dev zmean;
    karray1d_dev realnum;
    Kokkos::View<int*> bucket_index;
};

struct Bunch_params
{
    double z_mean;
    double z_left;
    double cell_size_z;
    double N_factor;
    int bucket;
};

class Impedance : public Collective_operator
{
private:

    const Impedance_options opts;
    std::string bunch_sim_id;

    // bunch properties over turns
    Bunch_props bps;

    // z_grid*3, in the fortran order for
    // zdensity, xmom, ymom
    karray1d_dev zbinning;
    karray1d_hst h_zbinning;

    // buffer for wake fields
    // z_grid*5, in the fortran order for
    // xwake_leading, xwake_trailing,
    // ywake_leading, ywake_trailing,
    // zwake0
    karray1d_dev wakes;
    karray1d_hst h_wakes;

    Wake_field wake_field;

private:

    void apply_impl(
            Bunch_simulator& simulator, 
            double time_step, 
            Logger& logger) override;

    void apply_bunch(
            Bunch& bunch, 
            double time_step, 
            Logger& logger);

    void construct_workspaces(
            Bunch_simulator const& sim);

    void store_bunches_data(
            Bunch_simulator const& sim);

    Bunch_params
    calculate_moments_and_partitions(
            Bunch const& bunch);

    void calculate_kicks(
            Bunch const& bunch,
            Bunch_params const& bp);

    void apply_impedance_kick(
            Bunch& bunch, 
            Bunch_params const& bp,
            double wake_factor);

public:

    Impedance(Impedance_options const& ops);
};


#if 0

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
    boost::shared_ptr<MArray1i > bin_partition_sptr;
    

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
    double bunch_z_left;
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
    MArray1i_ref &  get_bin_partition();
    MArray1d_ref const &  get_xmom() const;
    MArray1d_ref const &  get_ymom() const;
    MArray1d_ref const & get_zdensity() const;
    MArray1i_ref const &  get_bin_partition() const;
    
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
    using Collective_operator::apply;
    virtual
    void  apply(Bunch & bunch, double time_step, Step & step, int verbosity, Logger & logger);
    virtual void
    apply(Bunch_train & bunch_train, double time_step, Step & step, int verbosity,
            Train_diagnosticss const& per_operation_train_diagnosticss, Logger & logger);
           
	virtual void
    apply(Bunch_train & bunch_train, double time_step, Step & step, int verbosity,
            Train_diagnosticss const& per_operation_train_diagnosticss, 
            Propagate_actions * propagate_actions_ptr, Stepper & stepper, int step_count,  int turn, 
            Logger & logger);
};
#endif

#endif /* IMPEDANCE_H_ */



