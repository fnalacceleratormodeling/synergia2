#include "impedance.h"

#include "synergia/bunch/core_diagnostics.h"
#include "synergia/foundation/physical_constants.h"

#include "synergia/bunch/period.h"
#include "synergia/utils/kokkos_tools.h"

#include "synergia/utils/simple_timer.h"
//#include "synergia/utils/parallel_utils.h"


typedef Kokkos::TeamPolicy<> team_policy;
typedef typename team_policy::member_type team_member;

namespace
{

    struct alg_zeroer
    {
        karray1d_dev arr;

        KOKKOS_INLINE_FUNCTION
        void operator() (const int i) const
        { arr(i) = 0.0; }
    };

    void zero_karray(karray1d_dev const& arr)
    {
        alg_zeroer alg{arr};
        Kokkos::parallel_for(arr.extent(0), alg);
    }

    struct alg_z_binning
    {
        typedef double value_type[];

        const int value_count;
        const int z_grid;

        ConstParticles p;
        ConstParticleMasks masks;

        const double z_left;
        const double recip_h;

        alg_z_binning(
                ConstParticles const& p,
                ConstParticleMasks const& masks,
                int z_grid, double z_left, double h)
        : value_count(z_grid*3), z_grid(z_grid)
        , p(p), masks(masks)
        , z_left(z_left), recip_h(1.0/h)
        { }

        KOKKOS_INLINE_FUNCTION
        void operator() (const int i, value_type sum) const
        {
            if (masks(i))
            {
                int bin = (p(i, 4) - z_left) * recip_h;

                if (bin < 0 || bin >= z_grid)
                    throw std::runtime_error("impedance z-binning out of range");

                sum[z_grid*0 + bin] += 1;        // zdensity
                sum[z_grid*1 + bin] += p(i, 0);  // xmom
                sum[z_grid*2 + bin] += p(i, 2);  // ymom
            }
        }
    };

    struct alg_z_normalize
    {
        karray1d_dev binning;
        const int z_grid;

        KOKKOS_INLINE_FUNCTION
        void operator() (const int i) const
        {
            double zden = binning(z_grid*0 + i);

            if (zden)
            {
                binning(z_grid*1+i) /= zden;
                binning(z_grid*2+i) /= zden;
            }
        }
    };

    KOKKOS_INLINE_FUNCTION
    int get_zindex_for_wake(double z, double dz, int istart, double zstart)
    { 
        //if  (z< (-istart*istart*dz+zstart)) return -100;
        if ( z >= zstart ) return (static_cast<int>(floor(sqrt((z-zstart)/dz)))) + istart;
        else return (-static_cast<int>(floor(sqrt((zstart-z)/dz)))) +istart;
    }  

    struct alg_z_wake_reduce
    {
        typedef kt::array_type<double, 5> value_type;

        const int i;
        const int z_grid;

        const double cell_size_z;
        const double N_factor;

        const int size_wake;
        const int istart;
        const double zstart;
        const double delta_z;

        karray1d_dev const& zbins;
        karray1d_dev const& wf;

        alg_z_wake_reduce(int i, int z_grid,
                Bunch_params const& bp,
                Wake_field const& wf,
                karray1d_dev const& zbins )
            : i(i), z_grid(z_grid)
            , cell_size_z(bp.cell_size_z), N_factor(bp.N_factor)
            , size_wake(wf.size_wake), istart(wf.istart)
            , zstart(wf.zstart), delta_z(wf.delta_z)
            , zbins(zbins), wf(wf.terms)
        { }

        KOKKOS_INLINE_FUNCTION
        void operator() (const int j, value_type& sum) const
        {
            double* zdensity = &zbins(z_grid*0);
            double* xmom     = &zbins(z_grid*1);
            double* ymom     = &zbins(z_grid*2);

            double* z_coord  = &wf(size_wake*0);
            double* z_wake   = &wf(size_wake*1);
            double* xw_lead  = &wf(size_wake*2);
            double* xw_trail = &wf(size_wake*3);
            double* yw_lead  = &wf(size_wake*4);
            double* yw_trail = &wf(size_wake*5);

            double zji = (j - i) * cell_size_z;
            if (zji < z_coord[0]) return;

            // below it is assumed the wake function is stored using a quadratic grid
            int iz = get_zindex_for_wake(zji, delta_z, istart, zstart);

            if (iz+1 < size_wake)
            {
                double z1 = zji - z_coord[iz];
                double recip_z2 = 1.0 / (z_coord[iz+1] - z_coord[iz]);

                double xwl = xw_lead[iz]  + z1 * (xw_lead[iz+1] - xw_lead[iz])   * recip_z2;
                sum.arr[0] += zdensity[j] * N_factor * xmom[j] * xwl;

                double xwt = xw_trail[iz] + z1 * (xw_trail[iz+1] - xw_trail[iz]) * recip_z2;
                sum.arr[1] += zdensity[j] * N_factor * xwl;

                double ywl = yw_lead[iz]  + z1 * (yw_lead[iz+1] - yw_lead[iz])   * recip_z2;
                sum.arr[2] += zdensity[j] * N_factor * ymom[j] * ywl;

                double ywt = yw_trail[iz] + z1 * (yw_trail[iz+1] - yw_trail[iz]) * recip_z2;
                sum.arr[3] += zdensity[j] * N_factor * ywl;

                double zw = z_wake[iz] + z1 * (z_wake[iz+1] - z_wake[iz]) * recip_z2;
                sum.arr[4] += zdensity[j] * N_factor * zw;
            }
        }
    };


    struct alg_z_wake
    {
        Bunch_params bp;
        Wake_field wf;
        karray1d_dev zbins;
        karray1d_dev wakes;

        alg_z_wake( Bunch_params const& bp,
                Wake_field const& wf,
                karray1d_dev const& zbins,
                karray1d_dev const& wakes )
            : bp(bp), wf(wf), zbins(zbins), wakes(wakes)
        { }

        KOKKOS_INLINE_FUNCTION
        void operator() (const team_member& member) const
        {
            // league size is the z_grid
            int z_grid = member.league_size();

            // li is the index i in zji
            int li = member.league_rank();
            int ti = member.team_rank();

            typedef kt::array_type<double, 5> value_t;
            typedef kt::SumArray<double, 5> array_sum_res_t;

            value_t sum;

            alg_z_wake_reduce z_wake_reduce{
                li, z_grid, bp, wf, zbins
            };

            Kokkos::parallel_reduce(
                    TeamThreadRange(member, z_grid),
                    z_wake_reduce,
                    array_sum_res_t(sum) 
            );

            if (ti == 0)
            {
                wakes(z_grid*0+li) = sum.arr[0];
                wakes(z_grid*1+li) = sum.arr[1];
                wakes(z_grid*2+li) = sum.arr[2];
                wakes(z_grid*3+li) = sum.arr[3];
                wakes(z_grid*4+li) = sum.arr[4];
            }
        }
    };

}

Impedance::Impedance(Impedance_options const& opts)
    : Collective_operator("impedance", 1.0)
    , opts(opts)
    , wake_field(opts.wake_file, opts.wake_type)
{
}



void 
Impedance::apply_impl(Bunch_simulator& sim, 
        double time_step, Logger& logger)
{
    logger << "    Impedance\n";

    scoped_simple_timer timer("imp_total");

    // construct the workspace for a new bunch simulator
    if (bunch_sim_id != sim.id())
    {
        construct_workspaces(sim);
        bunch_sim_id = sim.id();
    }

    // pre-work
    store_bunches_data(sim);

    std::cout << "store bunch data\n";

    // apply to bunches
    for(size_t t=0; t<2; ++t)
    {
        for(size_t b=0; b<sim[t].get_bunch_array_size(); ++b)
        {
            apply_bunch(sim[t][b], time_step, logger);
        }
    }
}
 
void
Impedance::construct_workspaces(Bunch_simulator const& sim)
{
    zbinning = karray1d("zbinning", opts.z_grid*3);
    h_zbinning = Kokkos::create_mirror_view(zbinning);

    wakes = karray1d("wakes", opts.z_grid*5);
    h_wakes = Kokkos::create_mirror_view(wakes);
}

void
Impedance::store_bunches_data(Bunch_simulator const& sim)
{
    auto const& train = sim[0];
    auto num_bunches = train.get_num_bunches();
    auto num_local_bunches = train.get_num_local_bunches();

    // each bunch has 5 properties, x/y/z_mean, real_num, and bucket_idx
    karray1d vbi_buf("vbi_buf", num_bunches * 5);

    for (int i = 0; i < num_local_bunches; ++i)
    {
        auto const& bunch = train[i];
        auto means = Core_diagnostics::calculate_mean(bunch);

        int bucket_idx = bunch.get_bucket_index();
        int bunch_idx = bunch.get_bunch_index();

        if (opts.full_machine && (bucket_idx != bunch_idx))
            throw std::runtime_error(
                    "for full_machine the buckets have to be occupied in order");		     		    
        if (bunch.get_comm().rank() == 0)
        {
            vbi_buf[bunch_idx*5+0] = means[0];
            vbi_buf[bunch_idx*5+1] = means[2];
            vbi_buf[bunch_idx*5+2] = means[4];
            vbi_buf[bunch_idx*5+3] = bunch.get_real_num();
            vbi_buf[bunch_idx*5+4] = bucket_idx;
        }
    }

    int error = MPI_Allreduce( MPI_IN_PLACE,
                               (void*)vbi_buf.data(),
                               num_bunches * 5,
                               MPI_DOUBLE,
                               MPI_SUM,
                               train.get_comm() );

    if (error != MPI_SUCCESS) 
        throw std::runtime_error(
                "Impedance::store_bunches_data: MPI error in MPI_Allreduce");
		
    std::vector<Bunch_properties> vbi(num_bunches);
    for(int i=0; i<num_bunches; ++i)
    {
        vbi[i].x_mean = vbi_buf[i*5+0];
        vbi[i].y_mean = vbi_buf[i*5+1];
        vbi[i].z_mean = vbi_buf[i*5+2];
        vbi[i].realnum = vbi_buf[i*5+3];
        vbi[i].bucket_index = (int)vbi_buf[i*5+4];
    }

    stored_vbunches.push_front(vbi);

    if (stored_vbunches.size() > nstored_turns) 
        stored_vbunches.pop_back();
}

void
Impedance::apply_bunch(Bunch& bunch, 
        double time_step, Logger& logger)
{
    std::cout << "apply bunch\n";

    auto bp = calculate_moments_and_partitions(bunch);   
    std::cout << "cal moments\n";
    //t = simple_timer_show(t, "impedance_apply:  calculate_moments_and_partitions ");

    auto means = Core_diagnostics::calculate_mean(bunch);
    bp.z_mean = means[4];
    bp.N_factor = bunch.get_real_num() / bunch.get_total_num();
    std::cout << "cal z mean\n";
    //t = simple_timer_show(t, "impedance_apply: calculate_z_mean ");

    int bunch_bucket = bunch.get_bucket_index(); 
    calculate_kicks(bunch, bp);
    std::cout << "cal kicks\n";
    //t = simple_timer_show(t, "impedance_apply: calculate_kicks ");

#if 0
    double gamma = bunch.get_reference_particle().get_gamma();
    double beta = bunch.get_reference_particle().get_beta();
    double w_f = get_wake_factor()*time_step/(gamma*beta);  
    
    apply_impedance_kick(bunch,  w_f);
    t = simple_timer_show(t, "impedance apply:apply_impedance_kick ");
#endif
}


Bunch_params
Impedance::calculate_moments_and_partitions(Bunch const& bunch)
{
    // output cell_size_z, xmom, ymom, zdensity
    Bunch_params bp;
  
    auto bunchmin = Core_diagnostics::calculate_min(bunch);
    bp.z_left = bunchmin[2];   

    auto bunchmax = Core_diagnostics::calculate_max(bunch);
    double z_length = bunchmax[2] - bp.z_left;
    
    if (z_length <= 1.e-14) 
        throw std::runtime_error("z_length too small ");

    // 1e-14 is to make sure the max-z particle falls in the last bin
    bp.cell_size_z = z_length / double(opts.z_grid) + 1e-14;

    std::cout << "cellsize = " << bp.cell_size_z << "\n";

    // double h = z_length/(opts.z_grid-1.0); // AM why have I done that???
    double h = bp.cell_size_z;

    // TODO: need to verify if zeroing h_zbinning 
    // is needed for GPU version
    // ...

    // get binning results
    auto parts = bunch.get_local_particles();
    auto masks = bunch.get_local_particle_masks();
    alg_z_binning alg(parts, masks, opts.z_grid, bp.z_left, h);

    Kokkos::parallel_reduce("z_binning", 
            bunch.size(), alg, h_zbinning.data());

    std::cout << "zbinning\n";

    // MPI reduction to get global z-binning results
    int error = MPI_Allreduce(MPI_IN_PLACE, h_zbinning.data(), 
            opts.z_grid*3, MPI_DOUBLE, MPI_SUM, bunch.get_comm());

    if (error != MPI_SUCCESS)
        throw std::runtime_error("MPI error in Impedance reduce z_binning");

    // copy back to device
    Kokkos::deep_copy(zbinning, h_zbinning);

    // normalize
    alg_z_normalize alg2{zbinning, opts.z_grid};
    Kokkos::parallel_for(opts.z_grid, alg2);

    std::cout << "normalize\n";

    return bp;
}   


void Impedance::calculate_kicks(Bunch const& bunch, Bunch_params const& bp)
{
    std::cout << "cal kicks\n";

    int registered_turns = stored_vbunches.size();
    int numbunches;
    int num_trains = 0;

    if (registered_turns == 0)
        throw std::runtime_error(
                "registered_turns size cannot be zero, "
                "probably you propagate a bunch instead of a bunch_train");

    numbunches = (*stored_vbunches.begin()).size();

    if (opts.full_machine) 
    { 
        num_trains = int(opts.num_buckets / numbunches);  
        /// num_trains is relevant only when the full machine option is considered 
        /// a full machine consideres a num_train of bunches repeats with modulation wave wn[] 
        /// all buckets are full, but only numbunches bunches properties are stored 
        /// exemple: full_machine, all bunches identical, no wave across the machine: 
        /// num_trains=num_buckets, wn=[0,0,0], it's a one bunch simulation  
        /// example: full_machine, two bunch simulation, num_trains= num_buckets/2  

        if (std::abs(opts.num_buckets / float(numbunches) - num_trains) > 1e-8) 
            throw std::runtime_error( 
                    "full machine assumes repetitive numer of trains: "
                    "num_buckets should be divisible to numbunches");

        if (opts.wn[0] < 0 || opts.wn[0] >= num_trains 
                || opts.wn[1] < 0 || opts.wn[1] >= num_trains 
                || opts.wn[2] < 0 || opts.wn[2] >= num_trains) 
            throw std::runtime_error(
                    "full machine wave number cannot be smaller than zero "
                    "or larger than num_trains-1"); 
    }

    std::cout << "here\n";


    using Kokkos::TeamPolicy;
    using Kokkos::TeamThreadRange;

    // zbinning: zdensity, xmom, ymom
    // wakes: xw_lead, xw_trail, yw_lead, yw_trail, zwake
    alg_z_wake ft_z_wake{bp, wake_field, zbinning, wakes};
    Kokkos::parallel_for(TeamPolicy<>(opts.z_grid, 32), ft_z_wake);

    MPI_Abort(MPI_COMM_WORLD, 333);

   
#if 0
   
   // for (int i = offsets[rank]; i <offsets[rank]+counts[rank] ; ++i){
    for (int i = 0; i <counts[rank] ; ++i)
    {  
        int real_i=i+offsets[rank];

        xwake_leading_local[i]=0.; 
        xwake_trailing_local[i]=0. ;
        ywake_leading_local[i] =0.; 
        ywake_trailing_local[i]=0. ;
        zwake0_local[i] =0.;	

        // in-bunch impedance 
        for (int j = 0; j < z_grid; ++j)
        {
            double zji=(j-real_i)*cell_size_z;
    
            if (zji >= z_coord[0]) 
            {
                // below it is assumed the wake function is stored using a quadratic grid
                int iz = get_zindex_for_wake(zji, delta_z, istart, zstart);
                double xwl(0.), xwt(0.), ywl(0.), ywt(0.), zw(0.);

                if (iz+1 < zpoints) 
                { 
                    xwl = xw_lead[iz] + (zji-z_coord[iz]) * (xw_lead[iz+1]-xw_lead[iz]) / 
                        (z_coord[iz+1]-z_coord[iz]); 

                    xwt = xw_trail[iz] + (zji-z_coord[iz]) * (xw_trail[iz+1]-xw_trail[iz]) /
                        (z_coord[iz+1]-z_coord[iz]);

                    ywl=yw_lead[iz]+(zji-z_coord[iz])*(yw_lead[iz+1]-yw_lead[iz])/
                        (z_coord[iz+1]-z_coord[iz]); 

                    ywt=yw_trail[iz]+(zji-z_coord[iz])*(yw_trail[iz+1]-yw_trail[iz])/
                        (z_coord[iz+1]-z_coord[iz]);                  

                    zw=z_wake[iz]+(zji-z_coord[iz])*(z_wake[iz+1]-z_wake[iz])/
                        (z_coord[iz+1]-z_coord[iz]); 
                }

                xwake_leading_local[i]  +=zdensity[j]*N_factor*xmom[j]*xwl; 
                xwake_trailing_local[i]  += zdensity[j]*N_factor*xwt;
                ywake_leading_local[i]  += zdensity[j]*N_factor*ymom[j]*ywl; 
                ywake_trailing_local[i]  += zdensity[j]*N_factor*ywt;
                zwake0_local[i] += zdensity[j]*N_factor*zw;				
            }          
         } 

        // stored_vbunches.begin() stores the bunches info at 
        std::list< std::vector<Bunch_properties> >::const_iterator it=stored_vbunches.begin(); 

        // at the moment
        /// bucket 0 is in front of bucket 1, which is in front of bucket 2, etc...
        int mean_bin = static_cast<int>((bunch_z_mean-bunch_z_left)/cell_size_z);

        if ((mean_bin<0) || (mean_bin>=z_grid)){
            throw std::runtime_error(
                    "impedance: the index bin of beam min cannot be <0 or >z_grid, something is wrong  ");
        }

        double z_to_zmean = (mean_bin-real_i)*cell_size_z; 

        for (int ibunch= 0; ibunch<numbunches; ++ibunch)
        {            
            //  double xwl(0.), xwt(0.), ywl(0.), ywt(0.), zw(0.);
            int ibucket=(*it)[ibunch].bucket_index;

            if(ibucket<bunch_bucket) {///  same turn, the leading buckets effect    
                double  zji=z_to_zmean+bunch_spacing*(bunch_bucket-ibucket) +((*it)[ibunch].z_mean-bunch_z_mean);
            //  int iz=static_cast<int>(floor(sqrt((zji-z_coord[0])/(z_coord[1]-z_coord[0]))));  
                int iz=get_zindex_for_wake(zji, delta_z, istart, zstart);
                if ((iz+1 < zpoints) && (iz>0)) {
                    double xwl=xw_lead[iz]+(zji-z_coord[iz])*(xw_lead[iz+1]-xw_lead[iz])/(z_coord[iz+1]-z_coord[iz]);
                    double xwt=xw_trail[iz]+(zji-z_coord[iz])*(xw_trail[iz+1]-xw_trail[iz])/(z_coord[iz+1]-z_coord[iz]);
                    double ywl=yw_lead[iz]+(zji-z_coord[iz])*(yw_lead[iz+1]-yw_lead[iz])/(z_coord[iz+1]-z_coord[iz]);
                    double ywt=yw_trail[iz]+(zji-z_coord[iz])*(yw_trail[iz+1]-yw_trail[iz])/(z_coord[iz+1]-z_coord[iz]);                  
                    double zw=z_wake[iz]+(zji-z_coord[iz])*(z_wake[iz+1]-z_wake[iz])/(z_coord[iz+1]-z_coord[iz]); 

                    xwake_leading_local[i]  +=(*it)[ibunch].realnum*(*it)[ibunch].x_mean*xwl; 
                    xwake_trailing_local[i]  += (*it)[ibunch].realnum*xwt;
                    ywake_leading_local[i]  +=  (*it)[ibunch].realnum*(*it)[ibunch].y_mean*ywl;
                    ywake_trailing_local[i]  +=(*it)[ibunch].realnum*ywt;
                    zwake0_local[i] += (*it)[ibunch].realnum*zw;
                }
            }
        } // ibunch loop

        if (full_machine) { // it assumes that all the other trains are in front 
            for (int itrain= 1; itrain<num_trains; ++itrain){  //for one bunch simulation num_trains = num_buckets, all trains identical except displacement
                double wnx=cos(2.*mconstants::pi*wn[0]*itrain/double(num_trains));
                double wny=cos(2.*mconstants::pi*wn[1]*itrain/double(num_trains));
                double wnz=cos(2.*mconstants::pi*wn[2]*itrain/double(num_trains));
                for (int ibunch= 0; ibunch<numbunches; ++ibunch){		  
                    double  zji=z_to_zmean+bunch_spacing*numbunches*itrain+bunch_spacing*(bunch_bucket-ibunch)+
                    ((*it)[ibunch].z_mean*wnz-bunch_z_mean);
                    //  int iz=static_cast<int>(floor(sqrt((zji-z_coord[0])/(z_coord[1]-z_coord[0])))); 
                    int iz=get_zindex_for_wake(zji, delta_z, istart, zstart);
                    if ((iz+1 < zpoints) && (iz>0)) {
                        double xwl=xw_lead[iz]+(zji-z_coord[iz])*(xw_lead[iz+1]-xw_lead[iz])/(z_coord[iz+1]-z_coord[iz]);
                        double xwt=xw_trail[iz]+(zji-z_coord[iz])*(xw_trail[iz+1]-xw_trail[iz])/(z_coord[iz+1]-z_coord[iz]);
                        double ywl=yw_lead[iz]+(zji-z_coord[iz])*(yw_lead[iz+1]-yw_lead[iz])/(z_coord[iz+1]-z_coord[iz]);
                        double ywt=yw_trail[iz]+(zji-z_coord[iz])*(yw_trail[iz+1]-yw_trail[iz])/(z_coord[iz+1]-z_coord[iz]);                  
                        double zw=z_wake[iz]+(zji-z_coord[iz])*(z_wake[iz+1]-z_wake[iz])/(z_coord[iz+1]-z_coord[iz]);    

                        xwake_leading_local[i]  +=(*it)[ibunch].realnum*(*it)[ibunch].x_mean*wnx*xwl; 
                        xwake_trailing_local[i]  += (*it)[ibunch].realnum*xwt;
                        ywake_leading_local[i]  +=  (*it)[ibunch].realnum*(*it)[ibunch].y_mean*wny*ywl;
                        ywake_trailing_local[i]  +=(*it)[ibunch].realnum*ywt;
                        zwake0_local[i] += (*it)[ibunch].realnum*zw;
                    }
                }
            }
        }  // full_machine
        
        if (registered_turns>1) {
            ++it; ///previous turn, following buckets effect
            for (int ibunch= 0; ibunch<numbunches; ++ibunch){
                int ibucket=(*it)[ibunch].bucket_index;
                if(ibucket>=bunch_bucket) {///  following buckets effect    
                    double  zji=z_to_zmean+bunch_spacing*(bunch_bucket-ibucket)+orbit_length +((*it)[ibunch].z_mean-bunch_z_mean);
                    //  int iz=static_cast<int>(floor(sqrt((zji-z_coord[0])/(z_coord[1]-z_coord[0])))); 
                    int iz=get_zindex_for_wake(zji, delta_z, istart, zstart);
                    if ((iz+1 < zpoints) && (iz>0)) {
                        double xwl=xw_lead[iz]+(zji-z_coord[iz])*(xw_lead[iz+1]-xw_lead[iz])/(z_coord[iz+1]-z_coord[iz]);
                        double xwt=xw_trail[iz]+(zji-z_coord[iz])*(xw_trail[iz+1]-xw_trail[iz])/(z_coord[iz+1]-z_coord[iz]);
                        double ywl=yw_lead[iz]+(zji-z_coord[iz])*(yw_lead[iz+1]-yw_lead[iz])/(z_coord[iz+1]-z_coord[iz]);
                        double ywt=yw_trail[iz]+(zji-z_coord[iz])*(yw_trail[iz+1]-yw_trail[iz])/(z_coord[iz+1]-z_coord[iz]);                  
                        double zw=z_wake[iz]+(zji-z_coord[iz])*(z_wake[iz+1]-z_wake[iz])/(z_coord[iz+1]-z_coord[iz]);    
                        
                        xwake_leading_local[i]  +=(*it)[ibunch].realnum*(*it)[ibunch].x_mean*xwl; 
                        xwake_trailing_local[i]  += (*it)[ibunch].realnum*xwt;
                        ywake_leading_local[i]  +=  (*it)[ibunch].realnum*(*it)[ibunch].y_mean*ywl;
                        ywake_trailing_local[i]  +=(*it)[ibunch].realnum*ywt;
                        zwake0_local[i] += (*it)[ibunch].realnum*zw;                   
                    }
                }
            } // ibunch loop	 
        }//registered_turns>1         
    } // i loop
    t1 = simple_timer_show(t1, "impedance_calculate_kicks:  i loop ");
    /// it is not necessary to have a loop over i at larger distances, since the effect is negligible
    /// might not be  true for coasting beams!
    if (registered_turns>1) {
        double xwake_l=0.;
        double xwake_t=0.;
        double ywake_l=0.;
        double ywake_t=0.;
        double zwake_0=0.; 
              
        
        std::list< std::vector<Bunch_properties> >::const_iterator it;
        std::list< std::vector<Bunch_properties> >::const_iterator jt=stored_vbunches.begin();
        ++jt;       
        int iturn;
        for (it=jt,  iturn=1; it !=stored_vbunches.end(); ++it, ++iturn){
            for (int ibunch= 0; ibunch<numbunches; ++ibunch){
                int ibucket=(*it)[ibunch].bucket_index;                   
                if(((ibucket<bunch_bucket) && (iturn==1))///  finishing the previous turn, for the buckets ahead
                ||  (iturn>1))  /// previous turns effects                
                {
                    double  zji=bunch_spacing*(bunch_bucket-ibucket)+orbit_length*iturn +((*it)[ibunch].z_mean-bunch_z_mean);
                    //  int iz=static_cast<int>(floor(sqrt((zji-z_coord[0])/(z_coord[1]-z_coord[0]))));  
                    int iz=get_zindex_for_wake(zji, delta_z, istart, zstart);
                    if ((iz+1 < zpoints) && (iz>0)) {
                        double  xwl=xw_lead[iz]+(zji-z_coord[iz])*(xw_lead[iz+1]-xw_lead[iz])/(z_coord[iz+1]-z_coord[iz]);
                        double  xwt=xw_trail[iz]+(zji-z_coord[iz])*(xw_trail[iz+1]-xw_trail[iz])/(z_coord[iz+1]-z_coord[iz]);
                        double  ywl=yw_lead[iz]+(zji-z_coord[iz])*(yw_lead[iz+1]-yw_lead[iz])/(z_coord[iz+1]-z_coord[iz]);
                        double  ywt=yw_trail[iz]+(zji-z_coord[iz])*(yw_trail[iz+1]-yw_trail[iz])/(z_coord[iz+1]-z_coord[iz]);                  
                        double  zw=z_wake[iz]+(zji-z_coord[iz])*(z_wake[iz+1]-z_wake[iz])/(z_coord[iz+1]-z_coord[iz]);

                        xwake_l +=  (*it)[ibunch].realnum*(*it)[ibunch].x_mean*xwl;
                        xwake_t += (*it)[ibunch].realnum*xwt;
                        ywake_l +=  (*it)[ibunch].realnum*(*it)[ibunch].y_mean*ywl;
                        ywake_t += (*it)[ibunch].realnum*ywt; 
                        zwake_0 += (*it)[ibunch].realnum*zw;
                    }
                } 
            }

            if (full_machine) {
                for (int itrain= 1; itrain<num_trains; ++itrain){
                    double wnx=cos(2.*mconstants::pi*wn[0]*itrain/double(num_trains));
                    double wny=cos(2.*mconstants::pi*wn[1]*itrain/double(num_trains));
                    double wnz=cos(2.*mconstants::pi*wn[2]*itrain/double(num_trains));
                    for (int ibunch= 0; ibunch<numbunches; ++ibunch){
                        double  zji=orbit_length*iturn+bunch_spacing*numbunches*itrain+bunch_spacing*(bunch_bucket-ibunch)+
                        ((*it)[ibunch].z_mean*wnz-bunch_z_mean);
                        //    int iz=static_cast<int>(floor(sqrt((zji-z_coord[0])/(z_coord[1]-z_coord[0]))));
                        int iz=get_zindex_for_wake(zji, delta_z, istart, zstart);
                        if ((iz+1 < zpoints) && (iz>0)) {
                            double  xwl=xw_lead[iz]+(zji-z_coord[iz])*(xw_lead[iz+1]-xw_lead[iz])/(z_coord[iz+1]-z_coord[iz]);
                            double  xwt=xw_trail[iz]+(zji-z_coord[iz])*(xw_trail[iz+1]-xw_trail[iz])/(z_coord[iz+1]-z_coord[iz]);
                            double  ywl=yw_lead[iz]+(zji-z_coord[iz])*(yw_lead[iz+1]-yw_lead[iz])/(z_coord[iz+1]-z_coord[iz]);
                            double  ywt=yw_trail[iz]+(zji-z_coord[iz])*(yw_trail[iz+1]-yw_trail[iz])/(z_coord[iz+1]-z_coord[iz]);                  
                            double  zw=z_wake[iz]+(zji-z_coord[iz])*(z_wake[iz+1]-z_wake[iz])/(z_coord[iz+1]-z_coord[iz]);

                            xwake_l +=  (*it)[ibunch].realnum*(*it)[ibunch].x_mean*wnx*xwl;
                            xwake_t += (*it)[ibunch].realnum*xwt;
                            ywake_l +=  (*it)[ibunch].realnum*(*it)[ibunch].y_mean*wny*ywl;
                            ywake_t += (*it)[ibunch].realnum*ywt; 
                            zwake_0 += (*it)[ibunch].realnum*zw;
                        }		      
                    }//ibunch
                } ////itrain
            }  // full_machine
      } //iturn
      for (int i = 0; i <counts[rank] ; ++i){
            xwake_leading_local[i] +=  xwake_l;
            xwake_trailing_local[i] += xwake_t;
            ywake_leading_local[i] +=  ywake_l; 
            ywake_trailing_local[i] += ywake_t ;
            zwake0_local[i] += zwake_0;    
      }         
   }//registred_turns>1
    t1 = simple_timer_show(t1, "impedance_calculate_kicks: registred_turns loop   ");
   
     int error = MPI_Allgatherv(reinterpret_cast<void*>(xwake_leading_local.origin()),
		  counts[rank], MPI_DOUBLE,
		  reinterpret_cast<void*>(xwake_leading.origin()),
					  &counts[0], &offsets[0], MPI_DOUBLE, comm_sptr->get());
      if (error != MPI_SUCCESS) {
	  throw std::runtime_error(
	      "MPI error in Impedance_calculate_kicks:MPI_Allgatherv xwake_leading"); 
      } 
 

     error = MPI_Allgatherv(reinterpret_cast<void*>(xwake_trailing_local.origin()),
		  counts[rank], MPI_DOUBLE,
		  reinterpret_cast<void*>(xwake_trailing.origin()),
					  &counts[0], &offsets[0], MPI_DOUBLE, comm_sptr->get());
    if (error != MPI_SUCCESS) {
	  throw std::runtime_error(
	      "MPI error in Impedance_calculate_kicks:MPI_Allgatherv xwake_trailing"); 
     } 
   
    error = MPI_Allgatherv(reinterpret_cast<void*>(ywake_leading_local.origin()),
		  counts[rank], MPI_DOUBLE,
		  reinterpret_cast<void*>(ywake_leading.origin()),
					  &counts[0], &offsets[0], MPI_DOUBLE, comm_sptr->get());
    if (error != MPI_SUCCESS) {
	  throw std::runtime_error(
	      "MPI error in Impedance_calculate_kicks:MPI_Allgatherv ywake_leading"); 
     } 
 

    error = MPI_Allgatherv(reinterpret_cast<void*>(ywake_trailing_local.origin()),
		  counts[rank], MPI_DOUBLE,
		  reinterpret_cast<void*>(ywake_trailing.origin()),
					  &counts[0], &offsets[0], MPI_DOUBLE, comm_sptr->get());
    if (error != MPI_SUCCESS) {
	  throw std::runtime_error(
	      "MPI error in Impedance_calculate_kicks:MPI_Allgatherv ywake_trailing"); 
     } 
   
    error = MPI_Allgatherv(reinterpret_cast<void*>(zwake0_local.origin()),
		  counts[rank], MPI_DOUBLE,
		  reinterpret_cast<void*>(zwake0.origin()),
					  &counts[0], &offsets[0], MPI_DOUBLE, comm_sptr->get());
    if (error != MPI_SUCCESS) {
	  throw std::runtime_error(
	      "MPI error in Impedance_calculate_kicks:MPI_Allgatherv zwake0"); 
     }    
    t1 = simple_timer_show(t1, "impedance_calculate_kicks: mpi_allgather   ");
#endif
}

void 
Impedance::apply_impedance_kick(Bunch & bunch, double wake_factor)
{
#if 0
  
 MArray1i_ref const bin_partition(get_bin_partition());
 MArray1d_ref const  xwake_leading(get_xwake_leading());
 MArray1d_ref const xwake_trailing(get_xwake_trailing());
 MArray1d_ref const  ywake_leading(get_ywake_leading());
 MArray1d_ref const  ywake_trailing(get_ywake_trailing()); 
 MArray1d_ref const zwake0(get_zwake0());
 
 
 for (int part = 0; part < bunch.get_local_num(); ++part) {
        double xkick=0., ykick=0., zkick=0.;
        int bin=bin_partition[part];  // bin_partition(n) is the bin where you find the particle n 
//            if ((bin>=z_grid) || (bin<0))  { std::cout<<"bin="<<bin<<"z_grid="<<z_grid<<std::endl;
//         throw
//         std::runtime_error("something wrong with bining");} 

        xkick=xwake_leading[bin]+xwake_trailing[bin]*bunch.get_local_particles()[part][0];	
        ykick=ywake_leading[bin]+ywake_trailing[bin]*bunch.get_local_particles()[part][2];
        zkick = zwake0[bin];
    
       
        bunch.get_local_particles()[part][1] += wake_factor*xkick;   
        bunch.get_local_particles()[part][3]  += wake_factor*ykick;
        bunch.get_local_particles()[part][5]  += wake_factor*zkick;
      
      
    }
#endif
}





#if 0

Impedance::Impedance(
        std::string const & wake_file, 
        std::string const & wake_type, 
        int const  & zgrid, 
        double const & orbit_length, 
        double const & bunchsp, 
        int const nstored_turns, 
        bool full_machine, 
        std::vector<int > wn)
    : Collective_operator("impedance")
    , z_grid(zgrid)
    , nstored_turns(nstored_turns)
    , orbit_length(orbit_length)
    , bunch_spacing(bunchsp)
    , full_machine(full_machine)
    ,wn(wn)
{

	try {
		if (std::abs(orbit_length/bunchsp - int(orbit_length/bunchsp + 0.5))
				> 1e-8)
			throw std::runtime_error(
					"orbit length should divide exacty to bunch_spacing ");
	} catch (std::exception const& e) {
		std::cout << e.what() << " but the division is "
				<< std::setprecision(17) << orbit_length/bunchsp
				<< " with remainder " << std::setprecision(17)
				<< std::abs(orbit_length/bunchsp - int(orbit_length/bunchsp + 0.5))
				<< std::endl;
		MPI_Abort(MPI_COMM_WORLD, 137);
	}

	this->wake_field_sptr = Wake_field_sptr(
			new Wake_field(wake_file, wake_type));
	this->num_buckets = int(orbit_length/bunchsp + 0.5);
	construct();
}  

Impedance::Impedance(
        std::string const & wake_file, 
        std::string const & wake_type, 
        int const  & zgrid, 
        double const & orbit_length, 
        int const& num_buckets, 
        int const nstored_turns, 
        bool full_machine, 
        std::vector<int > wn):
Collective_operator("impedance"), z_grid(zgrid), nstored_turns(nstored_turns), 
		num_buckets(num_buckets), orbit_length(orbit_length),
              full_machine(full_machine), wn(wn)
{       
   this->wake_field_sptr=Wake_field_sptr(new Wake_field(wake_file, wake_type)); 
   this->bunch_spacing=orbit_length/num_buckets;
   construct();        
}  

void
Impedance::construct()
{
  try{

    wake_factor=-4.*mconstants::pi*pconstants::rp;   ///N.B. the wakefiled file reads W/(Z_0*L), Z_0=1/(epsilon_0*c)

    if (wn.size() !=3) {
	      wn.resize(3);
	      wn[0]=0;
	      wn[1]=0;
	      wn[2]=0;
       }

    stored_vbunches=std::list< std::vector<Bunch_properties> >();

   // xmom_sptr= boost::shared_ptr<Raw_MArray1d >(new Raw_MArray1d(boost::extents[z_grid]));
    xmom_sptr= boost::shared_ptr<MArray1d >(new MArray1d(boost::extents[z_grid]));
    ymom_sptr= boost::shared_ptr<MArray1d >(new MArray1d(boost::extents[z_grid]));
    zdensity_sptr= boost::shared_ptr<MArray1d >(new MArray1d(boost::extents[z_grid]));
    xwake_leading_sptr= boost::shared_ptr<MArray1d >(new MArray1d(boost::extents[z_grid]));
    xwake_trailing_sptr= boost::shared_ptr<MArray1d >(new MArray1d(boost::extents[z_grid]));
    ywake_leading_sptr= boost::shared_ptr<MArray1d >(new MArray1d(boost::extents[z_grid]));
    ywake_trailing_sptr= boost::shared_ptr<MArray1d >(new MArray1d(boost::extents[z_grid]));
    zwake0_sptr= boost::shared_ptr<MArray1d >(new MArray1d(boost::extents[z_grid])); 
    
    
  }
  catch (std::exception const& e){
        std::cout<<e.what()<<std::endl;   
        MPI_Abort(MPI_COMM_WORLD, 777);
  }  
} 



Impedance::Impedance(Impedance const& impedance)
{

   this->wake_field_sptr=impedance.wake_field_sptr;
   this->z_grid=impedance.z_grid;
   this->nstored_turns=impedance.nstored_turns;
   this->orbit_length=impedance.orbit_length;
   this->num_buckets=impedance.num_buckets;
   this->wake_factor=impedance.wake_factor;
   this->bunch_spacing=impedance.bunch_spacing;
   this->full_machine=impedance.full_machine;
   this->wn=impedance.wn;
   
 // the following data are not copied

//    this->stored_vbunches=std::list< std::vector<Bunch_properties> >(); 
//    this->xmom_sptr= boost::shared_ptr<MArray1d >(new MArray1d(boost::extents[z_grid]));
//    this->ymom_sptr= boost::shared_ptr<MArray1d >(new MArray1d(boost::extents[z_grid]));
//    this->zdensity_sptr= boost::shared_ptr<MArray1d >(new MArray1d(boost::extents[z_grid]));
//    this->xwake_leading_sptr= boost::shared_ptr<MArray1d >(new MArray1d(boost::extents[z_grid]));
//    this->xwake_trailing_sptr= boost::shared_ptr<MArray1d >(new MArray1d(boost::extents[z_grid]));
//    this->ywake_leading_sptr= boost::shared_ptr<MArray1d >(new MArray1d(boost::extents[z_grid]));
//    this->ywake_trailing_sptr= boost::shared_ptr<MArray1d >(new MArray1d(boost::extents[z_grid]));
//    this->zwake0_sptr= boost::shared_ptr<MArray1d >(new MArray1d(boost::extents[z_grid]));
  
   this->xmom_sptr=impedance.xmom_sptr; 
   this->ymom_sptr= impedance.ymom_sptr;
   this->zdensity_sptr=impedance. zdensity_sptr;
   this->xwake_leading_sptr=impedance.xwake_leading_sptr; 
   this->xwake_trailing_sptr=impedance.xwake_trailing_sptr;
   this->ywake_leading_sptr=impedance.ywake_leading_sptr;
   this->ywake_trailing_sptr=impedance.ywake_trailing_sptr;
   this->zwake0_sptr=impedance.zwake0_sptr;
   
}  


Impedance *
Impedance::clone()
{
    return new Impedance(*this);
}




void 
Impedance::set_z_grid(int const  & zgrid)
{
  this->z_grid=zgrid;
}   



Wake_field_sptr 
Impedance::get_wake_field_sptr() const
{
  return wake_field_sptr;
}  


int Impedance::get_z_grid() const { return z_grid;} 
double Impedance::get_orbit_length() const{ return orbit_length;} 
double Impedance::get_wake_factor() const { return wake_factor;} 
double Impedance::get_bunch_spacing() const { return bunch_spacing;}


MArray1d_ref &  Impedance::get_xmom() {return *xmom_sptr;}
MArray1d_ref &  Impedance::get_ymom() {return *ymom_sptr;}
MArray1d_ref &  Impedance::get_zdensity() {return *zdensity_sptr;}
MArray1i_ref & Impedance::get_bin_partition() {return *bin_partition_sptr;}


MArray1d_ref const &  Impedance::get_xmom() const {return *xmom_sptr;}
MArray1d_ref const &  Impedance::get_ymom() const {return *ymom_sptr;}
MArray1d_ref const &  Impedance::get_zdensity() const {return *zdensity_sptr;}
MArray1i_ref const & Impedance::get_bin_partition() const {return *bin_partition_sptr;}

MArray1d_ref &  Impedance::get_xwake_leading() {return *xwake_leading_sptr;}
MArray1d_ref &  Impedance::get_xwake_trailing() {return *xwake_trailing_sptr;}
MArray1d_ref &  Impedance::get_ywake_leading() {return *ywake_leading_sptr;}
MArray1d_ref &  Impedance::get_ywake_trailing() {return *ywake_trailing_sptr;}
MArray1d_ref &  Impedance::get_zwake0() {return *zwake0_sptr;}

MArray1d_ref const &  Impedance::get_xwake_leading() const {return *xwake_leading_sptr;}
MArray1d_ref const &  Impedance::get_xwake_trailing() const {return *xwake_trailing_sptr;}
MArray1d_ref const &  Impedance::get_ywake_leading() const {return *ywake_leading_sptr;}
MArray1d_ref const &  Impedance::get_ywake_trailing() const {return *ywake_trailing_sptr;}
MArray1d_ref const &  Impedance::get_zwake0() const {return *zwake0_sptr;}

bool Impedance::is_full_machine() const { return full_machine;}
int Impedance::get_nstored_turns() const { return nstored_turns;}
std::list< std::vector<Bunch_properties> > &
Impedance::get_stored_vbunches() {return stored_vbunches;}

int Impedance::get_num_buckets() const { return num_buckets;}
std::vector<int >  Impedance::get_train_wave() const { return wn;}



void
Impedance::calculate_moments_and_partitions(Bunch & bunch)
{
/// output cell_size_z, xmom, ymom, zdensity
  
    int rank(bunch.get_comm().get_rank());
    
  
    MArray1d bunchmin(Core_diagnostics::calculate_min(bunch));
    bunch_z_left=bunchmin[2];   
    double z_left=bunch_z_left;
    MArray1d bunchmax(Core_diagnostics::calculate_max(bunch)); 
    double z_length=bunchmax[2]-z_left;       
    cell_size_z= z_length/double(z_grid);
    
   // double h = z_length/(z_grid-1.0); // AM why have I done that???
    double h = cell_size_z;
    if (z_length<= 1.e-14 )   throw
                 std::runtime_error("z_length too small ");


    MArray1d_ref xmom(get_xmom());
    MArray1d_ref ymom(get_ymom());
    MArray1d_ref zdensity(get_zdensity());

    int lnum_part=bunch.get_local_num();
    bin_partition_sptr= boost::shared_ptr<MArray1i >(new MArray1i(boost::extents[lnum_part])); 
    MArray1i_ref bin_partition(get_bin_partition());
    
    MArray1d  local_zdensity(boost::extents[z_grid]);
    MArray1d  local_xmom(boost::extents[z_grid]);
    MArray1d  local_ymom(boost::extents[z_grid]);
     
      
    for (int i=0; i<z_grid;  ++i){
        local_zdensity[i]=0.0;
        local_xmom[i]=0.0;
        local_ymom[i]=0.0;
    }



     for (int part = 0;  part < bunch.get_local_num(); ++part) {
         int bin = static_cast<int>((bunch.get_local_particles()[part][4]-z_left)/h);
         if ((bin < z_grid) && (bin >= 0)) {
             local_zdensity[bin] += 1;
             local_xmom[bin] += bunch.get_local_particles()[part][0];
             local_ymom[bin] += bunch.get_local_particles()[part][2];
             bin_partition[part]=bin; //bin_partition(n) is the bin where you find the particle n
         }
         else if ((bin==z_grid) && std::abs(bunch.get_local_particles()[part][4]-z_length-z_left)<z_length*1.e-14) { 
             local_zdensity[bin-1] += 1; // put the edge particle in the last bin=z_grid-1
             bin_partition[part]=z_grid-1;
         } 
         else
         {   std::cout << "  z_left  "<<z_left<<"  rank= "<<rank<<std::endl;
          std::cout<<"bunch.get_local_particles()[part][4]="  <<bunch.get_local_particles()[part][4]<<"  rank= "<<rank<<std::endl; 
           std::cout<<"bunch.get_local_particles()[part]0,1,2,3,4,5="  <<bunch.get_local_particles()[part][0]<<
           "  "<<bunch.get_local_particles()[part][1]<<
           "  "<<bunch.get_local_particles()[part][2]<<
           "  "<<bunch.get_local_particles()[part][3]<<
           "  "<<bunch.get_local_particles()[part][4]<<
           "  "<<bunch.get_local_particles()[part][5]<<std::endl;
           
         std::cout<< " particle's id ="<<part<<std::endl; 
         std::cout << "  z_length  "<<z_length<<"  rank= "<<rank<<std::endl;
         std::cout << "(mbs.local_particles(4,n)-z_left)= "<<(bunch.get_local_particles()[part][4]-z_left)<<"  rank= "<<rank<<std::endl;
         std::cout << "bin: " << bin<<"  z_grid="<<z_grid<< "  h=" << h <<"  rank= "<<rank<<std::endl;                
         std::cout << "bunch.get_local_particles()[part][4]-z_length-z_left= "<<std::abs(bunch.get_local_particles()[part][4]-z_length-z_left)<<"  rank= "<<rank<<std::endl;
         throw
                 std::runtime_error("particles out of range in calculate_moments_and_partitions ");
         }

    
    }


    int error = MPI_Allreduce(reinterpret_cast<void*>(local_zdensity.origin()), 
			      reinterpret_cast<void*> (zdensity.origin()),
			      z_grid, MPI_DOUBLE, MPI_SUM, bunch.get_comm().get());
    if (error != MPI_SUCCESS) {
        throw std::runtime_error(
                "MPI error in Impedance zdensity");
    }
    
    
    MPI_Allreduce(reinterpret_cast<void*>(local_xmom.origin()),
                   reinterpret_cast<void*>(xmom.origin()),
                                           z_grid, MPI_DOUBLE, MPI_SUM, bunch.get_comm().get());
    if (error != MPI_SUCCESS) {
        throw std::runtime_error(
                "MPI error in Impedance xmom");
    }
  
   
    
    MPI_Allreduce(reinterpret_cast<void*>(local_ymom.origin()),
                   reinterpret_cast<void*>(ymom.origin()),
                                           z_grid, MPI_DOUBLE, MPI_SUM, bunch.get_comm().get());
    if (error != MPI_SUCCESS) {
        throw std::runtime_error(
                "MPI error in Impedance ymom");
    }    

//     //    dbg is set here 
//    int dbg = 0;    
    for (int k = 0; k < z_grid; ++k) {
       // std::cout<<"zdensity[k]="<<zdensity[k]<<std::endl;
        if (zdensity[k] != 0.0) {
         //   if (dbg) std::cout << "before bin: " << k << " zdensity(k): " << zdensity[k] << " xmom(k): " 
	//	  <<xmom[k] << " ymom(k): " <<  ymom[k] << std::endl;
            xmom[k] /= zdensity[k];
            ymom[k] /= zdensity[k];
	 //    if (dbg) std::cout << "after bin: " << k << " zdensity(k): " << zdensity[k] << " xmom(k): " 
	//	  <<xmom[k] << " ymom(k): " <<  ymom[k] << std::endl;
        } else {
            xmom[k] = 0.0;
            ymom[k] = 0.0;
        }
    }
    
}   









inline int get_zindex_for_wake(double z, double dz, int istart, double zstart)
{ //if  (z< (-istart*istart*dz+zstart)) return -100;
  if (z>=zstart){
    return (static_cast<int>(floor(sqrt((z-zstart)/dz)))) +istart;
  }else
  {
    return (-static_cast<int>(floor(sqrt((zstart-z)/dz)))) +istart;
  } 
}  




void Impedance::calculate_kicks(Commxx_sptr const & comm_sptr) {

	double t, t1;
	t = simple_timer_current();

	int zpoints = get_wake_field_sptr()->get_z_coord().size();
	double delta_z = get_wake_field_sptr()->get_delta_z();
	int istart = get_wake_field_sptr()->get_istart();
	double zstart = get_wake_field_sptr()->get_zstart();
	MArray1d_ref z_coord(get_wake_field_sptr()->get_z_coord());
	MArray1d_ref z_wake(get_wake_field_sptr()->get_z_wake());
	MArray1d_ref xw_lead(get_wake_field_sptr()->get_xw_lead());
	MArray1d_ref xw_trail(get_wake_field_sptr()->get_xw_trail());
	MArray1d_ref yw_lead(get_wake_field_sptr()->get_yw_lead());
	MArray1d_ref yw_trail(get_wake_field_sptr()->get_yw_trail());

	t = simple_timer_show(t,
			"impedance_calculate_kicks:  ref the  wake fields ");

	int registered_turns = stored_vbunches.size();
	int numbunches;
	int num_trains = 0;
	if (registered_turns == 0)
		throw std::runtime_error(
				"registered_turns size cannot be zero, probably you propagate a bunch instead of a bunch_train");

	numbunches = (*stored_vbunches.begin()).size();

	if ((full_machine) && (registered_turns != 0)) {
		num_trains = int(num_buckets / numbunches);  
      /// num_trains is relevant only when the full machine option is considered
      /// a full machine consideres a num_train of bunches repeats with modulation wave wn[]
     /// all buckets are full, but only numbunches bunches properties are stored
     /// exemple: full_machine, all bunches identical, no wave across the machine: num_trains=num_buckets, wn=[0,0,0], it's a one bunch simulation  
     /// example: full_machine, two bunch simulation, num_trains= num_buckets/2  

		if (std::abs(num_buckets / float(numbunches) - num_trains) > 1e-8)
			throw std::runtime_error(
					"full machine assumes repetitive numer of trains: num_buckets should be divisible to numbunches");
		if (wn[0] < 0 || wn[0] >= num_trains || wn[1] < 0 || wn[1] >= num_trains
				|| wn[2] < 0 || wn[2] >= num_trains)
			throw std::runtime_error(
					"full machine wave number cannot be smaller than zero or larger than num_trains-1");
	}
    
    
   // std::cout<<" registred turns= "<<registered_turns<<std::endl; 
   // std::cout<<" numbunches= "<<numbunches<<std::endl; 
   // std::cout<<" num_buckets= "<<num_buckets<<std::endl;
   //  std::cout<<" num_trains= "<<num_trains<<std::endl;
   //  std::cout<<"wn="<<wn[0]<<","<<wn[1]<<","<<wn[2]<<std::endl;
   
   MArray1d_ref const xmom(get_xmom());
   MArray1d_ref const ymom(get_ymom());
   MArray1d_ref const zdensity(get_zdensity());
   MArray1d_ref xwake_leading(get_xwake_leading());
   MArray1d_ref xwake_trailing(get_xwake_trailing());
   MArray1d_ref ywake_leading(get_ywake_leading());
   MArray1d_ref ywake_trailing(get_ywake_trailing()); 
   MArray1d_ref zwake0(get_zwake0());
   
   
 //    this->xwake_leading_sptr= boost::shared_ptr<MArray1d >(new MArray1d(boost::extents[z_grid]));
//    this->xwake_trailing_sptr= boost::shared_ptr<MArray1d >(new MArray1d(boost::extents[z_grid]));
//    this->ywake_leading_sptr= boost::shared_ptr<MArray1d >(new MArray1d(boost::extents[z_grid]));
//    this->ywake_trailing_sptr= boost::shared_ptr<MArray1d >(new MArray1d(boost::extents[z_grid]));
//    this->zwake0_sptr= boost::shared_ptr<MArray1d >(new MArray1d(boost::extents[z_grid]));
   
   int size=comm_sptr->get_size();
   std::vector<int >  offsets(size);
   std::vector<int >  counts(size);
   decompose_1d_raw(size, z_grid, offsets, counts);
   int rank=comm_sptr->get_rank();
   
   MArray1d xwake_leading_local(boost::extents[counts[rank]]);
   MArray1d xwake_trailing_local(boost::extents[counts[rank]]);
   MArray1d ywake_leading_local(boost::extents[counts[rank]]);
   MArray1d ywake_trailing_local(boost::extents[counts[rank]]);
   MArray1d zwake0_local(boost::extents[counts[rank]]); 
   
   
   
    t1 = simple_timer_current();
   // for (int i = offsets[rank]; i <offsets[rank]+counts[rank] ; ++i){
    for (int i = 0; i <counts[rank] ; ++i){  
        int real_i=i+offsets[rank];
        xwake_leading_local[i]=0.; 
        xwake_trailing_local[i]=0. ;
        ywake_leading_local[i] =0.; 
        ywake_trailing_local[i]=0. ;
        zwake0_local[i] =0.;	
      // in-bunch impedance 
        for (int j = 0; j < z_grid; ++j){
            double zji=(j-real_i)*cell_size_z;
    
            if (zji>=z_coord[0]) {
                // below it is assumed the wake function is stored using a quadratic grid           
                int iz=get_zindex_for_wake(zji, delta_z, istart, zstart);
                double xwl(0.), xwt(0.), ywl(0.), ywt(0.), zw(0.);
                    if (iz+1 < zpoints) {
                        xwl=xw_lead[iz]+(zji-z_coord[iz])*(xw_lead[iz+1]-xw_lead[iz])/(z_coord[iz+1]-z_coord[iz]);
                        xwt=xw_trail[iz]+(zji-z_coord[iz])*(xw_trail[iz+1]-xw_trail[iz])/(z_coord[iz+1]-z_coord[iz]);
                        ywl=yw_lead[iz]+(zji-z_coord[iz])*(yw_lead[iz+1]-yw_lead[iz])/(z_coord[iz+1]-z_coord[iz]);
                        ywt=yw_trail[iz]+(zji-z_coord[iz])*(yw_trail[iz+1]-yw_trail[iz])/(z_coord[iz+1]-z_coord[iz]);                  
                        zw=z_wake[iz]+(zji-z_coord[iz])*(z_wake[iz+1]-z_wake[iz])/(z_coord[iz+1]-z_coord[iz]); 
                    }
                xwake_leading_local[i]  +=zdensity[j]*N_factor*xmom[j]*xwl; 
                xwake_trailing_local[i]  += zdensity[j]*N_factor*xwt;
                ywake_leading_local[i]  += zdensity[j]*N_factor*ymom[j]*ywl; 
                ywake_trailing_local[i]  += zdensity[j]*N_factor*ywt;
                zwake0_local[i] += zdensity[j]*N_factor*zw;				
            }          
         }
        
        
        std::list< std::vector<Bunch_properties> >::const_iterator it=stored_vbunches.begin(); // stored_vbunches.begin() stores the bunches info at 
                                                                                              // at the moment
        /// bucket 0 is in front of bucket 1, which is in front of bucket 2, etc...
        int mean_bin=static_cast<int>((bunch_z_mean-bunch_z_left)/cell_size_z);
        if ((mean_bin<0) || (mean_bin>=z_grid)){
            throw std::runtime_error(
                    "impedance: the index bin of beam min cannot be <0 or >z_grid, something is wrong  ");
        }
        double z_to_zmean=(mean_bin-real_i)*cell_size_z; 
        for (int ibunch= 0; ibunch<numbunches; ++ibunch){            
            //  double xwl(0.), xwt(0.), ywl(0.), ywt(0.), zw(0.);
            int ibucket=(*it)[ibunch].bucket_index;
            if(ibucket<bunch_bucket) {///  same turn, the leading buckets effect    
                double  zji=z_to_zmean+bunch_spacing*(bunch_bucket-ibucket) +((*it)[ibunch].z_mean-bunch_z_mean);
            //  int iz=static_cast<int>(floor(sqrt((zji-z_coord[0])/(z_coord[1]-z_coord[0]))));  
                int iz=get_zindex_for_wake(zji, delta_z, istart, zstart);
                if ((iz+1 < zpoints) && (iz>0)) {
                    double xwl=xw_lead[iz]+(zji-z_coord[iz])*(xw_lead[iz+1]-xw_lead[iz])/(z_coord[iz+1]-z_coord[iz]);
                    double xwt=xw_trail[iz]+(zji-z_coord[iz])*(xw_trail[iz+1]-xw_trail[iz])/(z_coord[iz+1]-z_coord[iz]);
                    double ywl=yw_lead[iz]+(zji-z_coord[iz])*(yw_lead[iz+1]-yw_lead[iz])/(z_coord[iz+1]-z_coord[iz]);
                    double ywt=yw_trail[iz]+(zji-z_coord[iz])*(yw_trail[iz+1]-yw_trail[iz])/(z_coord[iz+1]-z_coord[iz]);                  
                    double zw=z_wake[iz]+(zji-z_coord[iz])*(z_wake[iz+1]-z_wake[iz])/(z_coord[iz+1]-z_coord[iz]); 

                    xwake_leading_local[i]  +=(*it)[ibunch].realnum*(*it)[ibunch].x_mean*xwl; 
                    xwake_trailing_local[i]  += (*it)[ibunch].realnum*xwt;
                    ywake_leading_local[i]  +=  (*it)[ibunch].realnum*(*it)[ibunch].y_mean*ywl;
                    ywake_trailing_local[i]  +=(*it)[ibunch].realnum*ywt;
                    zwake0_local[i] += (*it)[ibunch].realnum*zw;
                }
            }
        } // ibunch loop

        if (full_machine) { // it assumes that all the other trains are in front 
            for (int itrain= 1; itrain<num_trains; ++itrain){  //for one bunch simulation num_trains = num_buckets, all trains identical except displacement
                double wnx=cos(2.*mconstants::pi*wn[0]*itrain/double(num_trains));
                double wny=cos(2.*mconstants::pi*wn[1]*itrain/double(num_trains));
                double wnz=cos(2.*mconstants::pi*wn[2]*itrain/double(num_trains));
                for (int ibunch= 0; ibunch<numbunches; ++ibunch){		  
                    double  zji=z_to_zmean+bunch_spacing*numbunches*itrain+bunch_spacing*(bunch_bucket-ibunch)+
                    ((*it)[ibunch].z_mean*wnz-bunch_z_mean);
                    //  int iz=static_cast<int>(floor(sqrt((zji-z_coord[0])/(z_coord[1]-z_coord[0])))); 
                    int iz=get_zindex_for_wake(zji, delta_z, istart, zstart);
                    if ((iz+1 < zpoints) && (iz>0)) {
                        double xwl=xw_lead[iz]+(zji-z_coord[iz])*(xw_lead[iz+1]-xw_lead[iz])/(z_coord[iz+1]-z_coord[iz]);
                        double xwt=xw_trail[iz]+(zji-z_coord[iz])*(xw_trail[iz+1]-xw_trail[iz])/(z_coord[iz+1]-z_coord[iz]);
                        double ywl=yw_lead[iz]+(zji-z_coord[iz])*(yw_lead[iz+1]-yw_lead[iz])/(z_coord[iz+1]-z_coord[iz]);
                        double ywt=yw_trail[iz]+(zji-z_coord[iz])*(yw_trail[iz+1]-yw_trail[iz])/(z_coord[iz+1]-z_coord[iz]);                  
                        double zw=z_wake[iz]+(zji-z_coord[iz])*(z_wake[iz+1]-z_wake[iz])/(z_coord[iz+1]-z_coord[iz]);    

                        xwake_leading_local[i]  +=(*it)[ibunch].realnum*(*it)[ibunch].x_mean*wnx*xwl; 
                        xwake_trailing_local[i]  += (*it)[ibunch].realnum*xwt;
                        ywake_leading_local[i]  +=  (*it)[ibunch].realnum*(*it)[ibunch].y_mean*wny*ywl;
                        ywake_trailing_local[i]  +=(*it)[ibunch].realnum*ywt;
                        zwake0_local[i] += (*it)[ibunch].realnum*zw;
                    }
                }
            }
        }  // full_machine
        
        if (registered_turns>1) {
            ++it; ///previous turn, following buckets effect
            for (int ibunch= 0; ibunch<numbunches; ++ibunch){
                int ibucket=(*it)[ibunch].bucket_index;
                if(ibucket>=bunch_bucket) {///  following buckets effect    
                    double  zji=z_to_zmean+bunch_spacing*(bunch_bucket-ibucket)+orbit_length +((*it)[ibunch].z_mean-bunch_z_mean);
                    //  int iz=static_cast<int>(floor(sqrt((zji-z_coord[0])/(z_coord[1]-z_coord[0])))); 
                    int iz=get_zindex_for_wake(zji, delta_z, istart, zstart);
                    if ((iz+1 < zpoints) && (iz>0)) {
                        double xwl=xw_lead[iz]+(zji-z_coord[iz])*(xw_lead[iz+1]-xw_lead[iz])/(z_coord[iz+1]-z_coord[iz]);
                        double xwt=xw_trail[iz]+(zji-z_coord[iz])*(xw_trail[iz+1]-xw_trail[iz])/(z_coord[iz+1]-z_coord[iz]);
                        double ywl=yw_lead[iz]+(zji-z_coord[iz])*(yw_lead[iz+1]-yw_lead[iz])/(z_coord[iz+1]-z_coord[iz]);
                        double ywt=yw_trail[iz]+(zji-z_coord[iz])*(yw_trail[iz+1]-yw_trail[iz])/(z_coord[iz+1]-z_coord[iz]);                  
                        double zw=z_wake[iz]+(zji-z_coord[iz])*(z_wake[iz+1]-z_wake[iz])/(z_coord[iz+1]-z_coord[iz]);    
                        
                        xwake_leading_local[i]  +=(*it)[ibunch].realnum*(*it)[ibunch].x_mean*xwl; 
                        xwake_trailing_local[i]  += (*it)[ibunch].realnum*xwt;
                        ywake_leading_local[i]  +=  (*it)[ibunch].realnum*(*it)[ibunch].y_mean*ywl;
                        ywake_trailing_local[i]  +=(*it)[ibunch].realnum*ywt;
                        zwake0_local[i] += (*it)[ibunch].realnum*zw;                   
                    }
                }
            } // ibunch loop	 
        }//registered_turns>1         
    } // i loop
    t1 = simple_timer_show(t1, "impedance_calculate_kicks:  i loop ");
    /// it is not necessary to have a loop over i at larger distances, since the effect is negligible
    /// might not be  true for coasting beams!
    if (registered_turns>1) {
        double xwake_l=0.;
        double xwake_t=0.;
        double ywake_l=0.;
        double ywake_t=0.;
        double zwake_0=0.; 
              
        
        std::list< std::vector<Bunch_properties> >::const_iterator it;
        std::list< std::vector<Bunch_properties> >::const_iterator jt=stored_vbunches.begin();
        ++jt;       
        int iturn;
        for (it=jt,  iturn=1; it !=stored_vbunches.end(); ++it, ++iturn){
            for (int ibunch= 0; ibunch<numbunches; ++ibunch){
                int ibucket=(*it)[ibunch].bucket_index;                   
                if(((ibucket<bunch_bucket) && (iturn==1))///  finishing the previous turn, for the buckets ahead
                ||  (iturn>1))  /// previous turns effects                
                {
                    double  zji=bunch_spacing*(bunch_bucket-ibucket)+orbit_length*iturn +((*it)[ibunch].z_mean-bunch_z_mean);
                    //  int iz=static_cast<int>(floor(sqrt((zji-z_coord[0])/(z_coord[1]-z_coord[0]))));  
                    int iz=get_zindex_for_wake(zji, delta_z, istart, zstart);
                    if ((iz+1 < zpoints) && (iz>0)) {
                        double  xwl=xw_lead[iz]+(zji-z_coord[iz])*(xw_lead[iz+1]-xw_lead[iz])/(z_coord[iz+1]-z_coord[iz]);
                        double  xwt=xw_trail[iz]+(zji-z_coord[iz])*(xw_trail[iz+1]-xw_trail[iz])/(z_coord[iz+1]-z_coord[iz]);
                        double  ywl=yw_lead[iz]+(zji-z_coord[iz])*(yw_lead[iz+1]-yw_lead[iz])/(z_coord[iz+1]-z_coord[iz]);
                        double  ywt=yw_trail[iz]+(zji-z_coord[iz])*(yw_trail[iz+1]-yw_trail[iz])/(z_coord[iz+1]-z_coord[iz]);                  
                        double  zw=z_wake[iz]+(zji-z_coord[iz])*(z_wake[iz+1]-z_wake[iz])/(z_coord[iz+1]-z_coord[iz]);

                        xwake_l +=  (*it)[ibunch].realnum*(*it)[ibunch].x_mean*xwl;
                        xwake_t += (*it)[ibunch].realnum*xwt;
                        ywake_l +=  (*it)[ibunch].realnum*(*it)[ibunch].y_mean*ywl;
                        ywake_t += (*it)[ibunch].realnum*ywt; 
                        zwake_0 += (*it)[ibunch].realnum*zw;
                    }
                } 
            }

            if (full_machine) {
                for (int itrain= 1; itrain<num_trains; ++itrain){
                    double wnx=cos(2.*mconstants::pi*wn[0]*itrain/double(num_trains));
                    double wny=cos(2.*mconstants::pi*wn[1]*itrain/double(num_trains));
                    double wnz=cos(2.*mconstants::pi*wn[2]*itrain/double(num_trains));
                    for (int ibunch= 0; ibunch<numbunches; ++ibunch){
                        double  zji=orbit_length*iturn+bunch_spacing*numbunches*itrain+bunch_spacing*(bunch_bucket-ibunch)+
                        ((*it)[ibunch].z_mean*wnz-bunch_z_mean);
                        //    int iz=static_cast<int>(floor(sqrt((zji-z_coord[0])/(z_coord[1]-z_coord[0]))));
                        int iz=get_zindex_for_wake(zji, delta_z, istart, zstart);
                        if ((iz+1 < zpoints) && (iz>0)) {
                            double  xwl=xw_lead[iz]+(zji-z_coord[iz])*(xw_lead[iz+1]-xw_lead[iz])/(z_coord[iz+1]-z_coord[iz]);
                            double  xwt=xw_trail[iz]+(zji-z_coord[iz])*(xw_trail[iz+1]-xw_trail[iz])/(z_coord[iz+1]-z_coord[iz]);
                            double  ywl=yw_lead[iz]+(zji-z_coord[iz])*(yw_lead[iz+1]-yw_lead[iz])/(z_coord[iz+1]-z_coord[iz]);
                            double  ywt=yw_trail[iz]+(zji-z_coord[iz])*(yw_trail[iz+1]-yw_trail[iz])/(z_coord[iz+1]-z_coord[iz]);                  
                            double  zw=z_wake[iz]+(zji-z_coord[iz])*(z_wake[iz+1]-z_wake[iz])/(z_coord[iz+1]-z_coord[iz]);

                            xwake_l +=  (*it)[ibunch].realnum*(*it)[ibunch].x_mean*wnx*xwl;
                            xwake_t += (*it)[ibunch].realnum*xwt;
                            ywake_l +=  (*it)[ibunch].realnum*(*it)[ibunch].y_mean*wny*ywl;
                            ywake_t += (*it)[ibunch].realnum*ywt; 
                            zwake_0 += (*it)[ibunch].realnum*zw;
                        }		      
                    }//ibunch
                } ////itrain
            }  // full_machine
      } //iturn
      for (int i = 0; i <counts[rank] ; ++i){
            xwake_leading_local[i] +=  xwake_l;
            xwake_trailing_local[i] += xwake_t;
            ywake_leading_local[i] +=  ywake_l; 
            ywake_trailing_local[i] += ywake_t ;
            zwake0_local[i] += zwake_0;    
      }         
   }//registred_turns>1
    t1 = simple_timer_show(t1, "impedance_calculate_kicks: registred_turns loop   ");
   
     int error = MPI_Allgatherv(reinterpret_cast<void*>(xwake_leading_local.origin()),
		  counts[rank], MPI_DOUBLE,
		  reinterpret_cast<void*>(xwake_leading.origin()),
					  &counts[0], &offsets[0], MPI_DOUBLE, comm_sptr->get());
      if (error != MPI_SUCCESS) {
	  throw std::runtime_error(
	      "MPI error in Impedance_calculate_kicks:MPI_Allgatherv xwake_leading"); 
      } 
 

     error = MPI_Allgatherv(reinterpret_cast<void*>(xwake_trailing_local.origin()),
		  counts[rank], MPI_DOUBLE,
		  reinterpret_cast<void*>(xwake_trailing.origin()),
					  &counts[0], &offsets[0], MPI_DOUBLE, comm_sptr->get());
    if (error != MPI_SUCCESS) {
	  throw std::runtime_error(
	      "MPI error in Impedance_calculate_kicks:MPI_Allgatherv xwake_trailing"); 
     } 
   
    error = MPI_Allgatherv(reinterpret_cast<void*>(ywake_leading_local.origin()),
		  counts[rank], MPI_DOUBLE,
		  reinterpret_cast<void*>(ywake_leading.origin()),
					  &counts[0], &offsets[0], MPI_DOUBLE, comm_sptr->get());
    if (error != MPI_SUCCESS) {
	  throw std::runtime_error(
	      "MPI error in Impedance_calculate_kicks:MPI_Allgatherv ywake_leading"); 
     } 
 

    error = MPI_Allgatherv(reinterpret_cast<void*>(ywake_trailing_local.origin()),
		  counts[rank], MPI_DOUBLE,
		  reinterpret_cast<void*>(ywake_trailing.origin()),
					  &counts[0], &offsets[0], MPI_DOUBLE, comm_sptr->get());
    if (error != MPI_SUCCESS) {
	  throw std::runtime_error(
	      "MPI error in Impedance_calculate_kicks:MPI_Allgatherv ywake_trailing"); 
     } 
   
    error = MPI_Allgatherv(reinterpret_cast<void*>(zwake0_local.origin()),
		  counts[rank], MPI_DOUBLE,
		  reinterpret_cast<void*>(zwake0.origin()),
					  &counts[0], &offsets[0], MPI_DOUBLE, comm_sptr->get());
    if (error != MPI_SUCCESS) {
	  throw std::runtime_error(
	      "MPI error in Impedance_calculate_kicks:MPI_Allgatherv zwake0"); 
     }    
    t1 = simple_timer_show(t1, "impedance_calculate_kicks: mpi_allgather   ");
}

void 
Impedance::apply_impedance_kick(Bunch & bunch, double wake_factor)
{
  
 MArray1i_ref const bin_partition(get_bin_partition());
 MArray1d_ref const  xwake_leading(get_xwake_leading());
 MArray1d_ref const xwake_trailing(get_xwake_trailing());
 MArray1d_ref const  ywake_leading(get_ywake_leading());
 MArray1d_ref const  ywake_trailing(get_ywake_trailing()); 
 MArray1d_ref const zwake0(get_zwake0());
 
 
 for (int part = 0; part < bunch.get_local_num(); ++part) {
        double xkick=0., ykick=0., zkick=0.;
        int bin=bin_partition[part];  // bin_partition(n) is the bin where you find the particle n 
//            if ((bin>=z_grid) || (bin<0))  { std::cout<<"bin="<<bin<<"z_grid="<<z_grid<<std::endl;
//         throw
//         std::runtime_error("something wrong with bining");} 

        xkick=xwake_leading[bin]+xwake_trailing[bin]*bunch.get_local_particles()[part][0];	
        ykick=ywake_leading[bin]+ywake_trailing[bin]*bunch.get_local_particles()[part][2];
        zkick = zwake0[bin];
    
       
        bunch.get_local_particles()[part][1] += wake_factor*xkick;   
        bunch.get_local_particles()[part][3]  += wake_factor*ykick;
        bunch.get_local_particles()[part][5]  += wake_factor*zkick;
      
      
    }
}








void
Impedance::apply(Bunch & bunch, double time_step, Step & step, int verbosity, Logger & logger)
{
   double t;
   t = simple_timer_current();
   bunch.convert_to_state(Bunch::fixed_t_lab);        
   calculate_moments_and_partitions(bunch);   
   t = simple_timer_show(t, "impedance_apply:  calculate_moments_and_partitions ");
      
//       std::ofstream file;
//       file.open("zdensity.dat");
// 	for (int i = 0; i < z_grid; ++i){
// 	file<<i<<"   "<<(*zdensity_sptr)[i]<<"   "<<(*xmom_sptr)[i]<<"   "<<(*ymom_sptr)[i]<<std::endl;
// 	}  
//       file.close();
//       abort(); 
   N_factor=bunch.get_real_num()/bunch.get_total_num();
   bunch_z_mean=Core_diagnostics::calculate_z_mean(bunch);
  
   t = simple_timer_show(t, "impedance_apply: calculate_z_mean ");
   bunch_bucket=bunch.get_bucket_index(); 
   calculate_kicks(bunch.get_comm_sptr());
   t = simple_timer_show(t, "impedance_apply: calculate_kicks ");
   double gamma = bunch.get_reference_particle().get_gamma();
   double beta= bunch.get_reference_particle().get_beta();
   double w_f=get_wake_factor()*time_step/(gamma*beta);  
   
    
   apply_impedance_kick(bunch,  w_f);
   t = simple_timer_show(t, "impedance apply:apply_impedance_kick ");
   
}



void
Impedance::store_bunches_data(Bunch_train & bunch_train)
{       
  
    Bunches bunches(bunch_train.get_bunches());
    size_t num_bunches = bunch_train.get_size();   
    Bunch_properties bi;
    std::vector<Bunch_properties> vbi_local(0);
    std::vector<Bunch_properties> vbi(num_bunches);    
    for (int i = 0; i < num_bunches; ++i){
        if (bunches.at(i)->get_comm().has_this_rank()) {
	    Bunch_sptr bunch_sptr=bunches.at(i);
            bunch_sptr->convert_to_state(Bunch::fixed_t_lab);
	    MArray1d bunch_means=Core_diagnostics::calculate_mean(*bunch_sptr);
	    if (full_machine)  
	      if  (bunch_sptr->get_bucket_index() != i) 				
	                   throw std::runtime_error("for full_machine the buckets have to be occupied in order");		     		    
        bi.x_mean=bunch_means[0];
	    bi.y_mean=bunch_means[2];
	    bi.z_mean=bunch_means[4];
	    bi.realnum=bunch_sptr->get_real_num();
	    bi.bucket_index=bunch_sptr->get_bucket_index();  	
	    if  (bunch_sptr->get_comm().get_rank()==0)   vbi_local.push_back(bi);
                     ///only the rank 0 of every communicator sends the bi to all others
        }  
    }
	    
	     
	  
    MPI_Datatype Bunch_properties_type;
    MPI_Aint lb, extent;
    MPI_Type_get_extent(MPI_DOUBLE, &lb, &extent); 
    MPI_Datatype type[2] = {MPI_DOUBLE, MPI_INT};
    int blocklen[2] = {4,1};
    MPI_Aint disp[2];
    disp[0]=0;
    disp[1]=4*extent;
    MPI_Type_create_struct(2,blocklen, disp, type, &Bunch_properties_type);
    MPI_Type_commit(&Bunch_properties_type); 
				  
    int size_parent_comm=bunch_train.get_parent_comm_sptr()->get_size();	
    std::vector<int > counts(bunch_train.get_proc_counts_for_impedance());
    std::vector<int > offsets(bunch_train.get_proc_offsets_for_impedance());
    
    
    int error = MPI_Allgatherv(reinterpret_cast<void*>(&vbi_local[0]), vbi_local.size(), Bunch_properties_type,  
				    reinterpret_cast<void*>(&vbi[0]), &counts[0], &offsets[0], 
				    Bunch_properties_type,bunch_train.get_parent_comm_sptr()->get() );
    if (error != MPI_SUCCESS) {
      throw std::runtime_error("Impedance::store_bunches_data: MPI error in MPI_Allgatherv");
    } 
		
    MPI_Type_free(&Bunch_properties_type);
    stored_vbunches.push_front(vbi);          
    if (stored_vbunches.size()>nstored_turns) stored_vbunches.pop_back();
         

}


void
Impedance::apply(Bunch_train & bunch_train, double time_step, Step & step,
        int verbosity, Train_diagnosticss const& per_operation_diagnosticss,
        Logger & logger)
{ 
    store_bunches_data(bunch_train);
    Bunches bunches(bunch_train.get_bunches());
    size_t num_bunches = bunch_train.get_size();
    for (int i = 0; i < num_bunches; ++i)
        if (bunches.at(i)->get_comm().has_this_rank()) {
            apply(*bunches.at(i), time_step, step, verbosity,logger);
        }  
} 

void
Impedance::apply(Bunch_train & bunch_train, double time_step, Step & step, int verbosity,
            Train_diagnosticss const& per_operation_train_diagnosticss, 
            Propagate_actions * propagate_actions_ptr, Stepper & stepper, int step_count,  int turn, 
            Logger & logger)
{
   
      apply(bunch_train,time_step,step,verbosity,per_operation_train_diagnosticss,logger);  
}  


Impedance::Impedance()
{
}  

   
#endif
