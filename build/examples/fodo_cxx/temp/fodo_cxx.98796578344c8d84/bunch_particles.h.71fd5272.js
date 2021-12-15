var data = {lines:[
{"lineNum":"    1","line":"#ifndef BUNCH_PARTICLES_H"},
{"lineNum":"    2","line":"#define BUNCH_PARTICLES_H"},
{"lineNum":"    3","line":""},
{"lineNum":"    4","line":""},
{"lineNum":"    5","line":"#include \"synergia/foundation/trigon_traits.h\""},
{"lineNum":"    6","line":"#include \"synergia/utils/multi_array_typedefs.h\""},
{"lineNum":"    7","line":"#include \"synergia/utils/commxx.h\""},
{"lineNum":"    8","line":"#include \"synergia/utils/logger.h\""},
{"lineNum":"    9","line":""},
{"lineNum":"   10","line":"#include \"synergia/utils/hdf5_file.h\""},
{"lineNum":"   11","line":"#include \"synergia/utils/gsvector.h\""},
{"lineNum":"   12","line":""},
{"lineNum":"   13","line":"#include <cereal/cereal.hpp>"},
{"lineNum":"   14","line":""},
{"lineNum":"   15","line":"enum class ParticleGroup"},
{"lineNum":"   16","line":"{"},
{"lineNum":"   17","line":"    regular = 0,"},
{"lineNum":"   18","line":"    spectator = 1"},
{"lineNum":"   19","line":"};"},
{"lineNum":"   20","line":""},
{"lineNum":"   21","line":"typedef Kokkos::View<double*[7],"},
{"lineNum":"   22","line":"        Kokkos::LayoutLeft,"},
{"lineNum":"   23","line":"        Kokkos::DefaultExecutionSpace::memory_space > Particles;"},
{"lineNum":"   24","line":""},
{"lineNum":"   25","line":"typedef Kokkos::View<const double*[7],"},
{"lineNum":"   26","line":"        Kokkos::LayoutLeft,"},
{"lineNum":"   27","line":"        Kokkos::DefaultExecutionSpace::memory_space > ConstParticles;"},
{"lineNum":"   28","line":""},
{"lineNum":"   29","line":"typedef Kokkos::View<uint8_t*,"},
{"lineNum":"   30","line":"        Kokkos::DefaultExecutionSpace::memory_space > ParticleMasks;"},
{"lineNum":"   31","line":""},
{"lineNum":"   32","line":"typedef Kokkos::View<const uint8_t*,"},
{"lineNum":"   33","line":"        Kokkos::DefaultExecutionSpace::memory_space> ConstParticleMasks;"},
{"lineNum":"   34","line":""},
{"lineNum":"   35","line":"#if 0"},
{"lineNum":"   36","line":"typedef Kokkos::View<double*[7],"},
{"lineNum":"   37","line":"        Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace>> Particles;"},
{"lineNum":"   38","line":""},
{"lineNum":"   39","line":"typedef Kokkos::View<const double*[7],"},
{"lineNum":"   40","line":"        Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace>> ConstParticles;"},
{"lineNum":"   41","line":"#endif"},
{"lineNum":"   42","line":""},
{"lineNum":"   43","line":"typedef Particles::HostMirror HostParticles;"},
{"lineNum":"   44","line":"typedef ConstParticles::HostMirror ConstHostParticles;"},
{"lineNum":"   45","line":""},
{"lineNum":"   46","line":"#if 0"},
{"lineNum":"   47","line":"typedef Kokkos::View<uint8_t*> ParticleMasks;"},
{"lineNum":"   48","line":"typedef Kokkos::View<const uint8_t*> ConstParticleMasks;"},
{"lineNum":"   49","line":"#endif"},
{"lineNum":"   50","line":""},
{"lineNum":"   51","line":"typedef ParticleMasks::HostMirror HostParticleMasks;"},
{"lineNum":"   52","line":"typedef ConstParticleMasks::HostMirror ConstHostParticleMasks;"},
{"lineNum":"   53","line":""},
{"lineNum":"   54","line":"// serialization"},
{"lineNum":"   55","line":"namespace cereal"},
{"lineNum":"   56","line":"{"},
{"lineNum":"   57","line":"    // particles"},
{"lineNum":"   58","line":"    template<class AR>"},
{"lineNum":"   59","line":"    void save(AR & ar, Particles const& p)"},
{"lineNum":"   60","line":"    {"},
{"lineNum":"   61","line":"        ar(p.label(), p.stride(1));","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   62","line":"    }"},
{"lineNum":"   63","line":""},
{"lineNum":"   64","line":"    template<class AR>"},
{"lineNum":"   65","line":"    void load(AR & ar, Particles & p)"},
{"lineNum":"   66","line":"    {","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   67","line":"        std::string label;"},
{"lineNum":"   68","line":"        int slots;"},
{"lineNum":"   69","line":"        ar(label, slots);"},
{"lineNum":"   70","line":""},
{"lineNum":"   71","line":"#ifdef NO_PADDING"},
{"lineNum":"   72","line":"        auto alloc = Kokkos::view_alloc(label);"},
{"lineNum":"   73","line":"#else"},
{"lineNum":"   74","line":"        auto alloc = Kokkos::view_alloc(label, Kokkos::AllowPadding);"},
{"lineNum":"   75","line":"#endif"},
{"lineNum":"   76","line":""},
{"lineNum":"   77","line":"        p = Particles(alloc, slots);"},
{"lineNum":"   78","line":""},
{"lineNum":"   79","line":"        if (p.stride(1) != slots)","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   80","line":"            throw std::runtime_error(\"inconsistent padding while loading particles\");","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   81","line":"    }","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   82","line":""},
{"lineNum":"   83","line":"    // masks"},
{"lineNum":"   84","line":"    template<class AR>"},
{"lineNum":"   85","line":"    void save(AR & ar, ParticleMasks const& p)"},
{"lineNum":"   86","line":"    {"},
{"lineNum":"   87","line":"        ar(p.label(), p.extent(0));","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"   88","line":"    }"},
{"lineNum":"   89","line":""},
{"lineNum":"   90","line":"    template<class AR>"},
{"lineNum":"   91","line":"    void load(AR & ar, ParticleMasks & p)"},
{"lineNum":"   92","line":"    {","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   93","line":"        std::string label;"},
{"lineNum":"   94","line":"        int slots;"},
{"lineNum":"   95","line":""},
{"lineNum":"   96","line":"        ar(label, slots);"},
{"lineNum":"   97","line":"        p = ParticleMasks(label, slots);","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   98","line":"    }","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   99","line":"}"},
{"lineNum":"  100","line":""},
{"lineNum":"  101","line":""},
{"lineNum":"  102","line":"template<class PART>"},
{"lineNum":"  103","line":"class bunch_particles_t","class":"lineNoCov","hits":"0","possible_hits":"22",},
{"lineNum":"  104","line":"{"},
{"lineNum":"  105","line":"    using PG = ParticleGroup;"},
{"lineNum":"  106","line":""},
{"lineNum":"  107","line":"public:"},
{"lineNum":"  108","line":""},
{"lineNum":"  109","line":"    constexpr static const int particle_index_null = -1;"},
{"lineNum":"  110","line":""},
{"lineNum":"  111","line":"    using part_t = PART;"},
{"lineNum":"  112","line":""},
{"lineNum":"  113","line":"    using default_memspace ="},
{"lineNum":"  114","line":"        typename Kokkos::DefaultExecutionSpace::memory_space;"},
{"lineNum":"  115","line":""},
{"lineNum":"  116","line":"    using memspace ="},
{"lineNum":"  117","line":"        typename std::conditional<is_trigon<PART>::value,"},
{"lineNum":"  118","line":"                 Kokkos::HostSpace, default_memspace>::type;"},
{"lineNum":"  119","line":""},
{"lineNum":"  120","line":"    using parts_t ="},
{"lineNum":"  121","line":"        typename Kokkos::View<PART*[7],"},
{"lineNum":"  122","line":"                 Kokkos::LayoutLeft, memspace>;"},
{"lineNum":"  123","line":""},
{"lineNum":"  124","line":"    using masks_t ="},
{"lineNum":"  125","line":"        typename Kokkos::View<uint8_t*, memspace>;"},
{"lineNum":"  126","line":""},
{"lineNum":"  127","line":"    using host_parts_t ="},
{"lineNum":"  128","line":"        typename parts_t::HostMirror;"},
{"lineNum":"  129","line":""},
{"lineNum":"  130","line":"    using host_masks_t ="},
{"lineNum":"  131","line":"        typename masks_t::HostMirror;"},
{"lineNum":"  132","line":""},
{"lineNum":"  133","line":"    using const_parts_t ="},
{"lineNum":"  134","line":"        typename Kokkos::View<const PART*[7],"},
{"lineNum":"  135","line":"                 Kokkos::LayoutLeft, memspace>;"},
{"lineNum":"  136","line":""},
{"lineNum":"  137","line":"    using const_masks_t ="},
{"lineNum":"  138","line":"        typename Kokkos::View<const uint8_t*, memspace>;"},
{"lineNum":"  139","line":""},
{"lineNum":"  140","line":"    using const_host_parts_t ="},
{"lineNum":"  141","line":"        typename const_parts_t::HostMirror;"},
{"lineNum":"  142","line":""},
{"lineNum":"  143","line":"    using const_host_masks_t ="},
{"lineNum":"  144","line":"        typename const_masks_t::HostMirror;"},
{"lineNum":"  145","line":""},
{"lineNum":"  146","line":"    using gsv_t ="},
{"lineNum":"  147","line":"        typename std::conditional<is_trigon<PART>::value,"},
{"lineNum":"  148","line":"                 Vec<PART>, GSVector>::type;"},
{"lineNum":"  149","line":""},
{"lineNum":"  150","line":"    using exec_space ="},
{"lineNum":"  151","line":"        typename parts_t::execution_space;"},
{"lineNum":"  152","line":""},
{"lineNum":"  153","line":"private:"},
{"lineNum":"  154","line":""},
{"lineNum":"  155","line":"    /*"},
{"lineNum":"  156","line":"     * Local Particle Array Memory Layout:"},
{"lineNum":"  157","line":"     *"},
{"lineNum":"  158","line":"     *   P: regular particle,"},
{"lineNum":"  159","line":"     *   I: invalid (lost) particle"},
{"lineNum":"  160","line":"     *   R: reserved particle"},
{"lineNum":"  161","line":"     *"},
{"lineNum":"  162","line":"     *    part       mask"},
{"lineNum":"  163","line":"     *   +=====+    +=====+"},
{"lineNum":"  164","line":"     *   |  P  |    |  1  |"},
{"lineNum":"  165","line":"     *   +-----+    +-----+"},
{"lineNum":"  166","line":"     *   |  I  |    |  0  |"},
{"lineNum":"  167","line":"     *   +-----+    +-----+"},
{"lineNum":"  168","line":"     *   |  P  |    |  1  |"},
{"lineNum":"  169","line":"     *   +-----+    +-----+"},
{"lineNum":"  170","line":"     *   |  P  |    |  1  |"},
{"lineNum":"  171","line":"     *   +-----+    +-----+"},
{"lineNum":"  172","line":"     *   |  P  |    |  1  |"},
{"lineNum":"  173","line":"     *   +-----+    +-----+"},
{"lineNum":"  174","line":"     *   |  I  |    |  0  |"},
{"lineNum":"  175","line":"     *   +-----+    +-----+"},
{"lineNum":"  176","line":"     *   |  I  |    |  0  |"},
{"lineNum":"  177","line":"     *   +-----+    +-----+"},
{"lineNum":"  178","line":"     *   |  P  |    |  1  |"},
{"lineNum":"  179","line":"     *   +-----+    +-----+"},
{"lineNum":"  180","line":"     *   |  R  |    |  0  |"},
{"lineNum":"  181","line":"     *   +-----+    +-----+"},
{"lineNum":"  182","line":"     *   |  R  |    |  0  |"},
{"lineNum":"  183","line":"     *   +-----+    +-----+"},
{"lineNum":"  184","line":"     *   |  R  |    |  0  |"},
{"lineNum":"  185","line":"     *   +-----+    +-----+"},
{"lineNum":"  186","line":"     *   |  R  |    |  0  |"},
{"lineNum":"  187","line":"     *   +=====+    +=====+"},
{"lineNum":"  188","line":"     *"},
{"lineNum":"  189","line":"     *   num_valid    = 5   -- num_valid()"},
{"lineNum":"  190","line":"     *   num_active   = 8   -- num_active() / size()"},
{"lineNum":"  191","line":"     *   num_reserved = 12  -- num_reserved() / capacity()"},
{"lineNum":"  192","line":"     *"},
{"lineNum":"  193","line":"     *   when allocating an particle array, the padding flag of the"},
{"lineNum":"  194","line":"     *   Kokkos::View object is turned on. So the \'num_reserved\' is"},
{"lineNum":"  195","line":"     *   the actual array size in the first dimension."},
{"lineNum":"  196","line":"     *"},
{"lineNum":"  197","line":"     *   \'num_valid\' is the number of remainig local particles in"},
{"lineNum":"  198","line":"     *   the bunch"},
{"lineNum":"  199","line":"     *"},
{"lineNum":"  200","line":"     *   \'num_active\' is the number of local particles including lost"},
{"lineNum":"  201","line":"     *   ones. This one should be used when looping through the"},
{"lineNum":"  202","line":"     *   particles array"},
{"lineNum":"  203","line":"     *"},
{"lineNum":"  204","line":"     */"},
{"lineNum":"  205","line":""},
{"lineNum":"  206","line":"    // particle group (regular or spectator)"},
{"lineNum":"  207","line":"    ParticleGroup group;"},
{"lineNum":"  208","line":"    std::string label;"},
{"lineNum":"  209","line":""},
{"lineNum":"  210","line":"    // see the memory layout"},
{"lineNum":"  211","line":"    int n_valid;"},
{"lineNum":"  212","line":"    int n_active;"},
{"lineNum":"  213","line":"    int n_reserved;"},
{"lineNum":"  214","line":""},
{"lineNum":"  215","line":"    // sum of n_valid on all ranks"},
{"lineNum":"  216","line":"    int n_total;"},
{"lineNum":"  217","line":""},
{"lineNum":"  218","line":"    // number of particles discarded from the most recent aperture apply."},
{"lineNum":"  219","line":"    int n_last_discarded;"},
{"lineNum":"  220","line":""},
{"lineNum":"  221","line":"    // particle offset in cases where a single bunch span across"},
{"lineNum":"  222","line":"    // multiple ranks"},
{"lineNum":"  223","line":"    int poffset;"},
{"lineNum":"  224","line":""},
{"lineNum":"  225","line":"public:"},
{"lineNum":"  226","line":""},
{"lineNum":"  227","line":"    parts_t parts;"},
{"lineNum":"  228","line":"    masks_t masks;"},
{"lineNum":"  229","line":"    masks_t discards;"},
{"lineNum":"  230","line":""},
{"lineNum":"  231","line":"    host_parts_t hparts;"},
{"lineNum":"  232","line":"    host_masks_t hmasks;"},
{"lineNum":"  233","line":"    host_masks_t hdiscards;"},
{"lineNum":"  234","line":""},
{"lineNum":"  235","line":"public:"},
{"lineNum":"  236","line":""},
{"lineNum":"  237","line":"    // particles array will be allocated to the size of reserved_num"},
{"lineNum":"  238","line":"    // if reserved is less than total, it will be set to the total_num"},
{"lineNum":"  239","line":"    bunch_particles_t("},
{"lineNum":"  240","line":"            ParticleGroup pg,"},
{"lineNum":"  241","line":"            int total_num,"},
{"lineNum":"  242","line":"            int reserved_num,"},
{"lineNum":"  243","line":"            Commxx const& comm );"},
{"lineNum":"  244","line":""},
{"lineNum":"  245","line":"    int num_valid()          const { return n_valid; }","class":"lineNoCov","hits":"0","possible_hits":"43",},
{"lineNum":"  246","line":"    int num_active()         const { return n_active; }"},
{"lineNum":"  247","line":"    int num_reserved()       const { return n_reserved; }"},
{"lineNum":"  248","line":"    int num_total()          const { return n_total; }","class":"lineNoCov","hits":"0","possible_hits":"6",},
{"lineNum":"  249","line":"    int num_last_discarded() const { return n_last_discarded; }"},
{"lineNum":"  250","line":""},
{"lineNum":"  251","line":"    // getters with names more consistent with std containers"},
{"lineNum":"  252","line":"    int size()        const { return n_active; }","class":"lineNoCov","hits":"0","possible_hits":"31",},
{"lineNum":"  253","line":"    int capacity()    const { return n_reserved; }","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  254","line":"    int size_in_gsv() const { return ceil(1.0*n_active/gsv_t::size()); }","class":"lineNoCov","hits":"0","possible_hits":"48",},
{"lineNum":"  255","line":""},
{"lineNum":"  256","line":"    // copy particles/masks between host and device memories"},
{"lineNum":"  257","line":"    void checkin_particles() const"},
{"lineNum":"  258","line":"    { Kokkos::deep_copy(parts, hparts); Kokkos::deep_copy(masks, hmasks); }","class":"lineNoCov","hits":"0","possible_hits":"14",},
{"lineNum":"  259","line":""},
{"lineNum":"  260","line":"    void checkout_particles() const"},
{"lineNum":"  261","line":"    { Kokkos::deep_copy(hparts, parts); Kokkos::deep_copy(hmasks, masks); }","class":"lineNoCov","hits":"0","possible_hits":"15",},
{"lineNum":"  262","line":""},
{"lineNum":"  263","line":"    // change capacity (can only increase)"},
{"lineNum":"  264","line":"    // n is the new cap across the entire bunch"},
{"lineNum":"  265","line":"    void reserve(int n, Commxx const& comm);"},
{"lineNum":"  266","line":""},
{"lineNum":"  267","line":"    // change local capacity (can only increase)"},
{"lineNum":"  268","line":"    // n is the new local capacity"},
{"lineNum":"  269","line":"    void reserve_local(int n);"},
{"lineNum":"  270","line":""},
{"lineNum":"  271","line":"    // inject with"},
{"lineNum":"  272","line":"    void inject(bunch_particles_t const& o,"},
{"lineNum":"  273","line":"            karray1d_dev const& ref_st_diff,"},
{"lineNum":"  274","line":"            karray1d_dev const& tgt_st,"},
{"lineNum":"  275","line":"            karray1d_dev const& inj_st,"},
{"lineNum":"  276","line":"            double pdiff );"},
{"lineNum":"  277","line":""},
{"lineNum":"  278","line":"    // convert between fixed z lab and fixed t lab"},
{"lineNum":"  279","line":"    void convert_to_fixed_t_lab(double p_ref, double beta);"},
{"lineNum":"  280","line":"    void convert_to_fixed_z_lab(double p_ref, double beta);"},
{"lineNum":"  281","line":""},
{"lineNum":"  282","line":"#if 0"},
{"lineNum":"  283","line":"    void set_total_num(int num);"},
{"lineNum":"  284","line":"    void expand_local_num(int num, int added_lost);"},
{"lineNum":"  285","line":"#endif"},
{"lineNum":"  286","line":""},
{"lineNum":"  287","line":"    // update the valid num from the masks and return the old valid num"},
{"lineNum":"  288","line":"    int update_valid_num();"},
{"lineNum":"  289","line":""},
{"lineNum":"  290","line":"    // update total num across the ranks and returns the old total number"},
{"lineNum":"  291","line":"    int update_total_num(Commxx const& comm);"},
{"lineNum":"  292","line":""},
{"lineNum":"  293","line":"    // apply aperture operation"},
{"lineNum":"  294","line":"    template<typename AP>"},
{"lineNum":"  295","line":"    int apply_aperture(AP const& ap);"},
{"lineNum":"  296","line":""},
{"lineNum":"  297","line":"    // search/get particle(s)"},
{"lineNum":"  298","line":"    int search_particle(int pid, int last_idx) const;"},
{"lineNum":"  299","line":""},
{"lineNum":"  300","line":"    std::pair<karray1d_row, bool>"},
{"lineNum":"  301","line":"    get_particle(int idx) const;"},
{"lineNum":"  302","line":""},
{"lineNum":"  303","line":"    std::pair<karray2d_row, HostParticleMasks>"},
{"lineNum":"  304","line":"    get_particles_in_range(int idx, int num) const;"},
{"lineNum":"  305","line":""},
{"lineNum":"  306","line":"    karray2d_row"},
{"lineNum":"  307","line":"    get_particles_last_discarded() const;"},
{"lineNum":"  308","line":""},
{"lineNum":"  309","line":"    void check_pz2_positive();"},
{"lineNum":"  310","line":"    void print_particle(size_t idx, Logger& logger) const;"},
{"lineNum":"  311","line":""},
{"lineNum":"  312","line":"    // read from a hdf5 file. total_num of current bunch must be the same"},
{"lineNum":"  313","line":"    // as the one stored in the particle file"},
{"lineNum":"  314","line":"    void read_file_legacy(Hdf5_file const& file, Commxx const& comm);"},
{"lineNum":"  315","line":""},
{"lineNum":"  316","line":"    void read_file (Hdf5_file const& file, Commxx const& comm);"},
{"lineNum":"  317","line":"    void write_file(Hdf5_file const& file, int num_part, int offset, Commxx const& comm) const;"},
{"lineNum":"  318","line":""},
{"lineNum":"  319","line":"    // checkpoint save/load"},
{"lineNum":"  320","line":"    void save_checkpoint_particles(Hdf5_file & file, int idx) const;"},
{"lineNum":"  321","line":"    void load_checkpoint_particles(Hdf5_file & file, int idx);"},
{"lineNum":"  322","line":""},
{"lineNum":"  323","line":"    // assign ids cooperatively"},
{"lineNum":"  324","line":"    void assign_ids(int train_idx, int bunch_idx);"},
{"lineNum":"  325","line":""},
{"lineNum":"  326","line":"    // only available for trigons"},
{"lineNum":"  327","line":"    template<class U = PART>"},
{"lineNum":"  328","line":"    std::enable_if_t<is_trigon<U>::value, karray2d_row>"},
{"lineNum":"  329","line":"    get_jacobian(int idx) const"},
{"lineNum":"  330","line":"    {"},
{"lineNum":"  331","line":"        karray2d_row res(\"jacobian\", 6, 6);","class":"lineNoCov","hits":"0","possible_hits":"3",},
{"lineNum":"  332","line":"        for(int i=0; i<6; ++i)","class":"lineNoCov","hits":"0","possible_hits":"6",},
{"lineNum":"  333","line":"            for(int j=0; j<6; ++j)"},
{"lineNum":"  334","line":"                res(i,j) = hparts(idx, i).template get_subpower<1>().terms[j];","class":"lineNoCov","hits":"0","possible_hits":"3",},
{"lineNum":"  335","line":"        return res;"},
{"lineNum":"  336","line":"    }"},
{"lineNum":"  337","line":""},
{"lineNum":"  338","line":"private:"},
{"lineNum":"  339","line":""},
{"lineNum":"  340","line":"    void default_ids(int local_offset, Commxx const& comm);"},
{"lineNum":"  341","line":""},
{"lineNum":"  342","line":"    // serialization"},
{"lineNum":"  343","line":"    friend class cereal::access;"},
{"lineNum":"  344","line":""},
{"lineNum":"  345","line":"    template<class AR>"},
{"lineNum":"  346","line":"    void save(AR & ar) const"},
{"lineNum":"  347","line":"    {","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  348","line":"        ar(CEREAL_NVP(label));"},
{"lineNum":"  349","line":"        ar(CEREAL_NVP(n_valid));"},
{"lineNum":"  350","line":"        ar(CEREAL_NVP(n_active));"},
{"lineNum":"  351","line":"        ar(CEREAL_NVP(n_reserved));"},
{"lineNum":"  352","line":"        ar(CEREAL_NVP(n_total));"},
{"lineNum":"  353","line":""},
{"lineNum":"  354","line":"        ar(CEREAL_NVP(parts));"},
{"lineNum":"  355","line":"        ar(CEREAL_NVP(masks));"},
{"lineNum":"  356","line":"        ar(CEREAL_NVP(discards));"},
{"lineNum":"  357","line":"    }","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  358","line":""},
{"lineNum":"  359","line":"    template<class AR>"},
{"lineNum":"  360","line":"    void load(AR & ar)"},
{"lineNum":"  361","line":"    {","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  362","line":"        ar(CEREAL_NVP(label));","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  363","line":"        ar(CEREAL_NVP(n_valid));"},
{"lineNum":"  364","line":"        ar(CEREAL_NVP(n_active));"},
{"lineNum":"  365","line":"        ar(CEREAL_NVP(n_reserved));"},
{"lineNum":"  366","line":"        ar(CEREAL_NVP(n_total));"},
{"lineNum":"  367","line":""},
{"lineNum":"  368","line":"        ar(CEREAL_NVP(parts));","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  369","line":"        ar(CEREAL_NVP(masks));","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  370","line":"        ar(CEREAL_NVP(discards));","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  371","line":""},
{"lineNum":"  372","line":"        // construct the host array first"},
{"lineNum":"  373","line":"        hparts = Kokkos::create_mirror_view(parts);"},
{"lineNum":"  374","line":"        hmasks = Kokkos::create_mirror_view(masks);"},
{"lineNum":"  375","line":"        hdiscards = Kokkos::create_mirror_view(discards);"},
{"lineNum":"  376","line":"    }","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  377","line":"};"},
{"lineNum":"  378","line":""},
{"lineNum":"  379","line":"// instantiated with a double is still the default BunchParticles"},
{"lineNum":"  380","line":"typedef bunch_particles_t<double> BunchParticles;"},
{"lineNum":"  381","line":""},
{"lineNum":"  382","line":""},
{"lineNum":"  383","line":"// implementations"},
{"lineNum":"  384","line":"namespace bunch_particles_impl"},
{"lineNum":"  385","line":"{"},
{"lineNum":"  386","line":"    template<class AP>"},
{"lineNum":"  387","line":"    struct discard_applier","class":"lineNoCov","hits":"0","possible_hits":"19",},
{"lineNum":"  388","line":"    {"},
{"lineNum":"  389","line":"        typedef int value_type;"},
{"lineNum":"  390","line":""},
{"lineNum":"  391","line":"        AP ap;"},
{"lineNum":"  392","line":"        ConstParticles parts;"},
{"lineNum":"  393","line":"        ParticleMasks masks;"},
{"lineNum":"  394","line":"        ParticleMasks discards;"},
{"lineNum":"  395","line":""},
{"lineNum":"  396","line":"        KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  397","line":"        void operator() (const int i, int& discarded) const"},
{"lineNum":"  398","line":"        {","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  399","line":"            discards(i) = 0;","class":"lineNoCov","hits":"0","possible_hits":"6",},
{"lineNum":"  400","line":""},
{"lineNum":"  401","line":"            if (masks(i) && ap.discard(parts, masks, i))","class":"lineNoCov","hits":"0","possible_hits":"12",},
{"lineNum":"  402","line":"            {"},
{"lineNum":"  403","line":"                discards(i) = 1;","class":"lineNoCov","hits":"0","possible_hits":"6",},
{"lineNum":"  404","line":"                masks(i) = 0;","class":"lineNoCov","hits":"0","possible_hits":"6",},
{"lineNum":"  405","line":"                ++discarded;","class":"lineNoCov","hits":"0","possible_hits":"6",},
{"lineNum":"  406","line":"            }"},
{"lineNum":"  407","line":"        }","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  408","line":"    };"},
{"lineNum":"  409","line":"}"},
{"lineNum":"  410","line":""},
{"lineNum":"  411","line":""},
{"lineNum":"  412","line":"template<>"},
{"lineNum":"  413","line":"template<typename AP>"},
{"lineNum":"  414","line":"inline int bunch_particles_t<double>::apply_aperture(AP const& ap)"},
{"lineNum":"  415","line":"{","class":"lineNoCov","hits":"0","possible_hits":"4",},
{"lineNum":"  416","line":"    using namespace bunch_particles_impl;"},
{"lineNum":"  417","line":""},
{"lineNum":"  418","line":"    int ndiscarded = 0;","class":"lineNoCov","hits":"0","possible_hits":"5",},
{"lineNum":"  419","line":""},
{"lineNum":"  420","line":"    // go through each particle to see which one is been filtered out"},
{"lineNum":"  421","line":"    discard_applier<AP> da{ap, parts, masks, discards};","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  422","line":"    Kokkos::parallel_reduce(n_active, da, ndiscarded);","class":"lineNoCov","hits":"0","possible_hits":"5",},
{"lineNum":"  423","line":""},
{"lineNum":"  424","line":"    //std::cout << \"      discarded = \" << ndiscarded << \"\\n\";"},
{"lineNum":"  425","line":"    n_last_discarded = ndiscarded;","class":"lineNoCov","hits":"0","possible_hits":"5",},
{"lineNum":"  426","line":"    n_valid -= ndiscarded;","class":"lineNoCov","hits":"0","possible_hits":"5",},
{"lineNum":"  427","line":""},
{"lineNum":"  428","line":"    return ndiscarded;"},
{"lineNum":"  429","line":"}","class":"lineNoCov","hits":"0","possible_hits":"8",},
{"lineNum":"  430","line":""},
{"lineNum":"  431","line":"#include \"synergia/utils/parallel_utils.h\""},
{"lineNum":"  432","line":""},
{"lineNum":"  433","line":"namespace bunch_particles_impl"},
{"lineNum":"  434","line":"{"},
{"lineNum":"  435","line":"    struct pid_offset"},
{"lineNum":"  436","line":"    {"},
{"lineNum":"  437","line":"        static int offset;"},
{"lineNum":"  438","line":""},
{"lineNum":"  439","line":"        static int get(int request_num, Commxx const& comm)"},
{"lineNum":"  440","line":"        {"},
{"lineNum":"  441","line":"            MPI_Bcast((void *)&offset, 1, MPI_INT, 0, comm);","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  442","line":"            int old_offset = offset;"},
{"lineNum":"  443","line":"            int total_num;"},
{"lineNum":"  444","line":"            MPI_Reduce((void*)&request_num, (void*)&total_num, 1,","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  445","line":"                    MPI_INT, MPI_SUM, 0, comm);"},
{"lineNum":"  446","line":"            offset += total_num;","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  447","line":"            return old_offset;"},
{"lineNum":"  448","line":"        }"},
{"lineNum":"  449","line":"    };"},
{"lineNum":"  450","line":""},
{"lineNum":"  451","line":"    template<typename parts_t>"},
{"lineNum":"  452","line":"    struct pid_assigner","class":"lineNoCov","hits":"0","possible_hits":"9",},
{"lineNum":"  453","line":"    {"},
{"lineNum":"  454","line":"        parts_t parts;"},
{"lineNum":"  455","line":"        int64_t offset;"},
{"lineNum":"  456","line":""},
{"lineNum":"  457","line":"        KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  458","line":"        void operator() (const int i) const"},
{"lineNum":"  459","line":"        { parts(i, 6) = i + offset; }","class":"lineNoCov","hits":"0","possible_hits":"25",},
{"lineNum":"  460","line":"    };"},
{"lineNum":"  461","line":""},
{"lineNum":"  462","line":"    template<typename masks_t>"},
{"lineNum":"  463","line":"    struct particle_masks_initializer","class":"lineNoCov","hits":"0","possible_hits":"9",},
{"lineNum":"  464","line":"    {"},
{"lineNum":"  465","line":"        masks_t masks;"},
{"lineNum":"  466","line":"        const int num;"},
{"lineNum":"  467","line":""},
{"lineNum":"  468","line":"        KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  469","line":"        void operator() (const int i) const"},
{"lineNum":"  470","line":"        { masks(i) = i<num ? 1 : 0; }","class":"lineNoCov","hits":"0","possible_hits":"30",},
{"lineNum":"  471","line":"    };"},
{"lineNum":"  472","line":"}"},
{"lineNum":"  473","line":""},
{"lineNum":"  474","line":"// default ids only for double typed bunche_particle object"},
{"lineNum":"  475","line":"template<typename PART>"},
{"lineNum":"  476","line":"inline void"},
{"lineNum":"  477","line":"bunch_particles_t<PART>::default_ids("},
{"lineNum":"  478","line":"        int local_offset, Commxx const& comm)"},
{"lineNum":"  479","line":"{","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  480","line":"    using namespace bunch_particles_impl;"},
{"lineNum":"  481","line":""},
{"lineNum":"  482","line":"    int request_num = (comm.rank() == 0) ? n_total : 0;","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  483","line":"    int global_offset = pid_offset::get(request_num, comm);"},
{"lineNum":"  484","line":""},
{"lineNum":"  485","line":"    auto range = Kokkos::RangePolicy<exec_space>(0, n_active);","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  486","line":"    pid_assigner<parts_t> pia{parts, local_offset + global_offset};","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  487","line":"    Kokkos::parallel_for(range, pia);","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  488","line":"}","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  489","line":""},
{"lineNum":"  490","line":"template<typename PART>"},
{"lineNum":"  491","line":"inline"},
{"lineNum":"  492","line":"bunch_particles_t<PART>::bunch_particles_t("},
{"lineNum":"  493","line":"        ParticleGroup pg,"},
{"lineNum":"  494","line":"        int total, int reserved,"},
{"lineNum":"  495","line":"        Commxx const& comm)"},
{"lineNum":"  496","line":"    : group(pg)","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  497","line":"    , label(pg==PG::regular ? \"particles\" : \"spectators\")","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  498","line":"    , n_valid(0)","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  499","line":"    , n_active(0)"},
{"lineNum":"  500","line":"    , n_reserved(0)","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  501","line":"    , n_total(total)","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  502","line":"    , n_last_discarded(0)"},
{"lineNum":"  503","line":"    , poffset(0)","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  504","line":"    , parts()"},
{"lineNum":"  505","line":"    , masks()","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  506","line":"    , discards()","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  507","line":"    , hparts(Kokkos::create_mirror_view(parts))"},
{"lineNum":"  508","line":"    , hmasks(Kokkos::create_mirror_view(masks))"},
{"lineNum":"  509","line":"    , hdiscards(Kokkos::create_mirror_view(discards))"},
{"lineNum":"  510","line":"{","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  511","line":"    using namespace bunch_particles_impl;"},
{"lineNum":"  512","line":""},
{"lineNum":"  513","line":"    // minimum reserved"},
{"lineNum":"  514","line":"    if (reserved < total) reserved = total;","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  515","line":""},
{"lineNum":"  516","line":"    if (!comm.is_null() && reserved)","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  517","line":"    {"},
{"lineNum":"  518","line":"        int mpi_size = comm.size();"},
{"lineNum":"  519","line":"        int mpi_rank = comm.rank();"},
{"lineNum":"  520","line":""},
{"lineNum":"  521","line":"        std::vector<int> offsets_t(mpi_size);","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  522","line":"        std::vector<int> counts_t(mpi_size);"},
{"lineNum":"  523","line":"        decompose_1d(comm, total, offsets_t, counts_t);","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  524","line":""},
{"lineNum":"  525","line":"        std::vector<int> offsets_r(mpi_size);"},
{"lineNum":"  526","line":"        std::vector<int> counts_r(mpi_size);"},
{"lineNum":"  527","line":"        decompose_1d(comm, reserved, offsets_r, counts_r);","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  528","line":""},
{"lineNum":"  529","line":"        // local_num"},
{"lineNum":"  530","line":"        n_valid    = counts_t[mpi_rank];","class":"lineNoCov","hits":"0","possible_hits":"4",},
{"lineNum":"  531","line":"        n_active   = counts_t[mpi_rank];","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  532","line":"        n_reserved = counts_r[mpi_rank];","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  533","line":""},
{"lineNum":"  534","line":"        if (n_active % gsv_t::size())","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  535","line":"        {"},
{"lineNum":"  536","line":"            int padded = n_active + gsv_t::size() - n_active%gsv_t::size();","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  537","line":"            if (n_reserved < padded) n_reserved = padded;","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  538","line":"        }"},
{"lineNum":"  539","line":""},
{"lineNum":"  540","line":"        // local_num offset"},
{"lineNum":"  541","line":"        for(int i=0; i<mpi_rank; ++i)","class":"lineNoCov","hits":"0","possible_hits":"8",},
{"lineNum":"  542","line":"            poffset += counts_t[i];","class":"lineNoCov","hits":"0","possible_hits":"6",},
{"lineNum":"  543","line":""},
{"lineNum":"  544","line":"        // allocate"},
{"lineNum":"  545","line":"#ifdef NO_PADDING"},
{"lineNum":"  546","line":"        auto alloc = Kokkos::view_alloc(label);"},
{"lineNum":"  547","line":"#else"},
{"lineNum":"  548","line":"        auto alloc = Kokkos::view_alloc(label, Kokkos::AllowPadding);"},
{"lineNum":"  549","line":"#endif"},
{"lineNum":"  550","line":"        parts = parts_t(alloc, n_reserved);","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  551","line":""},
{"lineNum":"  552","line":"        // with possible paddings"},
{"lineNum":"  553","line":"        n_reserved = parts.stride(1);","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  554","line":""},
{"lineNum":"  555","line":"        masks = masks_t(label+\"_masks\", n_reserved);","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  556","line":"        discards = masks_t(label+\"_discards\", n_reserved);","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  557","line":""},
{"lineNum":"  558","line":"        hparts = Kokkos::create_mirror_view(parts);"},
{"lineNum":"  559","line":"        hmasks = Kokkos::create_mirror_view(masks);"},
{"lineNum":"  560","line":"        hdiscards = Kokkos::create_mirror_view(discards);"},
{"lineNum":"  561","line":""},
{"lineNum":"  562","line":"        // set default ids"},
{"lineNum":"  563","line":"        default_ids(offsets_t[mpi_rank], comm);","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  564","line":""},
{"lineNum":"  565","line":"        // valid particles"},
{"lineNum":"  566","line":"        auto range = Kokkos::RangePolicy<exec_space>(0, n_reserved);","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  567","line":"        particle_masks_initializer<masks_t> pmi{masks, n_valid};","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  568","line":"        Kokkos::parallel_for(range, pmi);","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  569","line":""},
{"lineNum":"  570","line":"        // sync host arrays with device arrays"},
{"lineNum":"  571","line":"        checkout_particles();"},
{"lineNum":"  572","line":"    }"},
{"lineNum":"  573","line":"}","class":"lineNoCov","hits":"0","possible_hits":"5",},
{"lineNum":"  574","line":""},
{"lineNum":"  575","line":""},
{"lineNum":"  576","line":"#endif"},
]};
var percent_low = 25;var percent_high = 75;
var header = { "command" : "fodo_cxx", "date" : "2021-12-15 09:03:42", "instrumented" : 87, "covered" : 0,};
var merged_data = [];
