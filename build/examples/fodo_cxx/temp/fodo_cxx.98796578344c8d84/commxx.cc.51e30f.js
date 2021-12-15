var data = {lines:[
{"lineNum":"    1","line":"#include \"commxx.h\""},
{"lineNum":"    2","line":"#include <stdexcept>"},
{"lineNum":"    3","line":"#include <climits>"},
{"lineNum":"    4","line":"#include <algorithm>"},
{"lineNum":"    5","line":""},
{"lineNum":"    6","line":""},
{"lineNum":"    7","line":"namespace"},
{"lineNum":"    8","line":"{"},
{"lineNum":"    9","line":"    struct comm_free"},
{"lineNum":"   10","line":"    {"},
{"lineNum":"   11","line":"        void operator() (MPI_Comm * comm)"},
{"lineNum":"   12","line":"        {"},
{"lineNum":"   13","line":"            if ( comm == nullptr) return;","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   14","line":"            if (*comm == MPI_COMM_NULL) return;","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   15","line":""},
{"lineNum":"   16","line":"            MPI_Comm_free(comm);","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"   17","line":"            delete comm;","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"   18","line":"        }"},
{"lineNum":"   19","line":"    };"},
{"lineNum":"   20","line":""},
{"lineNum":"   21","line":"    size_t hash(const char * s)"},
{"lineNum":"   22","line":"    {"},
{"lineNum":"   23","line":"        size_t h = 37062913;"},
{"lineNum":"   24","line":"        while (*s) h = h * 101 + (unsigned char) *s++;"},
{"lineNum":"   25","line":"        return h;"},
{"lineNum":"   26","line":"    }"},
{"lineNum":"   27","line":"}"},
{"lineNum":"   28","line":""},
{"lineNum":"   29","line":""},
{"lineNum":"   30","line":""},
{"lineNum":"   31","line":"#if 0"},
{"lineNum":"   32","line":"void"},
{"lineNum":"   33","line":"Commxx::construct(MPI_Comm const& parent_mpi_comm)"},
{"lineNum":"   34","line":"{"},
{"lineNum":"   35","line":"    MPI_Comm temp_comm;"},
{"lineNum":"   36","line":"    int error;"},
{"lineNum":"   37","line":""},
{"lineNum":"   38","line":"    if (ranks.size() > 0) {"},
{"lineNum":"   39","line":"        MPI_Group parent_group, group;"},
{"lineNum":"   40","line":"        error = MPI_Comm_group(parent_mpi_comm, &parent_group);"},
{"lineNum":"   41","line":"        if (error != MPI_SUCCESS) {"},
{"lineNum":"   42","line":"            throw std::runtime_error(\"MPI error in Commxx(MPI_Comm_group)\");"},
{"lineNum":"   43","line":"        }"},
{"lineNum":"   44","line":"        error = MPI_Group_incl(parent_group, ranks.size(), &ranks[0], &group);"},
{"lineNum":"   45","line":"        if (error != MPI_SUCCESS) {"},
{"lineNum":"   46","line":"            throw std::runtime_error(\"MPI error in Commxx(MPI_Group_incl)\");"},
{"lineNum":"   47","line":"        }"},
{"lineNum":"   48","line":"        error = MPI_Comm_create(parent_mpi_comm, group, &temp_comm);"},
{"lineNum":"   49","line":"        if (error != MPI_SUCCESS) {"},
{"lineNum":"   50","line":"            throw std::runtime_error(\"MPI error in Commxx(MPI_Comm_create)\");"},
{"lineNum":"   51","line":"        }"},
{"lineNum":"   52","line":"        error = MPI_Group_free(&parent_group);"},
{"lineNum":"   53","line":"        if (error != MPI_SUCCESS) {"},
{"lineNum":"   54","line":"            throw std::runtime_error("},
{"lineNum":"   55","line":"                    \"MPI error in Commxx(MPI_Group_free(parent))\");"},
{"lineNum":"   56","line":"        }"},
{"lineNum":"   57","line":"        error = MPI_Group_free(&group);"},
{"lineNum":"   58","line":"        if (error != MPI_SUCCESS) {"},
{"lineNum":"   59","line":"            throw std::runtime_error(\"MPI error in Commxx(MPI_Group_free)\");"},
{"lineNum":"   60","line":"        }"},
{"lineNum":"   61","line":""},
{"lineNum":"   62","line":"        has_this_rank_ = false;"},
{"lineNum":"   63","line":"        const int no_rank = -1;"},
{"lineNum":"   64","line":"        int this_rank = no_rank;"},
{"lineNum":"   65","line":"        if (!parent_sptr) {"},
{"lineNum":"   66","line":"            this_rank = Commxx().get_rank();"},
{"lineNum":"   67","line":"        } else {"},
{"lineNum":"   68","line":"            if (parent_sptr->has_this_rank()) {"},
{"lineNum":"   69","line":"                this_rank = parent_sptr->get_rank();"},
{"lineNum":"   70","line":"            }"},
{"lineNum":"   71","line":"        }"},
{"lineNum":"   72","line":"        if (this_rank != no_rank) {"},
{"lineNum":"   73","line":"            for (std::vector<int >::const_iterator it = ranks.begin(); it"},
{"lineNum":"   74","line":"                    != ranks.end(); ++it) {"},
{"lineNum":"   75","line":"                if ((*it) == this_rank) {"},
{"lineNum":"   76","line":"                    has_this_rank_ = true;"},
{"lineNum":"   77","line":"                }"},
{"lineNum":"   78","line":"            }"},
{"lineNum":"   79","line":"        }"},
{"lineNum":"   80","line":""},
{"lineNum":"   81","line":"    } else {"},
{"lineNum":"   82","line":"        temp_comm = parent_mpi_comm;"},
{"lineNum":"   83","line":"        has_this_rank_ = true;"},
{"lineNum":"   84","line":"    }"},
{"lineNum":"   85","line":""},
{"lineNum":"   86","line":"    if (per_host && has_this_rank_) {"},
{"lineNum":"   87","line":"        char name[MPI_MAX_PROCESSOR_NAME];"},
{"lineNum":"   88","line":"        int name_len;"},
{"lineNum":"   89","line":"        MPI_Get_processor_name(name, &name_len);"},
{"lineNum":"   90","line":""},
{"lineNum":"   91","line":"        int color = hash(name) % INT_MAX;"},
{"lineNum":"   92","line":""},
{"lineNum":"   93","line":"        int result = MPI_Comm_split(temp_comm, color, 0, &comm);"},
{"lineNum":"   94","line":"        if (result != MPI_SUCCESS) throw std::runtime_error("},
{"lineNum":"   95","line":"                \"MPI error in MPI_Comm_split\");"},
{"lineNum":"   96","line":"        if (ranks.size() > 0) {"},
{"lineNum":"   97","line":"            error = MPI_Comm_free(&temp_comm);"},
{"lineNum":"   98","line":"            if (error != MPI_SUCCESS) {"},
{"lineNum":"   99","line":"                throw std::runtime_error("},
{"lineNum":"  100","line":"                        \"MPI error in Commxx(MPI_Comm_free(temp_comm))\");"},
{"lineNum":"  101","line":"            }"},
{"lineNum":"  102","line":"        }"},
{"lineNum":"  103","line":"    } else {"},
{"lineNum":"  104","line":"        comm = temp_comm;"},
{"lineNum":"  105","line":"    }"},
{"lineNum":"  106","line":"}"},
{"lineNum":"  107","line":"#endif"},
{"lineNum":"  108","line":""},
{"lineNum":"  109","line":"const Commxx Commxx::World(comm_type::world);"},
{"lineNum":"  110","line":"const Commxx Commxx::Null(comm_type::null);"},
{"lineNum":"  111","line":""},
{"lineNum":"  112","line":"Commxx::Commxx(comm_type type)"},
{"lineNum":"  113","line":"    : comm(type==comm_type::null ? nullptr : new MPI_Comm(MPI_COMM_WORLD))","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  114","line":"    , parent_comm()"},
{"lineNum":"  115","line":"    , type(type)","class":"lineNoCov","hits":"0","possible_hits":"3",},
{"lineNum":"  116","line":"    , color(0)","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  117","line":"    , key(0)","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  118","line":"{","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  119","line":"    if (type == comm_type::regular)","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  120","line":"        throw std::runtime_error(\"only null or mpi_comm_world can be created\");","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  121","line":"}","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  122","line":""},
{"lineNum":"  123","line":"void Commxx::construct()"},
{"lineNum":"  124","line":"{","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  125","line":"    if (!parent_comm || parent_comm->is_null())","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  126","line":"        throw std::runtime_error(\"invalid parent communicator while in construct\");","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  127","line":""},
{"lineNum":"  128","line":"    MPI_Comm newcomm;"},
{"lineNum":"  129","line":"    MPI_Comm_split(*(parent_comm->comm), color, key, &newcomm);","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  130","line":""},
{"lineNum":"  131","line":"    // do not construct the null communicator"},
{"lineNum":"  132","line":"    // always take the ownership"},
{"lineNum":"  133","line":"    if (newcomm != MPI_COMM_NULL)","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  134","line":"        comm.reset(new MPI_Comm(newcomm), comm_free());","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  135","line":"}","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  136","line":""},
{"lineNum":"  137","line":"Commxx::Commxx(std::shared_ptr<const Commxx> const& parent, int color, int key)"},
{"lineNum":"  138","line":"    : comm()","class":"lineNoCov","hits":"0","possible_hits":"4",},
{"lineNum":"  139","line":"    , parent_comm(std::move(parent))","class":"lineNoCov","hits":"0","possible_hits":"4",},
{"lineNum":"  140","line":"    , type(comm_type::regular)","class":"lineNoCov","hits":"0","possible_hits":"5",},
{"lineNum":"  141","line":"    , color(color)","class":"lineNoCov","hits":"0","possible_hits":"4",},
{"lineNum":"  142","line":"    , key(key)","class":"lineNoCov","hits":"0","possible_hits":"5",},
{"lineNum":"  143","line":"{","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  144","line":"    construct();","class":"lineNoCov","hits":"0","possible_hits":"5",},
{"lineNum":"  145","line":"}","class":"lineNoCov","hits":"0","possible_hits":"5",},
{"lineNum":"  146","line":""},
{"lineNum":"  147","line":"int"},
{"lineNum":"  148","line":"Commxx::get_rank() const"},
{"lineNum":"  149","line":"{","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  150","line":"    if (type == comm_type::null)","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  151","line":"        throw std::runtime_error(\"Cannot get_rank() for a null commxx\");","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  152","line":""},
{"lineNum":"  153","line":"    int error, rank;"},
{"lineNum":"  154","line":"    error = MPI_Comm_rank(*comm, &rank);","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  155","line":"    if (error != MPI_SUCCESS) {","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  156","line":"        throw std::runtime_error(\"MPI error in MPI_Comm_rank\");","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  157","line":"    }"},
{"lineNum":"  158","line":"    return rank;","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  159","line":"}"},
{"lineNum":"  160","line":""},
{"lineNum":"  161","line":"int"},
{"lineNum":"  162","line":"Commxx::get_size() const"},
{"lineNum":"  163","line":"{","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  164","line":"    if (type == comm_type::null)","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  165","line":"        throw std::runtime_error(\"Cannot get_size() for a null commxx\");","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  166","line":""},
{"lineNum":"  167","line":"    int error, size;"},
{"lineNum":"  168","line":"    error = MPI_Comm_size(*comm, &size);","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  169","line":"    if (error != MPI_SUCCESS) {","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  170","line":"        throw std::runtime_error(\"MPI error in MPI_Comm_size\");","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  171","line":"    }"},
{"lineNum":"  172","line":"    return size;","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  173","line":"}"},
{"lineNum":"  174","line":""},
{"lineNum":"  175","line":"Commxx Commxx::parent() const"},
{"lineNum":"  176","line":"{","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  177","line":"    if (!parent_comm)","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  178","line":"        throw std::runtime_error(\"Invalid parent communicator\");","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  179","line":""},
{"lineNum":"  180","line":"    return *parent_comm;","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  181","line":"}"},
{"lineNum":"  182","line":""},
{"lineNum":"  183","line":"Commxx Commxx::dup() const"},
{"lineNum":"  184","line":"{","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  185","line":"    if (is_null()) throw std::runtime_error(\"dup from a null comm\");","class":"lineNoCov","hits":"0","possible_hits":"3",},
{"lineNum":"  186","line":"    return Commxx(shared_from_this(), 0, rank());","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  187","line":"}"},
{"lineNum":"  188","line":""},
{"lineNum":"  189","line":"Commxx Commxx::split(int color) const"},
{"lineNum":"  190","line":"{","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  191","line":"    if (is_null()) throw std::runtime_error(\"split from a null comm\");","class":"lineNoCov","hits":"0","possible_hits":"3",},
{"lineNum":"  192","line":"    return Commxx(shared_from_this(), color, rank());","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  193","line":"}"},
{"lineNum":"  194","line":""},
{"lineNum":"  195","line":"Commxx Commxx::split(int color, int key) const"},
{"lineNum":"  196","line":"{","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  197","line":"    if (is_null()) throw std::runtime_error(\"split from a null comm\");","class":"lineNoCov","hits":"0","possible_hits":"3",},
{"lineNum":"  198","line":"    return Commxx(shared_from_this(), color, key);","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  199","line":"}"},
{"lineNum":"  200","line":""},
{"lineNum":"  201","line":"Commxx Commxx::divide(int subgroup_size) const"},
{"lineNum":"  202","line":"{","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  203","line":"    if (is_null()) throw std::runtime_error(\"divide from a null comm\");","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  204","line":""},
{"lineNum":"  205","line":"    if (size() < subgroup_size) return dup();","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  206","line":"    if (size() % subgroup_size != 0)","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  207","line":"    {"},
{"lineNum":"  208","line":"        throw std::runtime_error(","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  209","line":"                \"Commxx::divide(): size must be divisible by subgroup_size\" );"},
{"lineNum":"  210","line":"    }"},
{"lineNum":"  211","line":""},
{"lineNum":"  212","line":"    int color = rank() / subgroup_size;","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  213","line":"    return split(color);","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  214","line":"}","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  215","line":""},
{"lineNum":"  216","line":"Commxx Commxx::group(std::vector<int> const & ranks) const"},
{"lineNum":"  217","line":"{","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  218","line":"    if (is_null()) throw std::runtime_error(\"group from a null comm\");","class":"lineNoCov","hits":"0","possible_hits":"3",},
{"lineNum":"  219","line":""},
{"lineNum":"  220","line":"    int r = rank();"},
{"lineNum":"  221","line":"    int color = (std::find(ranks.begin(), ranks.end(), r) != ranks.end())","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  222","line":"        ? 0 : MPI_UNDEFINED;"},
{"lineNum":"  223","line":""},
{"lineNum":"  224","line":"    return Commxx(shared_from_this(), color, r);","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  225","line":""},
{"lineNum":"  226","line":"#if 0"},
{"lineNum":"  227","line":"    MPI_Group grp;"},
{"lineNum":"  228","line":"    MPI_Comm_group(MPI_Comm(*this), &grp);"},
{"lineNum":"  229","line":""},
{"lineNum":"  230","line":"    MPI_Group subgrp;"},
{"lineNum":"  231","line":"    MPI_Group_incl(grp, ranks.size(), &ranks[0], &subgrp);"},
{"lineNum":"  232","line":""},
{"lineNum":"  233","line":"    MPI_Comm subcomm;"},
{"lineNum":"  234","line":"    MPI_Comm_create(MPI_Comm(*this), subgrp, &subcomm);"},
{"lineNum":"  235","line":""},
{"lineNum":"  236","line":"    MPI_Group_free(&grp);"},
{"lineNum":"  237","line":"    MPI_Group_free(&subgrp);"},
{"lineNum":"  238","line":""},
{"lineNum":"  239","line":"    return Commxx(subcomm, comm_create_kind::take_ownership);"},
{"lineNum":"  240","line":"#endif"},
{"lineNum":"  241","line":"}","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  242","line":""},
{"lineNum":"  243","line":"bool operator== (Commxx const & comm1, Commxx const & comm2)"},
{"lineNum":"  244","line":"{","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  245","line":"    // both null"},
{"lineNum":"  246","line":"    if (comm1.is_null() && comm2.is_null())","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  247","line":"        return true;"},
{"lineNum":"  248","line":""},
{"lineNum":"  249","line":"    // one is null, the other is not"},
{"lineNum":"  250","line":"    if (comm1.is_null() || comm2.is_null())"},
{"lineNum":"  251","line":"        return false;"},
{"lineNum":"  252","line":""},
{"lineNum":"  253","line":"    // both not null"},
{"lineNum":"  254","line":"    int result;"},
{"lineNum":"  255","line":"    MPI_Comm_compare((MPI_Comm)comm1, (MPI_Comm)comm2, &result);","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  256","line":"    return result == MPI_IDENT;","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  257","line":"}","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  258","line":""},
{"lineNum":"  259","line":"#if 0"},
{"lineNum":"  260","line":"Commxxs"},
{"lineNum":"  261","line":"generate_subcomms(Commxx_sptr parent_sptr, int count)"},
{"lineNum":"  262","line":"{"},
{"lineNum":"  263","line":"    Commxxs retval(0);"},
{"lineNum":"  264","line":"    std::vector<std::vector<int > > ranks(distribute_1d(*parent_sptr, count));"},
{"lineNum":"  265","line":"    for (int index = 0; index < count; ++index) {"},
{"lineNum":"  266","line":"        if (index == 0) {"},
{"lineNum":"  267","line":"            retval.push_back("},
{"lineNum":"  268","line":"                    Commxx_sptr(new Commxx(parent_sptr, ranks.at(index))));"},
{"lineNum":"  269","line":"        } else if ((ranks.at(index) == ranks.at(index - 1))) {"},
{"lineNum":"  270","line":"            retval.push_back(retval.at(index - 1));"},
{"lineNum":"  271","line":"        } else {"},
{"lineNum":"  272","line":"            retval.push_back("},
{"lineNum":"  273","line":"                    Commxx_sptr(new Commxx(parent_sptr, ranks.at(index))));"},
{"lineNum":"  274","line":"        }"},
{"lineNum":"  275","line":"    }"},
{"lineNum":"  276","line":"    return retval;"},
{"lineNum":"  277","line":"}"},
{"lineNum":"  278","line":""},
{"lineNum":"  279","line":"Commxx_sptr"},
{"lineNum":"  280","line":"make_optimal_spc_comm(Commxx_sptr comm_sptr, int optimal_number, bool equally_spread)"},
{"lineNum":"  281","line":"{"},
{"lineNum":"  282","line":"    int optimal_size=std::min(optimal_number, comm_sptr->get_size());"},
{"lineNum":"  283","line":"    std::vector<int > on_ranks(optimal_size);"},
{"lineNum":"  284","line":"    int start_rank;"},
{"lineNum":"  285","line":"     if (equally_spread){"},
{"lineNum":"  286","line":"        if ((comm_sptr->get_size() %  optimal_size) !=0)"},
{"lineNum":"  287","line":"\t  throw std::runtime_error(\"make_optimal_spc_comm, for equal_spread  the subsize is not a divider of comm size\");"},
{"lineNum":"  288","line":"\tstart_rank=comm_sptr->get_rank()/optimal_size;"},
{"lineNum":"  289","line":"     }"},
{"lineNum":"  290","line":"     else{"},
{"lineNum":"  291","line":"       start_rank=0;"},
{"lineNum":"  292","line":"     }"},
{"lineNum":"  293","line":"    for (int i=0; i<optimal_size;++i){"},
{"lineNum":"  294","line":"\t    on_ranks[i]=i+start_rank*optimal_size;"},
{"lineNum":"  295","line":"    }"},
{"lineNum":"  296","line":"    Commxx_sptr ret_comm_sptr(new Commxx(comm_sptr,on_ranks));"},
{"lineNum":"  297","line":"    return ret_comm_sptr;"},
{"lineNum":"  298","line":"}"},
{"lineNum":"  299","line":"#endif"},
]};
var percent_low = 25;var percent_high = 75;
var header = { "command" : "fodo_cxx", "date" : "2021-12-15 09:03:41", "instrumented" : 72, "covered" : 0,};
var merged_data = [];
