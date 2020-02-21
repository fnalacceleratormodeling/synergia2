#include "synergia/simulation/checkpoint.h"

#include "synergia/simulation/propagator.h"
#include "synergia/simulation/bunch_simulator.h"
#include "synergia/utils/json.h"



void 
syn::checkpoint_save( Propagator const& prop, 
                      Bunch_simulator const& sim )
{
    std::string prop_str = prop.dump();
    std::string sim_str = sim.dump();

    // collect sim_str to the root rank
    auto const& comm = sim.get_comm();

    const int root = 0;
    const int mpi_size = comm.size();
    const int mpi_rank = comm.rank();

    int len = sim_str.size();
    std::vector<int> lens(mpi_size, 0);

    // gather string size
    MPI_Gather(&len, 1, MPI_INT, lens.data(), 1, MPI_INT, root, comm);

    // accumulate sizes
    std::vector<int> displs(mpi_size, 0);
    for(int i=1; i<mpi_size; ++i) displs[i] = displs[i-1] + lens[i-1];

    // recv string buffer
    int total_len = displs[mpi_size-1] + lens[mpi_size-1];
    std::string sims_str(total_len, ' ');
    
    // gather strings
    MPI_Gatherv( sim_str.data(), sim_str.size(), MPI_CHAR,
            (void*)sims_str.data(), lens.data(), displs.data(), MPI_CHAR, 
            root, comm );

    // extract each string and parse into a JSON object
    if (mpi_rank == root)
    {
        syn::json cp = syn::json::object();

        cp["propagator"] = syn::json::parse(prop_str);
        cp["simulator"] = syn::json::array();

        for(int i=0; i<mpi_size; ++i)
        {
            auto begin = sims_str.begin() + displs[i];
            auto d = syn::json::parse(begin, begin + lens[i]);
            cp["simulator"].push_back(d);
        }

        // write states to file
        std::ofstream file("cp_state.json");
        if (!file.good()) throw std::runtime_error(
                "Error at creating checkpointing file");
        file << cp;
    }
}

std::pair<Propagator, Bunch_simulator>
syn::checkpoint_load()
{
    const int root = 0;
    const int mpi_size = Commxx::world_size();
    const int mpi_rank = Commxx::world_rank();

    size_t len;
    std::vector<char> buf;

    if (mpi_rank == root)
    {
        // read states from file
        std::ifstream file("cp_state.json");
        if (!file.good()) throw std::runtime_error(
                "Error at openning checkpointing file");

        file.seekg(0, std::ios::end);
        len = file.tellg();

        // broadcast the buffer length
        MPI_Bcast(&len, 1, MPI_UINT64_T, root, Commxx::World);

        // read the buffer
        buf.resize(len);
        file.seekg(0);
        file.read(&buf[0], len);

        // broadcast the buffer
        MPI_Bcast(&buf[0], len, MPI_BYTE, root, Commxx::World);
    }
    else
    {
        // receive the buffer length
        MPI_Bcast(&len, 1, MPI_UINT64_T, root, Commxx::World);

        // prepare buffer
        buf.resize(len);

        // receive the buffer
        MPI_Bcast(&buf[0], len, MPI_BYTE, root, Commxx::World);
    }

    // parse the json object
    auto cp = syn::json::parse(buf.begin(), buf.end());

    // recreate the objects
    return std::make_pair(
            Propagator::load_from_string(cp["propagator"].dump()),
            Bunch_simulator::load_from_string(cp["simulator"][mpi_rank].dump()) );
}

void syn::resume()
{
    Logger simlog(0, LoggerV::INFO_STEP);

    auto cp = syn::checkpoint_load();

    Propagator& propagator = cp.first;
    Bunch_simulator& sim = cp.second;

    propagator.propagate(sim, simlog);
}
