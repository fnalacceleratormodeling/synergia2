#include "synergia/utils/catch.hpp"

#include "synergia/utils/commxx.h"
#include "synergia/utils/cereal_files.h"



TEST_CASE("default construct", "[commxx]")
{
    REQUIRE_NOTHROW(Commxx());
}

TEST_CASE("construct with type", "[commxx]")
{
    REQUIRE_NOTHROW(Commxx(comm_type::world));
    REQUIRE_NOTHROW(Commxx(comm_type::null));

    // creating a comm_type::regular is not allowed
    REQUIRE_THROWS(Commxx(comm_type::regular));
}

TEST_CASE("static construct", "[commxx]")
{
    SECTION("null")
    {
        auto null = Commxx::Null;
        Commxx comm(comm_type::null);

        REQUIRE(null == comm);
    }

    SECTION("world")
    {
        auto world = Commxx::World;
        Commxx comm(comm_type::world);

        REQUIRE(world == comm);
    }
}

TEST_CASE("construct world", "[commxx]")
{
    SECTION("world")
    {
        Commxx comm(comm_type::world);
        CHECK(comm == MPI_COMM_WORLD);
    }

    SECTION("default world")
    {
        Commxx comm;
        CHECK(comm == MPI_COMM_WORLD);
    }
}

TEST_CASE("construct null", "[commxx]")
{
    Commxx comm(comm_type::null);
    CHECK(comm == MPI_COMM_NULL);
}

TEST_CASE("get_type")
{
    SECTION("world")
    {
        Commxx comm(comm_type::world);
        CHECK(comm.get_type() == comm_type::world);
    }

    SECTION("null")
    {
        Commxx comm(comm_type::null);
        CHECK(comm.get_type() == comm_type::null);
    }
}

TEST_CASE("compare")
{
    SECTION("null to null")
    {
        auto comm1 = Commxx::Null;
        auto comm2 = Commxx::Null;
        REQUIRE(comm1 == comm2);
    }

    SECTION("null to world")
    {
        auto comm1 = Commxx::Null;
        auto comm2 = Commxx::World;
        REQUIRE(comm1 != comm2);
    }

    SECTION("world to world")
    {
        auto comm1 = Commxx::World;
        auto comm2 = Commxx::World;
        REQUIRE(comm1 == comm2);
    }
}

TEST_CASE("get_size")
{
    SECTION("null")
    {
        Commxx comm(comm_type::null);
        REQUIRE_THROWS(comm.get_size());
        REQUIRE_THROWS(comm.size());
    }

    SECTION("world")
    {
        // get size from MPI_COMM_WORLD
        int world_size;
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);

        Commxx comm(comm_type::world);
        REQUIRE_NOTHROW(comm.get_size());
        REQUIRE_NOTHROW(comm.size());

        CHECK(world_size == comm.get_size());
        CHECK(world_size == comm.size());
    }
}

TEST_CASE("get_rank")
{
    SECTION("null")
    {
        Commxx comm(comm_type::null);
        REQUIRE_THROWS(comm.get_rank());
        REQUIRE_THROWS(comm.rank());
    }

    SECTION("world")
    {
        // get size from MPI_COMM_WORLD
        int world_rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

        Commxx comm(comm_type::world);
        REQUIRE_NOTHROW(comm.get_rank());
        REQUIRE_NOTHROW(comm.rank());

        REQUIRE(world_rank == comm.get_rank());
        REQUIRE(world_rank == comm.rank());
    }
}

TEST_CASE("is_null")
{
    SECTION("null")
    {
        Commxx comm(comm_type::null);
        REQUIRE(comm.is_null());
    }

    SECTION("world")
    {
        Commxx comm(comm_type::world);
        REQUIRE(!comm.is_null());
    }
}

TEST_CASE("is_root")
{
    SECTION("null")
    {
        Commxx comm(comm_type::null);
        REQUIRE(comm.is_root());
    }

    SECTION("world")
    {
        Commxx comm(comm_type::world);
        REQUIRE(comm.is_root());
    }
}

TEST_CASE("get_parent")
{
    SECTION("null")
    {
        Commxx comm(comm_type::null);
        REQUIRE_THROWS(comm.parent());
    }

    SECTION("world")
    {
        Commxx comm(comm_type::world);
        REQUIRE_THROWS(comm.parent());
    }
}

TEST_CASE("has_this_rank")
{
    SECTION("null")
    {
        Commxx comm(comm_type::null);
        REQUIRE(!comm.has_this_rank());
    }

    SECTION("world")
    {
        Commxx comm(comm_type::world);
        REQUIRE(comm.has_this_rank());
    }
}

TEST_CASE("dup")
{
    SECTION("null")
    {
        // cannot dup on a null communicator
        auto comm = Commxx::Null;
        REQUIRE_THROWS(comm.dup());
    }

    SECTION("world non-shared_ptr")
    {
        // dup on a non shared_ptr is not allowed
        auto comm = Commxx::World;
        REQUIRE_THROWS(comm.dup());
    }

    SECTION("world")
    {
        auto comm1_ptr = std::make_shared<Commxx>(comm_type::world);
        REQUIRE_NOTHROW(comm1_ptr->dup());

        auto comm2 = comm1_ptr->dup();
        auto const& comm1 = *comm1_ptr;

        // duplicated communicator is not the same as the original one
        CHECK(comm1 != comm2);

        // but the size and rank of each processor should be the same
        CHECK(comm1.size() == comm2.size());
        CHECK(comm1.rank() == comm2.rank());

        // compare the groups of two communicator
        MPI_Group g1;
        MPI_Group g2;

        MPI_Comm_group(comm1, &g1);
        MPI_Comm_group(comm2, &g2);

        int result;
        int res = MPI_Group_compare(g1, g2, &result);

        // the groups should be MPI_IDENT, meaning the member and
        // ordering are exactly the same
        REQUIRE(res == MPI_SUCCESS);
        CHECK(result == MPI_IDENT);

        MPI_Group_free(&g1);
        MPI_Group_free(&g2);
    }
}

TEST_CASE("split with color only")
{
    SECTION("null")
    {
        // cannot split a null comm
        auto comm = Commxx::Null;
        REQUIRE_THROWS(comm.split(0));
    }

    SECTION("world non-shared_ptr")
    {
        // split on a non shared_ptr is not allowed
        auto comm = Commxx::World;
        REQUIRE_THROWS(comm.split(0));
    }

    SECTION("world")
    {
        auto comm1_ptr = std::make_shared<Commxx>(comm_type::world);
        auto const& comm1 = *comm1_ptr;

        auto world_size = Commxx::world_size();
        auto world_rank = Commxx::world_rank();

        SECTION("no throw")
        {
            REQUIRE_NOTHROW(comm1_ptr->split(0));
        }

        SECTION("same color")
        {
            // this is effective just a dup
            auto comm2 = comm1_ptr->split(0);

            CHECK(comm2.size() == world_size);
            CHECK(comm2.rank() == world_rank);
        }

        SECTION("one color per rank")
        {
            int color = world_rank;
            auto comm2 = comm1_ptr->split(color);

            CHECK(comm2.size() == 1);
            CHECK(comm2.rank() == 0);
        }

        if (world_size == 3)
        {
            SECTION("split color [0 1 1]")
            {
                // workd_rank:       0  1  2
                // color:            0  1  1
                // size after split: 1  2  2
                // rank after split: 0  0  1
                int color[3] = {0, 1, 1};
                int size[3]  = {1, 2, 2};
                int rank[3]  = {0, 0, 1};

                auto comm2 = comm1_ptr->split(color[world_rank]);
                CHECK(comm2.size() == size[world_rank]);
                CHECK(comm2.rank() == rank[world_rank]);
            }
        }
        else if (world_size == 4)
        {
            SECTION("split color [0 0 1 1]")
            {
                // workd_rank:       0  1  2  3
                // color:            0  0  1  1
                // size after split: 2  2  2  2
                // rank after split: 0  1  0  1
                int color[4] = {0, 0, 1, 1};
                int size[4]  = {2, 2, 2, 2};
                int rank[4]  = {0, 1, 0, 1};

                auto comm2 = comm1_ptr->split(color[world_rank]);
                CHECK(comm2.size() == size[world_rank]);
                CHECK(comm2.rank() == rank[world_rank]);
            }

            SECTION("split color [0 1 1 1]")
            {
                // workd_rank:       0  1  2  3
                // color:            0  1  1  1
                // size after split: 1  3  3  3
                // rank after split: 0  0  1  2
                int color[4] = {0, 1, 1, 1};
                int size[4]  = {1, 3, 3, 3};
                int rank[4]  = {0, 0, 1, 2};

                auto comm2 = comm1_ptr->split(color[world_rank]);
                CHECK(comm2.size() == size[world_rank]);
                CHECK(comm2.rank() == rank[world_rank]);
            }
        }
    }
}

TEST_CASE("split with color and key")
{
    SECTION("null")
    {
        // cannot split a null comm
        auto comm = Commxx::Null;
        REQUIRE_THROWS(comm.split(0, 0));
    }

    SECTION("world non-shared_ptr")
    {
        // split on a non shared_ptr is not allowed
        auto comm = Commxx::World;
        REQUIRE_THROWS(comm.split(0, 0));
    }

    SECTION("world")
    {
        auto comm1_ptr = std::make_shared<Commxx>(comm_type::world);
        auto const& comm1 = *comm1_ptr;

        auto world_size = Commxx::world_size();
        auto world_rank = Commxx::world_rank();

        SECTION("no throw")
        {
            REQUIRE_NOTHROW(comm1_ptr->split(0, world_rank));
        }

        SECTION("same color, ascending key")
        {
            // this is effective just a dup
            auto comm2 = comm1_ptr->split(0, world_rank);

            CHECK(comm2.size() == world_size);
            CHECK(comm2.rank() == world_rank);
        }

        SECTION("same color, descending key")
        {
            // this is effective just a dup
            auto comm2 = comm1_ptr->split(0, world_size - world_rank - 1);

            CHECK(comm2.size() == world_size);
            CHECK(comm2.rank() == world_size - world_rank - 1);
        }

        SECTION("one color per rank")
        {
            int color = world_rank;
            auto comm2 = comm1_ptr->split(color, world_rank);

            CHECK(comm2.size() == 1);
            CHECK(comm2.rank() == 0);
        }

        if (world_size == 3)
        {
            SECTION("split color/key [0/0 1/2 1/1]")
            {
                // workd_rank:       0  1  2
                // color:            0  1  1
                // key:              0  2  1
                // size after split: 1  2  2
                // rank after split: 0  1  0
                int color[3] = {0, 1, 1};
                int key[3]   = {0, 2, 1};
                int size[3]  = {1, 2, 2};
                int rank[3]  = {0, 1, 0};

                auto comm2 = comm1_ptr->split(color[world_rank], key[world_rank]);
                CHECK(comm2.size() == size[world_rank]);
                CHECK(comm2.rank() == rank[world_rank]);
            }
        }
        else if (world_size == 4)
        {
            SECTION("split color/key [0/0 0/1 1/2 1/1]")
            {
                // workd_rank:       0  1  2  3
                // color:            0  0  1  1
                // key:              0  1  2  1
                // size after split: 2  2  2  2
                // rank after split: 0  1  1  0
                int color[4] = {0, 0, 1, 1};
                int key[4]   = {0, 1, 2, 1};
                int size[4]  = {2, 2, 2, 2};
                int rank[4]  = {0, 1, 1, 0};

                auto comm2 = comm1_ptr->split(color[world_rank], key[world_rank]);
                CHECK(comm2.size() == size[world_rank]);
                CHECK(comm2.rank() == rank[world_rank]);
            }

            SECTION("split color/key [0/0 1/3 1/2 1/1]")
            {
                // workd_rank:       0  1  2  3
                // color:            0  1  1  1
                // key:              0  3  2  1
                // size after split: 1  3  3  3
                // rank after split: 0  2  1  0
                int color[4] = {0, 1, 1, 1};
                int key[4]   = {0, 3, 2, 1};
                int size[4]  = {1, 3, 3, 3};
                int rank[4]  = {0, 2, 1, 0};

                auto comm2 = comm1_ptr->split(color[world_rank], key[world_rank]);
                CHECK(comm2.size() == size[world_rank]);
                CHECK(comm2.rank() == rank[world_rank]);
            }
        }
    }
}


#if 0
TEST_CASE("construct4")
{
    Commxx_sptr parent_sptr(new Commxx);
    std::vector<int > ranks(1);
    ranks.at(0) = 0;
    bool per_host = true;
    Commxx(parent_sptr, ranks, per_host);
}

TEST_CASE("get_rank")
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int commxx_rank = Commxx().get_rank();
    BOOST_CHECK_EQUAL(rank, commxx_rank);
}

TEST_CASE("get_size")
{
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int commxx_size = Commxx().get_size();
    BOOST_CHECK_EQUAL(size, commxx_size);
}

TEST_CASE("has_this_rank1")
{
    Commxx commxx;
    BOOST_CHECK(commxx.has_this_rank());
}

TEST_CASE("has_this_rank2")
{
    Commxx commxx(true);
    BOOST_CHECK(commxx.has_this_rank());
}

TEST_CASE("has_this_rank3")
{
    Commxx_sptr parent_sptr(new Commxx);
    std::vector<int > ranks(1);
    const int only_rank = 0;
    ranks.at(0) = only_rank;
    Commxx commxx(parent_sptr, ranks);

    BOOST_CHECK_EQUAL(commxx.has_this_rank(), (parent_sptr->get_rank() == only_rank));
}
#endif

#if 0
TEST_CASE("get")
{
    MPI_Comm comm = Commxx().get();
    int result;
    MPI_Comm_compare(MPI_COMM_WORLD, comm, &result);
    BOOST_CHECK(result == MPI_IDENT);
}

TEST_CASE("serialize1")
{
    const std::string serialize_file_name("commxx1.xml");
    {
        Commxx commxx;
        remove_serialization_directory();
        xml_save(commxx, get_serialization_path(serialize_file_name).c_str(),
                true);
    }

    {
        Commxx commxx_resumed;
        xml_load(commxx_resumed,
                get_serialization_path(serialize_file_name).c_str());
    }
}

TEST_CASE("serialize2")
{
    const std::string serialize_file_name("commxx2.xml");
    {
        Commxx commxx(true);
        remove_serialization_directory();
        xml_save(commxx, get_serialization_path(serialize_file_name).c_str(),
                true);
    }

    {
        Commxx commxx_resumed;
        xml_load(commxx_resumed,
                get_serialization_path(serialize_file_name).c_str());
    }
}

TEST_CASE("serialize3")
{
    const std::string serialize_file_name("commxx3.xml");
    {
        Commxx_sptr parent_sptr(new Commxx);
        std::vector<int > ranks(1);
        ranks.at(0) = 0;
        Commxx commxx(parent_sptr, ranks);
        remove_serialization_directory();
        xml_save(commxx, get_serialization_path(serialize_file_name).c_str(),
                true);
    }

    {
        Commxx commxx_resumed;
        xml_load(commxx_resumed,
                get_serialization_path(serialize_file_name).c_str());
    }
}

TEST_CASE("serialize4")
{
    const std::string serialize_file_name("commxx4.xml");
    {
        Commxx_sptr parent_sptr(new Commxx);
        std::vector<int > ranks(1);
        ranks.at(0) = 0;
        bool per_host = true;
        Commxx commxx(parent_sptr, ranks, per_host);
        remove_serialization_directory();
        xml_save(commxx, get_serialization_path(serialize_file_name).c_str(),
                true);
    }

    {
        Commxx commxx_resumed;
        xml_load(commxx_resumed,
                get_serialization_path(serialize_file_name).c_str());
    }
}

TEST_CASE("generate_subcomms_")
{
    Commxx_sptr parent_sptr(new Commxx);
    int world_size = parent_sptr->get_size();
    for (int size = 1; size<5; ++size) {
        Commxxs commxxs(generate_subcomms(parent_sptr,size));
        BOOST_CHECK_EQUAL(commxxs.size(), size);
        int includes_this_rank = 0;
        Commxx_sptr last_commxx_sptr;
        int unique_on_this_rank = 0;
        for (int i = 0; i< size; ++i) {
            if (commxxs.at(i)->has_this_rank()) {
                includes_this_rank += 1;
                if (commxxs.at(i) != last_commxx_sptr) {
                    ++unique_on_this_rank;
                }
                last_commxx_sptr = commxxs.at(i);
            }
        }
        if (world_size >= size) {
            BOOST_CHECK_EQUAL(includes_this_rank, 1);
        }
        BOOST_CHECK_EQUAL(unique_on_this_rank, 1);
    }
}
/*
TEST_CASE("make_optimal_spc_comm_world")
{
    Commxx_sptr parent_sptr(new Commxx);
    int world_size = parent_sptr->get_size();
    int optimal_number=3;
    Commxx_sptr comm_spc=make_optimal_spc_comm(parent_sptr, optimal_number);
    if (comm_spc->has_this_rank()){
       BOOST_CHECK_EQUAL(parent_sptr->get_rank(), comm_spc->get_rank());
    }
    if (parent_sptr->get_rank()>optimal_number-1) {
      BOOST_CHECK_EQUAL(comm_spc->has_this_rank(), false);
    }      
}*/

TEST_CASE("make_optimal_spc_comm_subcomm")
{
    Commxx_sptr parent_sptr(new Commxx);
    int optimal_number=3;
    for (int size = 1; size<5; ++size) {
       Commxxs commxxs(generate_subcomms(parent_sptr,size));
       for (int i = 0; i< size; ++i) {
	  Commxx_sptr comm_spc(make_optimal_spc_comm(commxxs[i]->get_parent_sptr(), optimal_number));
	  if (comm_spc->has_this_rank())
	        BOOST_CHECK_EQUAL(commxxs[i]->get_parent_sptr()->get_rank(), comm_spc->get_rank());
	  
	  if (commxxs[i]->get_parent_sptr()->get_rank()>optimal_number-1) 
	        BOOST_CHECK_EQUAL(comm_spc->has_this_rank(), false);

	  if (commxxs[i]->get_parent_sptr()->get_rank()==0) 
	       BOOST_CHECK_EQUAL(comm_spc->get_rank(),0);	       	    
       }	       
    }      
}

TEST_CASE("make_optimal_spc_comm_subcomm_equally_spread")
{
    Commxx_sptr parent_sptr(new Commxx);
    int world_size = parent_sptr->get_size();
    int size=2;
    for (int optimal_number=2; optimal_number<4; ++optimal_number){ 
       Commxxs commxxs(generate_subcomms(parent_sptr,size));
       if ((world_size%size==0) && ((world_size/size)%optimal_number==0)){
	  for (int i = 0; i< size; ++i) {
	      if (commxxs[i]->has_this_rank()){
		Commxx_sptr comm_spc(make_optimal_spc_comm(commxxs[i], optimal_number, true));		     
		BOOST_CHECK_EQUAL(comm_spc->has_this_rank(), 1);	
		BOOST_CHECK_EQUAL(commxxs[i]->get_rank()% optimal_number, comm_spc->get_rank());	     
	      }
	  }
       }
    }//optimal number     
}
#endif
