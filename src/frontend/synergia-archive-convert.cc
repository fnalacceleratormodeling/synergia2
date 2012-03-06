#include <iostream>
#include "synergia/utils/command_line_arg.h"
#include "synergia/lattice/lattice.h"
#include "synergia/simulation/propagator.h"
#include "synergia/utils/multi_array_typedefs.h"


struct Convert_options
{
    std::string input_file;
    std::string input_suffix;
    std::string output_file;
    std::string output_suffix;

    Convert_options(int argc, char **argv)
    {
        bool first = true;
        bool done = false;
        for (int i = 1; i < argc; ++i) {
            if (std::string(argv[i]) == "--help") {
                usage(0);
            } else {
                if (first) {
                    input_file = argv[i];
                    input_suffix = get_suffix(input_file);
                    first = false;
                } else {
                    if (!done) {
                        output_file = argv[i];
                        output_suffix = get_suffix(output_file);
                        done = true;
                    } else {
                        std::cerr << "Unknown argument " << argv[i]
                                << std::endl;
                        usage(1);
                    }
                }
            }
        }
        if (!done) {
            usage(1);
        }
        if (input_suffix == output_suffix) {
            std::cerr
                    << "synergia-archive-convert: will not \"convert\" equal types\n";
            exit(1);
        }
    }

    std::string
    get_suffix(std::string const& filename)
    {
        std::string suffix;
        int idx = filename.rfind('.');
        if (idx == std::string::npos) {
            suffix = "";
        } else {
            suffix = filename.substr(idx + 1);
        }
        if ((suffix != "bina") && (suffix != "txta") && (suffix != "xml")) {
            std::cerr << "Unknown suffix in file name \"" << filename << "\""
                    << std::endl;
            usage(1);
        }
        return suffix;
    }

    void
    usage(int retval)
    {
        std::cout
                << "usage: synergia-archive-convert <archive> <new archive>\n";
        std::cout << "       archive types are inferred from file suffixes:\n";
        std::cout << "           .bina : binary\n";
        std::cout << "           .txta : text\n";
        std::cout << "           .xml  : xml\n";
        std::cout << "  options:\n";
        std::cout << "    --help: this message\n";
        exit(retval);
    }
};

template<typename T, typename A_in, typename A_out>
    void
    real_convert(std::string const& filename_in,
            std::string const& filename_out)
    {
        T object;
        archive_load<T, A_in > (object, filename_in);
        archive_save<T, A_out > (object, filename_out);
    }

template<typename A>
    std::string
    get_object_typename(std::string const& filename)
    {
        std::ifstream input_stream(filename.c_str());
        if (!input_stream.good()) {
            std::cerr << "synergia-archive-convert: unable to open "
                    << filename << std::endl;
            exit(1);
        }
        A input_archive(input_stream);
        int num_objects;
        input_archive >> BOOST_SERIALIZATION_NVP(num_objects);
        std::string object_typename;
        input_archive >> BOOST_SERIALIZATION_NVP(object_typename);
        input_stream.close();

        return object_typename;
    }

template<typename A_in, typename A_out>
    void
    convert(std::string const& filename_in, std::string const& filename_out)
    {
        std::string object_typename(get_object_typename<A_in > (filename_in));

        if (object_typename == typeid(Lattice).name()) real_convert<Lattice,
                A_in, A_out > (filename_in, filename_out);
        else if (object_typename == typeid(Propagator).name()) real_convert<
                Propagator, A_in, A_out > (filename_in, filename_out);
        else if (object_typename == typeid(Propagator::State).name()) real_convert<
                Propagator::State, A_in, A_out > (filename_in, filename_out);
        else if (object_typename == typeid(MArray1d).name()) real_convert<
                MArray1d, A_in, A_out > (filename_in, filename_out);
        else if (object_typename == typeid(MArray2d).name()) real_convert<
                MArray2d, A_in, A_out > (filename_in, filename_out);
        else {
            std::cerr
                    << "synergia-archive-convert: do not know how to convert objects of type \""
                    << object_typename << "\"\n";
            exit(1);
        }
    }

void
run(Convert_options &opts)
{
    if (opts.input_suffix == "bina") {
        if (opts.output_suffix == "txta") {
            convert<boost::archive::binary_iarchive,
                    boost::archive::text_oarchive > (opts.input_file,
                    opts.output_file);
        } else {
            if (opts.output_suffix == "xml") {
                convert<boost::archive::binary_iarchive,
                        boost::archive::xml_oarchive > (opts.input_file,
                        opts.output_file);
            }
        }
    }
    if (opts.input_suffix == "txta") {
        if (opts.output_suffix == "bina") {
            convert<boost::archive::text_iarchive,
                    boost::archive::binary_oarchive > (opts.input_file,
                    opts.output_file);
        } else {
            if (opts.output_suffix == "xml") {
                convert<boost::archive::text_iarchive,
                        boost::archive::xml_oarchive > (opts.input_file,
                        opts.output_file);
            }
        }
    }
    if (opts.input_suffix == "xml") {
        if (opts.output_suffix == "txta") {
            convert<boost::archive::xml_iarchive, boost::archive::text_oarchive > (
                    opts.input_file, opts.output_file);
        } else {
            if (opts.output_suffix == "bina") {
                convert<boost::archive::xml_iarchive,
                        boost::archive::binary_oarchive > (opts.input_file,
                        opts.output_file);
            }
        }
    }
}

int
main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    Convert_options opts(argc, argv);
    run(opts);
    MPI_Finalize();
    return 0;
}
