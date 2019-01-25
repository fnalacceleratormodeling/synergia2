#include <iostream>
#include "synergia/utils/command_line_arg.h"
#include "synergia/simulation/resume.h"

#include "ramp_actions.h"

void
usage(int retval)
{
    std::cout << "usage: synergia-resume [options] [checkpoint directory]\n";
    std::cout << "  options:\n";
    std::cout << "    --help: this message\n";
    std::cout
            << "    --new-dir=<dir>: directory name to use for subsequent checkpointing\n";
    std::cout << "    --period=<period>: period for subsequent checkpointing\n";
    std::cout << "    --max=<num>: maximum number of turns for this run\n";
    std::cout << "    --turns=<num>: total number of turns to run\n";
    std::cout << "    --concurrent-io=<num>: number of processors that simultaneously write checkpoints\n";
    std::cout << "    --verbosity=<num>: verbosity for this run\n";
    exit(retval);
}

struct Resume_options
{
    std::string directory;
    std::string new_checkpoint_directory;
    int checkpoint_period;
    int num_turns;
    int max_turns;
    int concurrent_io;
    int verbosity;

    static const int unspecified_int = -1;
    static const char unspecified_str[];

    Resume_options(int argc, char **argv) :
        directory(Propagator::default_checkpoint_dir),
                new_checkpoint_directory(unspecified_str),
                checkpoint_period(unspecified_int), num_turns(unspecified_int),
                max_turns(unspecified_int),
                concurrent_io(unspecified_int),
                verbosity(unspecified_int)
    {
        for (int i = 1; i < argc; ++i) {
            if (std::string(argv[i]) == "--help") {
                usage(0);
            } else {
                Command_line_arg arg(argv[i]);
                if (arg.is_equal_pair()) {
                    if (arg.get_lhs() == "--new-dir") {
                        new_checkpoint_directory = arg.extract_value<
                                std::string > ();
                    } else if (arg.get_lhs() == "--period") {
                        checkpoint_period = arg.extract_value<int > ();
                    } else if (arg.get_lhs() == "--turns") {
                    	num_turns = arg.extract_value<int > ();
                    } else if (arg.get_lhs() == "--max") {
                        max_turns = arg.extract_value<int > ();
                    } else if (arg.get_lhs() == "--concurrent-io") {
                        concurrent_io = arg.extract_value<int > ();
                     } else if (arg.get_lhs() == "--verbosity") {
                        verbosity = arg.extract_value<int > ();
                    } else {
                        std::cerr << "Unknown argument " << arg.get_lhs()
                                << std::endl;
                        usage(1);
                    }
                } else {
                    if (argv[i][0] == '-') {
                        std::cerr << "Unknown argument " << argv[i]
                                << std::endl;
                        usage(1);
                    }
                    directory = argv[i];
                }
            }
        }
    }
};

const char Resume_options::unspecified_str[] = "";

void
run(Resume_options &opts)
{
    Resume resume(opts.directory);

    if (opts.checkpoint_period != Resume_options::unspecified_int) {
        resume.set_checkpoint_period(opts.checkpoint_period);
    }
    if (opts.new_checkpoint_directory != Resume_options::unspecified_str) {
        resume.set_new_checkpoint_dir(opts.new_checkpoint_directory);
    }
    bool new_num_turns = (opts.num_turns != Resume_options::unspecified_int);
    bool new_max_turns = (opts.max_turns != Resume_options::unspecified_int);
    bool new_verbosity = (opts.verbosity != Resume_options::unspecified_int);
    if (opts.concurrent_io != Resume_options::unspecified_int) {
    	resume.set_concurrent_io(opts.concurrent_io);
    }
    resume.propagate(new_num_turns, opts.num_turns, new_max_turns, opts.max_turns, new_verbosity,
            opts.verbosity);
}

int
main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    Resume_options opts(argc, argv);
    run(opts);
    MPI_Finalize();
    return 0;
}
