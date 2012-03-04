#include <iostream>
#include "synergia/utils/command_line_arg.h"
#include "synergia/simulation/resume.h"

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
    std::cout << "\n";
    exit(retval);
}

struct Resume_options
{
    std::string directory;
    std::string new_checkpoint_directory;
    int checkpoint_period;
    int max_turns;

    static const int unspecified_int = -1;
    static const char unspecified_str[];

    Resume_options(int argc, char **argv) :
        directory(Propagator::default_checkpoint_dir),
                new_checkpoint_directory(unspecified_str),
                checkpoint_period(unspecified_int), max_turns(unspecified_int)
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
                    } else if (arg.get_lhs() == "--max") {
                        max_turns = arg.extract_value<int > ();
                    } else {
                        std::cout << "Unknown argument " << arg.get_lhs()
                                << std::endl;
                        usage(1);
                    }
                } else {
                    if (argv[i][0] == '-') {
                        std::cout << "Unknown argument " << argv[i]
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
#if 0
    std::cout << "about to load from " << opts.directory << std::endl;
    Resume resume(opts.directory);
    std::cout << "jfa minor success\n";

    if (opts.checkpoint_period != Resume_options::unspecified_int) {
        resume.set_checkpoint_period(opts.checkpoint_period);
    }
    if (opts.new_checkpoint_directory != Resume_options::unspecified_str) {
        resume.set_new_checkpoint_dir(opts.new_checkpoint_directory);
    }
    if (opts.max_turns == Resume_options::unspecified_int) {
        std::cout << "jfa 00\n";
        resume.propagate();
    } else {
        resume.propagate(opts.max_turns);
    }
#endif
    Propagator propagator;
    remove_serialization_directory();
    symlink_serialization_directory(Propagator::default_checkpoint_dir);
    binary_load(
            propagator,
            get_combined_path(Propagator::default_checkpoint_dir, "propagator.bina").c_str());
    unlink_serialization_directory();
    propagator.resume(Propagator::default_checkpoint_dir);

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
