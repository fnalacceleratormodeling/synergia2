#include <iostream>
#include <stdexcept>
#include "synergia/simulation/propagator.h"
#include "synergia/simulation/resume.h"
#include "booster_options.h"

void
run()
{

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);  
  
  
  
  std::cout<<"run start on rank: "<<rank<<std::endl; 
  
  Options_resume opts;
  opts.print();
  
  Resume resume(Propagator::default_checkpoint_dir);

  resume.set_checkpoint_period(opts.checkpointperiod);
  resume.set_concurrent_io(opts.concurrentio);
  resume.set_new_checkpoint_dir(opts.directory);
  resume.propagate(opts.new_num_turns, opts.num_turns, opts.new_maxturns, 
                   opts.maxturns, opts.new_verbosity, opts.verbosity);
 
}

int
main(int argc, char **argv)
{
   MPI_Init(&argc, &argv);
   run();
   MPI_Finalize();
   return 0;
}
