#include "mytimer.h"
#include <iostream>

double 
time()
{
  timeb timestruct;
  ftime(&timestruct);
  return timestruct.time + timestruct.millitm/1000.0;
}
  
double timer_t0;
void
reset_timer()
{
  timer_t0 = time();
}

void
nottimer(std::string message)
{
  double t1 = time();
  std::cout << "timer " << message << ": " << t1 - timer_t0 << std::endl;
  timer_t0 = time();
}

void
timer(std::string message) { return; }
