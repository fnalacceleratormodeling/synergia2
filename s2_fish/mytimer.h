#ifndef MYTIMER_H
#include <string>
#include <sys/timeb.h>
double time();
void reset_timer();
void timer(std::string message);
#define MYTIMER_H
#endif
