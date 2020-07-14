#ifndef PERIOD_H_
#define PERIOD_H_


class Bunch;

void apply_longitudinal_bucket_barrier(Bunch& bunch, double length);
void apply_longitudinal_periodicity(Bunch& bunch, double length); 
void apply_zcut(Bunch& bunch, double length);

#endif /* PERIOD_H_ */
