#ifndef DIGITS_H_
#define DIGITS_H_

long int
digits(long int val)
{
    long int retval = 1;
    long int base = 1;
    for (int i = 0; i < 9; ++i) {
        base *= 10;
        if (val > base) {
            retval += 1;
        }
    }
    return retval;
}

#endif /* DIGITS_H_ */
