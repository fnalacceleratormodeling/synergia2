#include "range.h"

Range::Range()
{
    first_is_keyword = true;
    first_keyword = begin;
    last_is_keyword = true;
    last_keyword = end;
    step = 1;
}

Range::Range(int only)
{
    first_is_keyword = false;
    last_is_keyword = false;
    first = only;
    last = only;
    step = 1;
}

Range::Range(keyword only)
{
    first_is_keyword = true;
    last_is_keyword = true;
    first_keyword = only;
    last_keyword = only;
    step = 1;
}

Range::Range(int first, int last, int step)
{
    first_is_keyword = false;
    this->first = first;
    last_is_keyword = false;
    this->last = last;
    this->step = step;
}

Range::Range(keyword first, int last, int step)
{
    first_is_keyword = true;
    first_keyword = first;
    last_is_keyword = false;
    this->last = last;
    this->step = step;
}

Range::Range(int first, keyword last, int step)
{
    first_is_keyword = false;
    this->first = first;
    last_is_keyword = true;
    last_keyword = last;
    this->step = step;
}

Range::Range(keyword first, keyword last, int step)
{
    first_is_keyword = true;
    first_keyword = first;
    last_is_keyword = true;
    last_keyword = last;
    this->step = step;
}

//~ void
//~ Range::show(int lower, int upper) const
//~ {
//~ if (first_is_keyword) {std::cout << keyword_to_int(first_keyword, lower, upper);}
//~ else {std::cout << first;}
//~ std::cout << ":";
//~ if (last_is_keyword) {std::cout << keyword_to_int(last_keyword, lower, upper);}
//~ else {std::cout << last;}
//~ if (step != 1) {
//~ std::cout << ":" << step;
//~ }
//~ }

bool
Range::is_unit_length() const
{
    if (first_is_keyword && last_is_keyword &&
            (first_keyword == last_keyword)) {
        return true;
    } else if ((! first_is_keyword) && (! last_is_keyword) &&
               (first == last)) {
        return true;
    } else {
        return false;
    }
}

int
Range::get_first(int min, int max) const
{
    if (first_is_keyword) {
        return keyword_to_int(first_keyword,min,max);
    } else {
        return first;
    }
}

int
Range::get_last(int min, int max) const
{
    if (last_is_keyword) {
        return keyword_to_int(last_keyword,min,max);
    } else {
        return last;
    }
}

int 
Range::get_step() const
{
    return step;
}

int
Range::keyword_to_int(keyword word, int lower, int upper) const
{
    if (word == begin) {
        return lower;
    } else {
        return upper;
    }
}
