#ifndef HAVE_RANGE_H
class Range
{
public:
    enum keyword {begin, end};
    Range();
    Range(int only);
    Range(keyword only);
    Range(int first, int last, int step = 1);
    Range(keyword first, int last, int step = 1);
    Range(int first, keyword last, int step = 1);
    Range(keyword first, keyword last, int step = 1);

    bool is_unit_length() const;
    int get_first(int min, int max) const;
    int get_last(int min, int max) const;
    int get_step() const;

private:
    int first, last, step;
    keyword first_keyword, last_keyword;
    bool first_is_keyword, last_is_keyword;
    int keyword_to_int(keyword word, int lower, int upper) const;
};

#define HAVE_RANGE_H  true
#endif // HAVE_RANGE_H
