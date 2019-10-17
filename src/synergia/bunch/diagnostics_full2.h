#ifndef DIAGNOSTICS_FULL2_H_
#define DIAGNOSTICS_FULL2_H_

#include "synergia/bunch/diagnostics.h"

/// Diagnostics_full2 provides the full set of statistical
/// quantities to be calculated for a Bunch up to the second moments.
class Diagnostics_full2 : public Diagnostics
{
public:

    constexpr static const char* diag_type = "diagnostics_full2";
    constexpr static const bool  diag_write_serial = true;


private:

    bool first_write = true;

    double s;
    double s_n;

    int repetition;
    int num_particles;

    double real_num_particles;
    double pz;

    karray1d mean;
    karray1d std;
    karray1d min;
    karray1d max;
    karray2d_row mom2;
    karray2d_row corr;

    double emitx, emity, emitz, emitxy, emitxyz;

private:

    void do_update(Bunch const& bunch) override;
    void do_write (Bunch const& bunch) override;

public:

    /// Create a Diagnostics_full2 object
    /// @param bunch the Bunch
    /// @param filename filename for output
    /// @param local_dir local directory to use for temporary scratch
    Diagnostics_full2(
            std::string const& filename, 
            std::string const& local_dir = "" );
};
#endif /* DIAGNOSTICS_FULL2_H_ */
