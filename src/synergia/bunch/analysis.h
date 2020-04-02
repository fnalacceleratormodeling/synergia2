#ifndef ANALYSIS_H_
#define ANALYSIS_H_

#include <string>
#include <boost/shared_ptr.hpp>
#include "synergia/utils/multi_array_typedefs.h"

#include "synergia/bunch/bunch.h"


/// Read back the files produced by Diagnositics_particles 
///  and obtain some variable of interest. 
/// Currently implemented: just the Betatron Tune, based on the so-called "particle" files.  
///

class Analysis
{

private:
    size_t maxNumberOfTurn; // The number of turn performed in the simulation. User intput
    std::string tokenname; // Token, not the full file name  
    std::string filename; // the file name currently opened, volatile. 
     
    size_t numTurns;
    size_t numParticles1rstBunch;
    bool gotHorizontalTunes;
    bool gotVerticalTunes;
    size_t numXTunesFound;
    MArray1d XTunes;
    size_t numYTunesFound;
    MArray1d YTunes;
    size_t minimumRequiredTurnNum;  // The minimum number of turns to get a meaningfull tune 
    MArray3d coords; // X and Y coordinates, for both transverse planes
                     // [TurnNumber][ParticlesNumber][2], X or Y

public:
    explicit Analysis(std::string const& filename, size_t maxTurn=0);

    double get_betatron_averaged_tune(bool isHorizontal);
    double get_betatron_tune_spread(bool isHorizontal);
    std::vector<double> get_betatron_tunes(bool isHorizontal);
    std::vector<double> get_XCoords_forTunes(size_t selectedParticle) const ;
    std::vector<double> get_YCoords_forTunes(size_t selectedParticle) const ;
    size_t get_num_betatron_tune_found(bool isHorizontal) const {
       if (isHorizontal) return numXTunesFound; else return numYTunesFound; return 0;
    } 
 
    void set_minimum_required_turn_Num(size_t n) {minimumRequiredTurnNum=n;}
    size_t get_minimum_required_turn_Num() const {return minimumRequiredTurnNum;}
    
    int get_num_turns() const {return numTurns;}
    int get_num_part_first_bunch() const {return numParticles1rstBunch;}
    std::vector<double> get_transverse_action_for_particle(bool isH,  size_t selectedParticle,
                                                           double alpha, double beta) const;
    std::vector<double> get_transverse_action_for_bunch(bool isH,  size_t turnNumber, 
                                                        double alpha, double beta) const;
      
private:
    void uploadCoords();
    void compute_betatron_tunes(bool isHorizontal);      

};

#endif
