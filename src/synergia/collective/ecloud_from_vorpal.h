#ifndef ECLOUD_FROM_VORPAL_H
#define  ECLOUD_FROM_VORPAL_H
#include "synergia/simulation/operator.h"
#include "synergia/simulation/step.h"
#include "synergia/bunch/bunch.h"
#include "synergia/utils/commxx.h"
#include "synergia/collective/ecloud_field_vorpal2d.h"

class Ecloud_from_vorpal : public Collective_operator
{
private:
   std::string file_name_archive; // the UNix path name for where the archive file that contains the transverse fields representation.  
   Commxx_sptr comm_sptr; // A single commnicator. Used only once at instantiation, but kept here for possible futur use.
   ECloudEFieldVORPAL2D e_field; // E-Field represented as a series of Chebychev coefficicent Both X and Y component included. 
   std::string field_name; // Just a name, good practice to name things.
   std::vector<std::string> subjectedDevices; // a list of CHEF device name where this e-cloud is present. OBSOLETE, Done from above, I am told. 
   
public: 
    // Constructor takes a comm pointer, the file name where the Field is, and a a device Name. 
    // if none provided, assume the entire ring is subjected to the e-cloud. 
    Ecloud_from_vorpal(Commxx_sptr comm_sptr, const std::string &file_name_archive, const std::string aDeviceName=std::string("all")); 
    Ecloud_from_vorpal(); // Empty collective operator 
    
    virtual void
      apply(Bunch & bunch, double time_step, Step & step, int verbosity,
            Logger & logger);
    virtual
    ~Ecloud_from_vorpal();
    
   
    inline std::string get_field_name() { return field_name;}
    inline std::string get_file_name_archive() { return field_name;}
    inline size_t get_number_devices() { return subjectedDevices.size();}
    inline std::string get_device_name(size_t k=0) {if (subjectedDevices.size() == 0) return std::string("None");
                                                    return subjectedDevices[k]; }
    inline void add_device(const std::string &device) { subjectedDevices.push_back(device);}						    					   
};
#endif // ECLOUD_FROM_VORPAL_H
BOOST_CLASS_EXPORT_KEY(Ecloud_from_vorpal)

typedef boost::shared_ptr<Ecloud_from_vorpal > Ecloud_from_vorpal_sptr; // syndoc:include
