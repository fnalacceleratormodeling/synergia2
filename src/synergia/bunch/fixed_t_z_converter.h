#ifndef FIXED_T_Z_CONVERTER_H_

#include "synergia/utils/serialization.h"
class Bunch;

/// Fixed_t_z_converter is a virtual base class for converters
/// between the bunch fixed-z representation and the bunch
/// fixed-t representation

/// There are three converters one can choose:
///   1. Fixed_t_z_alex: the most general, considers the horizontal motion (beta_x_i and
///       beta_y_i  or particles not zeros) and the reference time and z are the ones for the
///       reference particle.
///   2.  Fixed_t_z_synergia20 corresponds to the old version of synergia, obtained from
///       Fixed_t_z_alex by making beta_x=0, beta_y=0, (1-beta*beta_z_i) =gamma^-2. Notice
///       that this approximation condisers different velocities along z direction for particles
///       i.e. in general  beta_z_i != beta.
///   3.  Fixed_t_z_zeroth: zeroth order approximation, obtained from  Fixed_t_z_alex
///       by making beta_x_i=beta_y_i=0 and beta_z_i = beta.

class Fixed_t_z_converter
{
public:
    /// Convert from the fixed-z state in the accelerator frame to the fixed-t state in the accelerator frame.
    virtual void
    from_z_lab_to_t_lab(Bunch &bunch) = 0;
    /// Convert from the fixed-t state in the accelerator frame to the fixed-z state in the accelerator frame.
    virtual void
    from_t_lab_to_z_lab(Bunch &bunch) = 0;
    /// Convert from the fixed-z state in the accelerator frame to the fixed-t state in the beam frame.
    virtual void
    from_z_lab_to_t_bunch(Bunch &bunch) = 0;
    /// Convert from the fixed-t state in the beam frame to the fixed-z state in the accelerator frame.
    virtual void
    from_t_bunch_to_z_lab(Bunch &bunch) = 0;

    virtual void
    from_t_lab_to_t_bunch(Bunch &bunch) = 0;

    virtual void
    from_t_bunch_to_t_lab(Bunch &bunch)=0;

//
//     virtual void
//     fixed_z_lab_to_z_bunch(Bunch &bunch) = 0;
//


    /// Convert from the fixed-t state to the fixed-z state.
    // virtual void
    // fixed_t_to_fixed_z(Bunch &bunch) = 0;

    /// Convert from the fixed-z state to the fixed-t state.
    //  virtual void
    //  fixed_z_to_fixed_t(Bunch &bunch) = 0;

    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version)
        {
        }

    virtual
    ~Fixed_t_z_converter()
    {
    }
    ;
};
BOOST_CLASS_EXPORT_KEY(Fixed_t_z_converter)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(Fixed_t_z_converter);

/// Fixed_t_z_zeroth implements a fixed-t-fixed-z converter using
/// the simplest approximation: longitudinal coordinates are transformed,
/// but transverse coordinates are unaffected.
class Fixed_t_z_zeroth : public Fixed_t_z_converter
{
public:
    /// Convert from the fixed-z state in the accelerator frame to the fixed-t state in the accelerator frame.
    void
    from_z_lab_to_t_lab(Bunch &bunch);
    /// Convert from the fixed-t state in the accelerator frameto the fixed-z state in the accelerator frame.
    void
    from_t_lab_to_z_lab(Bunch &bunch);
    /// Convert from the fixed-z state in the accelerator frame to the fixed-t state in the beam frame.
    void
    from_z_lab_to_t_bunch(Bunch &bunch);
    /// Convert from the fixed-t state in the beam frame to the fixed-z state in the accelerator frame.
    void
    from_t_bunch_to_z_lab(Bunch &bunch);

    void
    from_t_lab_to_t_bunch(Bunch &bunch);

    void
    from_t_bunch_to_t_lab(Bunch &bunch);

    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version)
        {
            ar &
            BOOST_SERIALIZATION_BASE_OBJECT_NVP(Fixed_t_z_converter);
        }
};
BOOST_CLASS_EXPORT_KEY(Fixed_t_z_zeroth)

/// Fixed_t_z_zeroth implements a fixed-t-fixed-z converter using
/// the ballistic approximation: longitudinal coordinates are transformed,
/// then transverse coordinates are transformed using the ballistic
/// approximation, i.e., as though the particles were traveling independently
/// in free space.
/// UNIMPLEMENTED
class Fixed_t_z_ballistic : public Fixed_t_z_converter
{
public:
    /// Convert from the fixed-z state in the accelerator frame to the fixed-t state in the accelerator frame.
    void
    from_z_lab_to_t_lab(Bunch &bunch)
    {
    }
    ;
    /// Convert from the fixed-t state in the accelerator frame to the fixed-z state in the accelerator frame.
    void
    from_t_lab_to_z_lab(Bunch &bunch)
    {
    }
    ;
    /// Convert from the fixed-z state in the accelerator frame to the fixed-t state in the beam frame.
    void
    from_z_lab_to_t_bunch(Bunch &bunch);
    /// Convert from the fixed-t state in the beam frame to the fixed-z state in the accelerator frame.
    void
    from_t_bunch_to_z_lab(Bunch &bunch);
    void
    from_t_lab_to_t_bunch(Bunch &bunch){};
    void
    from_t_bunch_to_t_lab(Bunch &bunch){};

    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version)
        {
            ar &
            BOOST_SERIALIZATION_BASE_OBJECT_NVP(Fixed_t_z_converter);
        }
};
BOOST_CLASS_EXPORT_KEY(Fixed_t_z_ballistic)

class Fixed_t_z_alex : public Fixed_t_z_converter
{
public:
    /// Convert from the fixed-z state in the accelerator frame to the fixed-t state in the accelerator frame.
    void
    from_z_lab_to_t_lab(Bunch &bunch);
    /// Convert from the fixed-t state in the accelerator frame to the fixed-z state in the accelerator frame.
    void
    from_t_lab_to_z_lab(Bunch &bunch);
    /// Convert from the fixed-z state in the accelerator frame to the fixed-t state in the beam frame.
    void
    from_z_lab_to_t_bunch(Bunch &bunch);
    /// Convert from the fixed-t state in the beam frame to the fixed-z state in the accelerator frame.
    void
    from_t_bunch_to_z_lab(Bunch &bunch);
     void
     from_t_lab_to_t_bunch(Bunch &bunch);

     void
     from_t_bunch_to_t_lab(Bunch &bunch);

    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version)
        {
            ar &
            BOOST_SERIALIZATION_BASE_OBJECT_NVP(Fixed_t_z_converter);
        }
};
BOOST_CLASS_EXPORT_KEY(Fixed_t_z_alex)

/// transformation as in the old synergia....
class Fixed_t_z_synergia20 : public Fixed_t_z_converter
{
public:
    /// Convert from the fixed-z state in the accelerator frame to the fixed-t state in the accelerator frame.
    void
    from_z_lab_to_t_lab(Bunch &bunch);
    /// Convert from the fixed-t state in the accelerator frame to the fixed-z state in the accelerator frame.
    void
    from_t_lab_to_z_lab(Bunch &bunch);
    /// Convert from the fixed-z state in the accelerator frame to the fixed-t state in the beam frame.
    void
    from_z_lab_to_t_bunch(Bunch &bunch);
    /// Convert from the fixed-t state in the beam frame to the fixed-z state in the accelerator frame.
    void
    from_t_bunch_to_z_lab(Bunch &bunch);
    void
    from_t_lab_to_t_bunch(Bunch &bunch);
    void
    from_t_bunch_to_t_lab(Bunch &bunch);

    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version)
        {
            ar &
            BOOST_SERIALIZATION_BASE_OBJECT_NVP(Fixed_t_z_converter);
        }
};
BOOST_CLASS_EXPORT_KEY(Fixed_t_z_synergia20)

#define FIXED_T_Z_CONVERTER_H_

#endif /* FIXED_T_Z_CONVERTER_H_ */
