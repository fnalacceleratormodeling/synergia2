#ifndef FIXED_T_Z_CONVERTER_H_


class Bunch;

/// Fixed_t_z_converter is a virtual base class for converters
/// between the bunch fixed-z representation and the bunch
/// fixed-t representation
class Fixed_t_z_converter
{
public:
 
    /// Convert from the fixed-z state in the accelerator frame to the fixed-t state in the accelerator frame.
    virtual void
    from_zacc_to_tacc(Bunch &bunch) = 0;
    /// Convert from the fixed-t state in the accelerator frameto the fixed-z state in the accelerator frame.
    virtual void
    from_tacc_to_zacc(Bunch &bunch) = 0;
    /// Convert from the fixed-z state in the accelerator frame to the fixed-t state in the beam frame.
    virtual void
    from_zacc_to_tbeam(Bunch &bunch) = 0;
     /// Convert from the fixed-t state in the beam frame to the fixed-z state in the accelerator frame.
    virtual void
    from_tbeam_to_zacc(Bunch &bunch) = 0;
    
    
///    not implemented yet
//     virtual void
//     fixed_zacc_to_zbeam(Bunch &bunch) = 0;   
//     
//     virtual void
//     fixed_zacc_to_zbeam(Bunch &bunch) = 0; 
//     
    
    
    /// Convert from the fixed-t state to the fixed-z state.
   // virtual void
   // fixed_t_to_fixed_z(Bunch &bunch) = 0;

    /// Convert from the fixed-z state to the fixed-t state.
  //  virtual void
  //  fixed_z_to_fixed_t(Bunch &bunch) = 0;

    virtual
    ~Fixed_t_z_converter()
    {
    }
    ;
};

/// Fixed_t_z_zeroth implements a fixed-t-fixed-z converter using
/// the simplest approximation: longitudinal coordinates are transformed,
/// but transverse coordinates are unaffected.
class Fixed_t_z_zeroth : public Fixed_t_z_converter
{
public:
     /// Convert from the fixed-z state in the accelerator frame to the fixed-t state in the accelerator frame.
    void
    from_zacc_to_tacc(Bunch &bunch){};
    /// Convert from the fixed-t state in the accelerator frameto the fixed-z state in the accelerator frame.
    void
    from_tacc_to_zacc(Bunch &bunch){};
    /// Convert from the fixed-z state in the accelerator frame to the fixed-t state in the beam frame.
    void
    from_zacc_to_tbeam(Bunch &bunch);
    /// Convert from the fixed-t state in the beam frame to the fixed-z state in the accelerator frame.
    void
    from_tbeam_to_zacc(Bunch &bunch);
    
};

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
    from_zacc_to_tacc(Bunch &bunch){};
    /// Convert from the fixed-t state in the accelerator frameto the fixed-z state in the accelerator frame.
    void
    from_tacc_to_zacc(Bunch &bunch){};
    /// Convert from the fixed-z state in the accelerator frame to the fixed-t state in the beam frame.
    void
    from_zacc_to_tbeam(Bunch &bunch);
    /// Convert from the fixed-t state in the beam frame to the fixed-z state in the accelerator frame.
    void
    from_tbeam_to_zacc(Bunch &bunch);
    
};

class Fixed_t_z_alex : public Fixed_t_z_converter
{
public:  
    /// Convert from the fixed-z state in the accelerator frame to the fixed-t state in the accelerator frame.
    void
    from_zacc_to_tacc(Bunch &bunch);
    /// Convert from the fixed-t state in the accelerator frameto the fixed-z state in the accelerator frame.
    void
    from_tacc_to_zacc(Bunch &bunch);
    /// Convert from the fixed-z state in the accelerator frame to the fixed-t state in the beam frame.
    void
    from_zacc_to_tbeam(Bunch &bunch);
    /// Convert from the fixed-t state in the beam frame to the fixed-z state in the accelerator frame.
    void
    from_tbeam_to_zacc(Bunch &bunch);
};


/// transformation as in the old synergia....
class Fixed_t_z_synergia20 : public Fixed_t_z_converter
{
public:
   /// Convert from the fixed-z state in the accelerator frame to the fixed-t state in the accelerator frame.
    void
    from_zacc_to_tacc(Bunch &bunch);
    /// Convert from the fixed-t state in the accelerator frameto the fixed-z state in the accelerator frame.
    void
    from_tacc_to_zacc(Bunch &bunch);
    /// Convert from the fixed-z state in the accelerator frame to the fixed-t state in the beam frame.
    void
    from_zacc_to_tbeam(Bunch &bunch);
    /// Convert from the fixed-t state in the beam frame to the fixed-z state in the accelerator frame.
    void
    from_tbeam_to_zacc(Bunch &bunch);

    
};


#define FIXED_T_Z_CONVERTER_H_


#endif /* FIXED_T_Z_CONVERTER_H_ */
