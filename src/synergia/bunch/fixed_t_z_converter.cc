#include "synergia/utils/multi_array_typedefs.h"
#include "synergia/foundation/reference_particle.h"
#include "fixed_t_z_converter.h"
#include "bunch.h"
#include "period.h"

#include <iostream>
#include <cmath>

//*********************Fixed_t_z_zeroth*****************************************

void
Fixed_t_z_zeroth::from_z_lab_to_t_bunch(Bunch &bunch)
{
    double gamma = bunch.get_reference_particle().get_gamma();
    double beta = bunch.get_reference_particle().get_beta();
    double m = bunch.get_mass();
    double p_ref = bunch.get_reference_particle().get_momentum();
    MArray2d_ref particles = bunch.get_local_particles();
    for (int part = 0; part < bunch.get_local_num(); ++part) {
        // z in beam rest frame
        particles[part][Bunch::z] = -gamma * beta * particles[part][Bunch::cdt];

        // total momentum in accelerator frame
        double p = p_ref + particles[part][Bunch::dpop] * p_ref;
        // E/c in accelerator frame
        double Eoc = std::sqrt(p * p + m * m);
        // p_{x,y,z} in accelerator frame
        double px = particles[part][Bunch::xp] * p_ref;
        double py = particles[part][Bunch::yp] * p_ref;
        double pz2 = p * p - px * px - py * py;
        if (pz2 < 0.0) {
            throw std::runtime_error(
                    "Fixed_t_z_zeroth::fixed_z_to_fixed_t: particle has negative pz^2");
        }
        double pz = std::sqrt(pz2);
        // zp = pz/p_{ref}^{total}
        particles[part][Bunch::zp] = gamma * (pz - beta * Eoc) / p_ref;

        // n.b. in the zeroth approximation, the transformation from
        //      t' = gamma cdt to t' = 0
        //      is a no-op.
    }
}

void
Fixed_t_z_zeroth::from_t_bunch_to_z_lab(Bunch &bunch)
{
    double gamma = bunch.get_reference_particle().get_gamma();
    double beta = bunch.get_reference_particle().get_beta();
    double m = bunch.get_mass();
    double p_ref = bunch.get_reference_particle().get_momentum();
    MArray2d_ref particles = bunch.get_local_particles();
    for (int part = 0; part < bunch.get_local_num(); ++part) {
        // ct in accelerator frame
        particles[part][Bunch::cdt] = -particles[part][Bunch::z]
                / (gamma * beta);

        // p'_{x,y,z} in beam frame
        double pxp = particles[part][Bunch::xp] * p_ref;
        double pyp = particles[part][Bunch::yp] * p_ref;
        double pzp = particles[part][Bunch::zp] * p_ref;
        double p_perp2 = pxp * pxp + pyp * pyp;
        // E'/c in beam frame
        double Epoc = std::sqrt(p_perp2 + pzp * pzp + m * m);
        double pz = gamma * (pzp + beta * Epoc);
        // dpop = (p - p_ref)/p_ref
        double p = std::sqrt(p_perp2 + pz * pz);
        particles[part][Bunch::dpop] = (p - p_ref) / p_ref;
    }
}


void
Fixed_t_z_zeroth::from_z_lab_to_t_lab(Bunch &bunch)
{
    double m = bunch.get_mass();
    double beta = bunch.get_reference_particle().get_beta();
    double p_ref = bunch.get_reference_particle().get_momentum();
    MArray2d_ref particles = bunch.get_local_particles();
    for (int part = 0; part < bunch.get_local_num(); ++part) {
          // total momentum in accelerator frame
        double p = p_ref + particles[part][Bunch::dpop] * p_ref;
         // p_{x,y,z} in accelerator frame
        double px = particles[part][Bunch::xp] * p_ref;
        double py = particles[part][Bunch::yp] * p_ref;
        double pz2 = p * p - px * px - py * py;
        if (pz2 < 0.0) {
            throw std::runtime_error(
                    "Fixed_t_z_zeroth::fixed_z_to_fixed_t: particle has negative pz^2");
        }
        double pz = std::sqrt(pz2);
        particles[part][Bunch::z] = - particles[part][Bunch::cdt]*beta;

        // zp = pz/p_{ref}^{total}
         particles[part][Bunch::zp] =pz/p_ref;

    }

}

void
Fixed_t_z_zeroth::from_t_lab_to_z_lab(Bunch &bunch)
{

    double m = bunch.get_mass();
    double beta = bunch.get_reference_particle().get_beta();
    double p_ref = bunch.get_reference_particle().get_momentum();
    MArray2d_ref particles = bunch.get_local_particles();

    for (int part = 0; part < bunch.get_local_num(); ++part) {
          // p'_{x,y,z} in beam rest frame
        double px = particles[part][Bunch::xp]* p_ref;
        double py = particles[part][Bunch::yp]* p_ref;
        double pz = particles[part][Bunch::zp]* p_ref;
        double p = std::sqrt(px * px + py * py + pz * pz);
        particles[part][Bunch::dpop] = (p - p_ref) / p_ref;
        particles[part][Bunch::cdt] = - particles[part][Bunch::z]/beta;

    }
}

void
Fixed_t_z_zeroth::from_t_lab_to_t_bunch(Bunch &bunch)
{
     double gamma=bunch.get_reference_particle().get_gamma();
     double beta = bunch.get_reference_particle().get_beta();
     double m = bunch.get_mass();
     double p_ref = bunch.get_reference_particle().get_momentum();
     MArray2d_ref particles = bunch.get_local_particles();
     for (int part = 0; part < bunch.get_local_num(); ++part) {
         double px  = particles[part][Bunch::xp] * p_ref;
         double py  = particles[part][Bunch::yp] * p_ref;
         double pz  = particles[part][Bunch::zp]  * p_ref;
         double p=std::sqrt(px * px + py * py + pz * pz);
         double Eoc   = std::sqrt(p * p + m * m);

         particles[part][Bunch::z] =  gamma * particles[part][Bunch::z];
         particles[part][Bunch::zp] = gamma*(pz - beta * Eoc)/p_ref;
    }
}

void
Fixed_t_z_zeroth::from_t_bunch_to_t_lab(Bunch &bunch)
{
    double gamma = bunch.get_reference_particle().get_gamma();
    double beta = bunch.get_reference_particle().get_beta();
    double m = bunch.get_mass();
    double p_ref = bunch.get_reference_particle().get_momentum();
    MArray2d_ref particles = bunch.get_local_particles();
    for (int part = 0; part < bunch.get_local_num(); ++part) {
          // p'_{x,y,z} in beam rest frame
        double pxp = particles[part][Bunch::xp]* p_ref;
        double pyp = particles[part][Bunch::yp]* p_ref;
        double pzp = particles[part][Bunch::zp]* p_ref;

        // E'/c in beam rest frame
        double Epoc = std::sqrt(pxp*pxp + pyp*pyp + pzp*pzp + m * m);
        double pz = gamma * (pzp + beta * Epoc);
        particles[part][Bunch::z] = particles[part][Bunch::z]/gamma;
        particles[part][Bunch::dpop] = pz / p_ref;
    }
}
BOOST_CLASS_EXPORT_IMPLEMENT(Fixed_t_z_zeroth)
//*********************Fixed_t_z_zeroth end*****************************************

void
Fixed_t_z_ballistic::from_t_bunch_to_z_lab(Bunch &bunch)
{
    std::cout << "stub: ballistic fixed_t_to_fixed_z\n";
}

void
Fixed_t_z_ballistic::from_z_lab_to_t_bunch(Bunch &bunch)
{
    std::cout << "stub: ballistic fixed_z_to_fixed_t\n";
}

BOOST_CLASS_EXPORT_IMPLEMENT(Fixed_t_z_ballistic)

//*********************Fixed_t_z_alex*********************************************

void
Fixed_t_z_alex::from_z_lab_to_t_bunch(Bunch &bunch)
{
    double gamma=bunch.get_reference_particle().get_gamma();
    double gamma_inv = 1.0/gamma;
    double beta = bunch.get_reference_particle().get_beta();
    double m = bunch.get_mass();
    double p_ref = bunch.get_reference_particle().get_momentum();
    MArray2d_ref particles = bunch.get_local_particles();
    for (int part = 0; part < bunch.get_local_num(); ++part) {
         // p_{x,y} in accelerator frame
        double px = particles[part][Bunch::xp] * p_ref;
        double py = particles[part][Bunch::yp] * p_ref;
        // total momentum in accelerator frame
        double p = p_ref + particles[part][Bunch::dpop] * p_ref;
        double pz2 = p * p - px * px - py * py;
        if (pz2 < 0.0) {
            throw std::runtime_error(
                    "Fixed_t_z_zeroth::fixed_z_to_fixed_t: particle has negative pz^2");
        }
        double pz = std::sqrt(pz2);
         // E/c in accelerator frame
        double Eoc = std::sqrt(p * p + m * m);

        double  betaz_part=pz/Eoc;
        double  betax_part=px/Eoc;
        double  betay_part=py/Eoc;
        double  cdt=particles[part][Bunch::cdt];
         // x,y,z in beam rest frame
        particles[part][Bunch::z] = - cdt*gamma_inv*betaz_part/(1.0-beta*betaz_part);
        particles[part][Bunch::x] = particles[part][Bunch::x] -  cdt*betax_part/(1.0-beta*betaz_part);
        particles[part][Bunch::y] = particles[part][Bunch::y] -  cdt*betay_part/(1.0-beta*betaz_part);

        // zp = pz/p_{ref}^{total}, pz in the beam frame,  p_ref in the accelerator frame
       particles[part][Bunch::zp] = gamma*(pz - beta * Eoc)/p_ref;
    }
}

void
Fixed_t_z_alex::from_t_bunch_to_z_lab(Bunch &bunch)
{
    double gamma = bunch.get_reference_particle().get_gamma();
    double beta = bunch.get_reference_particle().get_beta();
    double m = bunch.get_mass();
    double p_ref = bunch.get_reference_particle().get_momentum();
    MArray2d_ref particles = bunch.get_local_particles();
    for (int part = 0; part < bunch.get_local_num(); ++part) {
          // p'_{x,y,z} in beam rest frame
        double pxp = particles[part][Bunch::xp]* p_ref;
        double pyp = particles[part][Bunch::yp]* p_ref;
        double pzp = particles[part][Bunch::zp]* p_ref;

        // E'/c in beam rest frame
        double Epoc = std::sqrt(pxp*pxp + pyp*pyp + pzp*pzp + m * m);
        // E/c in accelerator frame
        double Eoc =gamma * (Epoc+ beta*pzp);
        double pz = gamma * (pzp + beta * Epoc);



        double  betaz_part=pz/Eoc;
        double  betax_part=pxp/Eoc;
        double  betay_part=pyp/Eoc;
        double  cdt=- particles[part][Bunch::z]*gamma*(1.0-beta*betaz_part)/betaz_part;
        particles[part][Bunch::cdt] = cdt;
        particles[part][Bunch::x] = particles[part][Bunch::x] + cdt*betax_part/(1.0-beta*betaz_part);
        particles[part][Bunch::y] = particles[part][Bunch::y] + cdt*betay_part/(1.0-beta*betaz_part);

        double p = std::sqrt(pxp * pxp + pyp * pyp + pz * pz);
        particles[part][Bunch::dpop] = (p - p_ref) / p_ref;

    }

}

void
Fixed_t_z_alex::from_z_lab_to_t_lab(Bunch &bunch)
{

     double m = bunch.get_mass();
     double p_ref = bunch.get_reference_particle().get_momentum();
     MArray2d_ref particles = bunch.get_local_particles();
     for (int part = 0; part < bunch.get_local_num(); ++part) {
           // total momentum in accelerator frame
         double p = p_ref + particles[part][Bunch::dpop] * p_ref;
          // p_{x,y,z} in accelerator frame
         double px = particles[part][Bunch::xp] * p_ref;
         double py = particles[part][Bunch::yp] * p_ref;
         double  cdt=particles[part][Bunch::cdt];
         double pz2 = p * p - px * px - py * py;
         if (pz2 < 0.0) {
             throw std::runtime_error(
                     "Fixed_t_z_alex::from_z_lab_to_t_lab: particle has negative pz^2");
         }
         double pz = std::sqrt(pz2);
         double Eoc = std::sqrt(p * p + m * m);

         double  betaz_part=pz/Eoc;
         double  betax_part=px/Eoc;
         double  betay_part=py/Eoc;
           // x,y,z, zp in beam rest frame
         particles[part][Bunch::z] = - cdt*betaz_part;
         particles[part][Bunch::x] = particles[part][Bunch::x] - cdt*betax_part;
         particles[part][Bunch::y] = particles[part][Bunch::y] - cdt*betay_part;
         // zp = pz/p_{ref}^{total}
         particles[part][Bunch::zp] = pz/p_ref;

     }
}


void
Fixed_t_z_alex::from_t_lab_to_z_lab(Bunch &bunch)
{
     double m = bunch.get_mass();
     double p_ref = bunch.get_reference_particle().get_momentum();
     MArray2d_ref particles = bunch.get_local_particles();
     for (int part = 0; part < bunch.get_local_num(); ++part) {
         double px  = particles[part][Bunch::xp] * p_ref;
         double py  = particles[part][Bunch::yp] * p_ref;
         double pz  = particles[part][Bunch::zp]  * p_ref;
         double p=std::sqrt(px * px + py * py + pz * pz);
         double Eoc   = std::sqrt(p * p + m * m);
         double  betaz_part=pz/Eoc;
         double  betax_part=px/Eoc;
         double  betay_part=py/Eoc;
         double cdt = -particles[part][Bunch::z]/betaz_part;
         particles[part][Bunch::cdt] = cdt;
         particles[part][Bunch::x] = particles[part][Bunch::x] + cdt*betax_part;
         particles[part][Bunch::y] = particles[part][Bunch::y] +  cdt*betay_part;
         particles[part][Bunch::dpop]=(p-p_ref)/p_ref;
     }

}


void
Fixed_t_z_alex::from_t_lab_to_t_bunch(Bunch &bunch)
{
     double gamma=bunch.get_reference_particle().get_gamma();
     double gamma_inv = 1.0/gamma;
     double beta = bunch.get_reference_particle().get_beta();
     double m = bunch.get_mass();
     double p_ref = bunch.get_reference_particle().get_momentum();
     MArray2d_ref particles = bunch.get_local_particles();
     for (int part = 0; part < bunch.get_local_num(); ++part) {
         double px  = particles[part][Bunch::xp] * p_ref;
         double py  = particles[part][Bunch::yp] * p_ref;
         double pz  = particles[part][Bunch::zp]  * p_ref;
         double p=std::sqrt(px * px + py * py + pz * pz);
         double Eoc   = std::sqrt(p * p + m * m);
         double  betaz_part=pz/Eoc;
         double  betax_part=px/Eoc;
         double  betay_part=py/Eoc;
         double dz=particles[part][Bunch::z];
         particles[part][Bunch::z] =  gamma_inv * particles[part][Bunch::z]/(1.-beta* betaz_part);
         particles[part][Bunch::x]= particles[part][Bunch::x]+dz*beta*betax_part/(1.-beta* betaz_part);
         particles[part][Bunch::y]= particles[part][Bunch::y]+dz*beta*betay_part/(1.-beta* betaz_part);
         particles[part][Bunch::zp] = gamma*(pz - beta * Eoc)/p_ref;
    }
}

void
Fixed_t_z_alex::from_t_bunch_to_t_lab(Bunch &bunch)
{
    double gamma = bunch.get_reference_particle().get_gamma();
    double beta = bunch.get_reference_particle().get_beta();
    double m = bunch.get_mass();
    double p_ref = bunch.get_reference_particle().get_momentum();
    MArray2d_ref particles = bunch.get_local_particles();
    for (int part = 0; part < bunch.get_local_num(); ++part) {
          // p'_{x,y,z} in beam rest frame
        double pxp = particles[part][Bunch::xp]* p_ref;
        double pyp = particles[part][Bunch::yp]* p_ref;
        double pzp = particles[part][Bunch::zp]* p_ref;

        // E'/c in beam rest frame
        double Epoc = std::sqrt(pxp*pxp + pyp*pyp + pzp*pzp + m * m);
        // E/c in accelerator frame
        double Eoc =gamma * (Epoc+ beta*pzp);
        double pz = gamma * (pzp + beta * Epoc);
        double  betaz_part=pz/Eoc;
        double  betax_part=pxp/Eoc;
        double  betay_part=pyp/Eoc;
        double dz = particles[part][Bunch::z]*gamma*(1.0-beta*betaz_part); // in t_lab frame
        particles[part][Bunch::z] = dz;
        particles[part][Bunch::x] = particles[part][Bunch::x]-dz*beta*betax_part/(1.0-beta*betaz_part);
        particles[part][Bunch::y] = particles[part][Bunch::y]-dz*beta*betay_part/(1.0-beta*betaz_part);
        particles[part][Bunch::dpop] = pz / p_ref;
    }
}
BOOST_CLASS_EXPORT_IMPLEMENT(Fixed_t_z_alex)

//*********************Fixed_t_z_alex end*********************************************

//*********************Fixed_t_z_synergia20*****************************************


void
Fixed_t_z_synergia20::from_z_lab_to_t_bunch(Bunch &bunch)
{
    double gamma=bunch.get_reference_particle().get_gamma();
    double gamma_inv = 1.0/gamma;
    double beta = bunch.get_reference_particle().get_beta();
    double m = bunch.get_mass();
    double p_ref = bunch.get_reference_particle().get_momentum();
    MArray2d_ref particles = bunch.get_local_particles();
    for (int part = 0; part < bunch.get_local_num(); ++part) {
         // p_{x,y} in accelerator frame
        double px = particles[part][Bunch::xp] * p_ref;
        double py = particles[part][Bunch::yp] * p_ref;
        // total momentum in accelerator frame
        double p = p_ref + particles[part][Bunch::dpop] * p_ref;
        double pz2 = p * p - px * px - py * py;
        if (pz2 < 0.0) {
            throw std::runtime_error(
                    "Fixed_t_z_zeroth::fixed_z_to_fixed_t: particle has negative pz^2");
        }
        double pz = std::sqrt(pz2);
         // E/c in accelerator frame
        double Eoc = std::sqrt(p * p + m * m);
        double  betaz_part=pz/Eoc;

        particles[part][Bunch::z] = -particles[part][Bunch::cdt]*gamma*betaz_part;
       // zp = pz/p_{ref}^{total}, pz in the beam frame,  p_ref in the accelerator frame
        particles[part][Bunch::zp] = gamma*(pz - beta * Eoc)/p_ref;
    }

}


void Fixed_t_z_synergia20::from_t_bunch_to_z_lab(Bunch &bunch)
{
    double gamma = bunch.get_reference_particle().get_gamma();
    double beta = bunch.get_reference_particle().get_beta();
    double m = bunch.get_mass();
    double p_ref = bunch.get_reference_particle().get_momentum();
    MArray2d_ref particles = bunch.get_local_particles();
    for (int part = 0; part < bunch.get_local_num(); ++part) {
          // p'_{x,y,z} in beam rest frame
        double pxp = particles[part][Bunch::xp]* p_ref;
        double pyp = particles[part][Bunch::yp]* p_ref;
        double pzp = particles[part][Bunch::zp]* p_ref;

        // E'/c in beam rest frame
        double Epoc = std::sqrt(pxp*pxp + pyp*pyp + pzp*pzp + m * m);
        // E/c in accelerator frame
        double Eoc =gamma * (Epoc+ beta*pzp);
        double pz = gamma * (pzp + beta * Epoc);

        double  betaz_part=pz/Eoc;
        particles[part][Bunch::cdt] = - particles[part][Bunch::z]/gamma/betaz_part;

        double p = std::sqrt(pxp * pxp + pyp * pyp + pz * pz);
        particles[part][Bunch::dpop] = (p - p_ref) / p_ref;

    }
}

void
Fixed_t_z_synergia20::from_z_lab_to_t_lab(Bunch &bunch)
{

    double m = bunch.get_mass();
    double p_ref = bunch.get_reference_particle().get_momentum();
    MArray2d_ref particles = bunch.get_local_particles();
    for (int part = 0; part < bunch.get_local_num(); ++part) {
          // total momentum in accelerator frame
        double p = p_ref + particles[part][Bunch::dpop] * p_ref;
         // p_{x,y,z} in accelerator frame
        double px = particles[part][Bunch::xp] * p_ref;
        double py = particles[part][Bunch::yp] * p_ref;
        double pz2 = p * p - px * px - py * py;
        if (pz2 < 0.0) {
            throw std::runtime_error(
                    "Fixed_t_z_zeroth::fixed_z_to_fixed_t: particle has negative pz^2");
        }
        double pz = std::sqrt(pz2);
         // E/c in accelerator frame
        double Eoc = std::sqrt(p * p + m * m);

        double  betaz_part=pz/Eoc;
        particles[part][Bunch::z] = - particles[part][Bunch::cdt]*betaz_part;

        // zp = pz/p_{ref}^{total}
         particles[part][Bunch::zp] =pz/p_ref;

    }
}

void Fixed_t_z_synergia20::from_t_lab_to_z_lab(Bunch &bunch)
{

    double m = bunch.get_mass();
    double p_ref = bunch.get_reference_particle().get_momentum();
    MArray2d_ref particles = bunch.get_local_particles();

   for (int part = 0; part < bunch.get_local_num(); ++part) {
          // p'_{x,y,z} in beam rest frame
        double px = particles[part][Bunch::xp]* p_ref;
        double py = particles[part][Bunch::yp]* p_ref;
        double pz = particles[part][Bunch::zp]* p_ref;

        double Eoc = std::sqrt(px*px + py*py + pz*pz + m * m);
        double p = std::sqrt(px * px + py * py + pz * pz);
        particles[part][Bunch::dpop] = (p - p_ref) / p_ref;

        double  betaz_part=pz/Eoc;
        particles[part][Bunch::cdt] = - particles[part][Bunch::z]/betaz_part;
    }
}



void
Fixed_t_z_synergia20::from_t_lab_to_t_bunch(Bunch &bunch)
{
     double gamma=bunch.get_reference_particle().get_gamma();
     double beta = bunch.get_reference_particle().get_beta();
     double m = bunch.get_mass();
     double p_ref = bunch.get_reference_particle().get_momentum();
     MArray2d_ref particles = bunch.get_local_particles();
     for (int part = 0; part < bunch.get_local_num(); ++part) {
         double px  = particles[part][Bunch::xp] * p_ref;
         double py  = particles[part][Bunch::yp] * p_ref;
         double pz  = particles[part][Bunch::zp]  * p_ref;
         double p=std::sqrt(px * px + py * py + pz * pz);
         double Eoc   = std::sqrt(p * p + m * m);

         particles[part][Bunch::z] =  gamma * particles[part][Bunch::z];

         particles[part][Bunch::zp] = gamma*(pz - beta * Eoc)/p_ref;
    }
}
void
Fixed_t_z_synergia20::from_t_bunch_to_t_lab(Bunch &bunch)
{
    double gamma = bunch.get_reference_particle().get_gamma();
    double beta = bunch.get_reference_particle().get_beta();
    double m = bunch.get_mass();
    double p_ref = bunch.get_reference_particle().get_momentum();
    MArray2d_ref particles = bunch.get_local_particles();
    for (int part = 0; part < bunch.get_local_num(); ++part) {
          // p'_{x,y,z} in beam rest frame
        double pxp = particles[part][Bunch::xp]* p_ref;
        double pyp = particles[part][Bunch::yp]* p_ref;
        double pzp = particles[part][Bunch::zp]* p_ref;
        // E'/c in beam rest frame
        double Epoc = std::sqrt(pxp*pxp + pyp*pyp + pzp*pzp + m * m);
        double pz = gamma * (pzp + beta * Epoc);
        particles[part][Bunch::z] = particles[part][Bunch::z]/gamma;
        particles[part][Bunch::dpop] = pz / p_ref;

    }
}
BOOST_CLASS_EXPORT_IMPLEMENT(Fixed_t_z_synergia20)


//*********************Fixed_t_z_synergia20 end*****************************************
