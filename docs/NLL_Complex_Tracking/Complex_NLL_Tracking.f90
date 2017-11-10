
      subroutine NonlinearLensPropagatorCmplx(knll,cnll,coord)      
    !*****************************************************************
    ! The following subroutine computes the nonlinear momentum kick
    ! across a thin lens associated with a single short segment of the
    ! nonlinear magnetic insert described in V. Danilov and S. Nagaitsev,
    ! PRSTAB 13, 084002 (2010), Sect. V.A.  The arguments are as follows:
    !         knll - integrated strength of the lens (m)
    !         cnll - distance of singularities from the origin (m)
    !         coord = (x [m], px/p0, y [m], py/p0)
    ! This implementation is based on expressions in "Nonlinear Lens 
    ! Tracking in the IOTA Complex Potential," C. Mitchell, Feb. 9, 2017.
    ! Variable definitions are chosen to be consistent with TRACK_EXT.
    ! C. Mitchell 2/9/2017
    !*****************************************************************
        implicit none
!        include 'mpif.h'
        double precision, intent(in) :: knll,cnll
        double precision, dimension(4), intent(inout) :: coord
        double complex:: dF,Fderivative
        double precision :: x,y,kick,dPx,dPy
        x = coord(1)/cnll                          !Dimensionless horizontal coord
        y = coord(3)/cnll                          !Dimensionless vertical coord
        kick = -knll/cnll                          !Dimensionless kick strength
      !Avoid the branch cuts:
        if((y==0.d0).and.(dabs(x).ge.1.d0)) then
          write(*,*) "Error:  NonlinearLensPropagatorCmplx propagates across branch cuts!"
          return
        else
          dF = Fderivative(x,y)
          dPx = kick*real(dF)
          dPy = -kick*aimag(dF)
        endif
      !Momentum update
        coord(2)=coord(2)+dPx
        coord(4)=coord(4)+dPy
      end subroutine NonlinearLensPropagatorCmplx 

  
     subroutine InvariantPotentials(x,y,Hinv,Iinv)
   !*****************************************************************
   ! This subroutine computes the dimensionless potentials that 
   ! define the contribution of the NLI to the two invariants 
   ! (H,I) for the IOTA ring.
   ! The arguments are as follows:
   !       (x,y) - normalized dimensionless coordinates
   !       Hinv - the vector potential describing the spatial
   !              dependence of the first invariant H.
   !       Iinv - function describing the spatial dependence
   !              of the second invariant I.
   ! This implementation is based on expressions in "Nonlinear Lens 
   ! Tracking in the IOTA Complex Potential," C. Mitchell, Feb. 9, 2017.
   ! C. Mitchell, 2/9/2017.
   !*****************************************************************
     implicit none
     double precision, intent(in):: x,y
     double precision, intent(out):: Hinv,Iinv
     double complex:: zeta,zetaconj,Hpotential,Ipotential
     double complex:: croot,carcsin
     zeta = dcmplx(x,y)
     zetaconj = conjg(zeta)
     Hpotential = zeta/croot(zeta)
     Ipotential = (zeta+zetaconj)/croot(zeta)   
     Hpotential = Hpotential*carcsin(zeta)
     Ipotential = Ipotential*carcsin(zeta)
     Hinv = real(Hpotential)
     Iinv = real(Ipotential)
     end subroutine 


     function Fpotential(x,y)
   !****************************************
   ! Computes the dimensionless complex
   ! potential Az+i*Psi of the IOTA
   ! nonlinear insert.
   !****************************************
     implicit none
     double precision, intent(in):: x,y
     double complex:: Fpotential,zeta
     double complex:: croot,carcsin
     zeta = dcmplx(x,y)
     Fpotential = zeta/croot(zeta)
     Fpotential = Fpotential*carcsin(zeta)
     end function

     function Fderivative(x,y)
   !****************************************
   ! Computes the derivative of the
   ! dimensionless complex potential for
   ! the IOTA nonlinear insert.
   !****************************************
     implicit none
     double precision, intent(in):: x,y
     double complex:: Fderivative,zeta
     double complex:: denom,croot,carcsin
     zeta = dcmplx(x,y)
     denom = croot(zeta)
     Fderivative = zeta/denom**2
     Fderivative = Fderivative+carcsin(zeta)/denom**3
     end function

     function carcsin(z)
   !******************************************
   ! Computes the complex function arcsin(z)
   ! using the principal branch.
   !******************************************
     implicit none
     double complex, intent(in):: z
     double complex:: carcsin,im1,croot
     im1 = dcmplx(0.0d0,1.0d0)
     carcsin = im1*z+croot(z)
     carcsin = -im1*cdlog(carcsin)
     end function

      function croot(z)
   !*******************************************
   ! Computes the complex function sqrt(1-z^2)
   ! using the principal branch.
   !*******************************************
      implicit none
      double complex, intent(in):: z
      double complex:: croot,re1
      re1 = dcmplx(1.0d0,0.0d0)
      croot = re1-z**2
      croot = cdsqrt(croot)
      end function croot
