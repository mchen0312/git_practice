!------------------------------------------------------------------
! hunt_resp.f90
! an implementation of Huntingford, et al's maintenance respiration
! formulation as will be used in the ACME land module
!------------------------------------------------------------------

subroutine mresp(Rdark,r0,r1,r2,nl,Tg,Ta)
implicit none

!---------- Description ----------------------------------------
! takes four pft dependent parameters from an expoential fit to
! data (citation) r0,r1,r2, and nl
! as well as two environmental variables: Tg (10 day avg. T, in C) 
! and Ta, instantaneous leaf temperature (in C)
!---------------------------------------------------------------

real*8 :: Rdark,r0,r1,r2,nl,Tg,Ta

!f2py intent(in) r0,r1,r2,nl,Tg,Ta
!f2py intent(out) Rdark
!f2py depend(Rdark) r0,r1,r2,nl,Tg,Ta

Rdark = (r0 + r1 * nl - r2 * Tg) * &
         exp( 0.1012_8 * (Ta - 25.0_8) - 0.0005_8 * (Ta**2 - 25.0_8**2) )

end subroutine mresp
