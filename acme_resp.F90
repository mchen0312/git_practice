!---------------------------------------------------------------------------
! acme_resp.F90
! a simplification of the CLM/ACME leaf maintenance respiration code
! particularly in regard to the leaf to canopy scaling factor: vcmaxcint
!---------------------------------------------------------------------------

subroutine amresp(lmr_z,tl,lnc,elai)
implicit none

!---------------Description-------------------------------------------------
! takes leaf nitrogen content (lnc), leaf temperature (tl) and leaf area
! index (lai) as inputs to leaf maintentance respiration (lmr_z) 
! currently only for a single big leaf canopy layer
! Note: expects input T[K] and lnc [(gN leaf)/m^2]
!---------------------------------------------------------------------------

real*8 :: lmr_z,tl,lnc,elai
real*8 :: lmr25,lmr25top,nscaler,extkn,lmrc,lmrse,lmrhd,lmrha,rgas,tfrz 

!f2py intent(in) lnc,tl,elai
!f2py intent(out) lmr_z
!f2py depend(lmr_z) lnc,tl,elai,extkn,lmrse,lmrhd,lmrha

tfrz = 273.15_8
rgas = 6.02214e26_8 * 1.38065e-23_8

lmrha   = 46390._8
lmrhd   = 150650._8
lmrse   = 490._8
lmrc    = 1._8 + exp( (-lmrhd+lmrse*(tfrz+25._8)) / (rgas*1.e-3_8*(tfrz+25._8)) )

! substitute a simplified eqn. for vcmaxcint, no sun/sha, no burial w/snow
! nscaler = vcmaxcint
extkn = 0.30_8
nscaler = (1._8 - exp(-extkn*elai)) / extkn

lmr25top = 2.525e-6_8 * (1.5_8 ** ((25._8 - 20._8)/10._8))
lmr25top = lmr25top * lnc / 12.e-06_8

lmr25 = lmr25top * nscaler

lmr_z = lmr25 * &
        exp( lmrha / (rgas*1.e-3_8*(tfrz+25._8)) * (1._8 - (tfrz+25._8)/tl) ) * &
        lmrc / ( 1._8 + exp( (-lmrhd+lmrse*tl) / (rgas*1.e-3_8*tl) ) )

end subroutine amresp
