! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Routine to calculate spectral snow albedos for JULES, based on
! the Marshall (1989) parametrization of the Wiscombe and Warren (1980)
! model. Influence of contaminants in snow has not been included - see
! UM vn4.5 deck FTSA1A.

! *********************************************************************
MODULE albsnow_mod
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ALBSNOW_MOD'
CONTAINS
SUBROUTINE albsnow (land_pts,nsurft,surft_index,surft_pts,                     &
                    cosz_gb,rgrain_surft,snowdepth_surft,alb_snow,             &
                    nsnow_surft,soot_mass,sootm_tot,sice_surft,                &
                    sliq_surft,snow_surft,snow_grnd_surft)

USE jules_surface_types_mod,  ONLY: ntype
USE jules_snow_mod,           ONLY: r0, amax, cansnowtile, nsmax
USE jules_surface_mod,        ONLY: l_aggregate

USE parkind1,                 ONLY: jprb, jpim
USE yomhook,                  ONLY: lhook, dr_hook

USE um_types, ONLY: real_jlslsm

USE ereport_mod, ONLY: ereport

IMPLICIT NONE

! Subroutine arguments
!   Scalar arguments with intent(in):
INTEGER, INTENT(IN) ::                                                         &
 land_pts                                                                      &
                             ! Number of land points.
,nsurft                      ! Number of surface tiles.

!   Array arguments with intent(in):
INTEGER, INTENT(IN) ::                                                         &
 surft_pts(ntype)                                                              &
                             ! Number tile points.
,surft_index(land_pts,ntype)                                                   &
                             ! Index of tile points.
,nsnow_surft(land_pts,nsurft)! Number of snow layers on ground on tiles

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
 cosz_gb(land_pts)                                                             &
                             ! Cosine of the zenith angle.
,rgrain_surft(land_pts,nsurft)                                                 &
                             ! Snow grain size (microns).
,snowdepth_surft(land_pts,nsurft)                                              &
                             ! Depth of snow on the ground (m).
,soot_mass(land_pts,nsurft,nsmax)                                              &
                             ! Snow layer mass of soot on surface tiles (kg/m2)
,sootm_tot(land_pts,nsurft)                                                    &
                             ! Total soot mass on surface tiles (kg/m2)
,sice_surft(land_pts,nsurft,nsmax)                                             &
                             ! Snow layer ice mass on tiles (kg/m2)
,sliq_surft(land_pts,nsurft,nsmax)                                             &
                             ! Snow layer liquid mass on tiles (kg/m2)
,snow_surft(land_pts,nsurft)                                                   &
                             ! Lying snow on tiles
                             ! (in the canopy or on the ground) (kg/m2)
,snow_grnd_surft(land_pts,nsurft)
                             ! Snow on the ground beneath the canopy (kg/m2)

!   Array arguments with intent(out):
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
 alb_snow(land_pts,ntype,4)! Snow albedo.
!                              (*,*,1) - Direct beam visible
!                              (*,*,2) - Diffuse visible
!                              (*,*,3) - Direct beam near-IR
!                              (*,*,4) - Diffuse near-IR

! Local scalars:
REAL(KIND=real_jlslsm) ::                                                      &
 reff                                                                          &
                             ! Zenith effective grain size.
,sigma                                                                         &
                             ! Scaled soot content
,snow_mass_top                                                                 &
                             ! Mass of snow in the top snow layer (kg/m2)
,snow_mass_ground                                                              &
                             ! Total mass of snow beneath the canopy on a
                             ! given surface tile (kg/m2)
,soot_conc
                             ! Soot concentration in snow (kg/kg)

! Precalculated values for optimisation
REAL(KIND=real_jlslsm) ::                                                      &
  sqrt_r0, sigma_power

INTEGER ::                                                                     &
 band,j,l,n                ! Loop counters.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ALBSNOW'

INTEGER                     :: errcode  ! Error code
CHARACTER(LEN=880)          :: ERRMSG   ! Error message

!End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!Pre-calculate for optimisation purposes
sqrt_r0 = SQRT(r0)

DO n = 1,nsurft
  DO l = 1,land_pts
    IF (snowdepth_surft(l,n) > 0.0) THEN
      reff = rgrain_surft(l,n) * ( 1.0 + 0.77 * (cosz_gb(l) - 0.65) )**2

      alb_snow(l,n,1) = amax(1) - 0.002 * (SQRT(reff)              - sqrt_r0)
      alb_snow(l,n,2) = amax(1) - 0.002 * (SQRT(rgrain_surft(l,n)) - sqrt_r0)
      alb_snow(l,n,3) = amax(2) - 0.09  * LOG(reff / r0)
      alb_snow(l,n,4) = amax(2) - 0.09  * LOG(rgrain_surft(l,n) / r0)
    END IF
  END DO
END DO

! Adjust visible snow albedos for soot content
DO n = 1,nsurft
  DO j = 1,surft_pts(n)
    l = surft_index(j,n)
    IF (sootm_tot(l,n) < 0.0) THEN
      errcode = 1
      ERRMSG = 'negative soot mass in snow albedo calculations'
      CALL ereport(routinename, errcode, ERRMSG)
    END IF
    IF (cansnowtile(n)) THEN
      snow_mass_ground = snow_grnd_surft(l,n)
    ELSE
      snow_mass_ground = snow_surft(l,n)
    END IF
    IF (snow_mass_ground /= snow_mass_ground) THEN
      errcode = 1
      ERRMSG = 'NaN values in snow mass fields'
      CALL ereport(routinename, errcode, ERRMSG)
    END IF
    IF ( snow_mass_ground > 0.0 ) THEN
      IF ( nsnow_surft(l,n) > 0 ) THEN
        snow_mass_top = sice_surft(l,n,1) + sliq_surft(l,n,1)
        soot_conc = soot_mass(l,n,1) / snow_mass_top
      ELSE
        soot_conc = sootm_tot(l,n) / snow_mass_ground
      END IF

      sigma = soot_conc * rgrain_surft(l,n) / 0.0017

      IF ( sigma  >   1.0 ) THEN
        sigma_power     = sigma**0.46
        alb_snow(l,n,1) = 0.07 +                                               &
                          0.5 * (alb_snow(l,n,1) - 0.07) / (sigma_power)
        alb_snow(l,n,2) = 0.07 +                                               &
                          0.5 * (alb_snow(l,n,2) - 0.07) / (sigma_power)
      ELSE
        sigma_power     = sigma**0.6
        alb_snow(l,n,1) = alb_snow(l,n,1) -                                    &
                          0.5 * (alb_snow(l,n,1) - 0.07) * (sigma_power)
        alb_snow(l,n,2) = alb_snow(l,n,2) -                                    &
                          0.5 * (alb_snow(l,n,2) - 0.07) * (sigma_power)
      END IF
    END IF
  END DO
END DO

! Adjust near-IR snow albedos for soot content
DO n = 1,nsurft
  DO j = 1,surft_pts(n)
    l = surft_index(j,n)
    ! Calculate soot concentration based on snow on the ground, not in the
    ! canopy. As with snow density, soot mass is not stored in the canopy so
    ! top snow layer values are used.
    IF (cansnowtile(n)) THEN
      snow_mass_ground = snow_grnd_surft(l,n)
    ELSE
      snow_mass_ground = snow_surft(l,n)
    END IF
    IF ( snow_mass_ground  >   0.0 ) THEN
      IF ( nsnow_surft(l,n) > 0 ) THEN
        snow_mass_top = sice_surft(l,n,1) + sliq_surft(l,n,1)
        soot_conc = soot_mass(l,n,1) / snow_mass_top
      ELSE
        soot_conc = sootm_tot(l,n) / snow_mass_ground
      END IF
      sigma = soot_conc * rgrain_surft(l,n) / 0.004

      IF ( sigma > 1.0 ) THEN
        sigma_power     = sigma**0.6
        alb_snow(l,n,3) = 0.06 +                                               &
                          0.5 * (alb_snow(l,n,3) - 0.06) / (sigma_power)
        alb_snow(l,n,4) = 0.06 +                                               &
                          0.5 * (alb_snow(l,n,4) - 0.06) / (sigma_power)
      ELSE
        sigma_power     = sigma**0.7
        alb_snow(l,n,3) = alb_snow(l,n,3) -                                    &
                          0.5 * (alb_snow(l,n,3) - 0.06) * (sigma_power)
        alb_snow(l,n,4) = alb_snow(l,n,4) -                                    &
                          0.5 * (alb_snow(l,n,4) - 0.06) * (sigma_power)
      END IF
    END IF
  END DO
END DO

IF ( l_aggregate ) THEN
  DO band = 1,4
    DO n = 2,ntype
      DO j = 1,surft_pts(n)
        l = surft_index(j,n)
        alb_snow(l,n,band) = alb_snow(l,1,band)
      END DO
    END DO
  END DO
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE albsnow
END MODULE albsnow_mod
