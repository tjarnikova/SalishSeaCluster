MODULE p4zmort
   !!======================================================================
   !!                         ***  MODULE p4zmort  ***
   !! TOP :   SMELT Compute the mortality terms for phytoplankton
   !!======================================================================
   !! History :   2015 (E. Olson) adapted from 1-d SOG and PISCES
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   'key_pisces'                               SMELT / PISCES bio-model
   !!----------------------------------------------------------------------
   !!   p4z_mort       :   Compute the mortality terms for phytoplankton
   !!   p4z_mort_init  :   Initialize the mortality params for phytoplankton
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  SMELT Source Minus Sink variables
   !USE p4zsink         !  vertical flux of particulate matter due to sinking
   !USE p4zprod         !  Primary productivity 
   USE prtctl_trc      !  print control for debugging
   USE p4zprod, ONLY: zz_rate_Si_ratio_diat, zz_rate_Si_ratio_myri, zz_rate_Si_ratio_nano
   USE iom

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_mort    
   PUBLIC   p4z_mort_init    

   !! * Shared module variables
   REAL(wp) ::  zz_rate_pico_Rm        ! picophyt natural mort
   REAL(wp) ::  zz_frac_waste_FNM_NH   ! waste frac from picophyt natural mort to NH
   REAL(wp) ::  zz_frac_waste_FNM_DON  ! waste frac from picophyt natural mort to DON
   REAL(wp) ::  zz_frac_waste_FNM_PON  ! waste frac from picophyt natural mort to PON
   REAL(wp) ::  zz_frac_waste_FNM_Ref  ! waste frac from picophyt natural mort to refractory
   REAL(wp) ::  zz_frac_waste_FNM_Bsi  ! waste frac from picophyt natural mort to bSi
   REAL(wp) ::  zz_rate_myri_Rm        ! picophyt natural mort
   REAL(wp) ::  zz_frac_waste_NNM_NH   ! waste frac from picophyt natural mort to NH
   REAL(wp) ::  zz_frac_waste_NNM_DON  ! waste frac from picophyt natural mort to DON
   REAL(wp) ::  zz_frac_waste_NNM_PON  ! waste frac from picophyt natural mort to PON
   REAL(wp) ::  zz_frac_waste_NNM_Ref  ! waste frac from picophyt natural mort to refractory
   REAL(wp) ::  zz_frac_waste_NNM_Bsi  ! waste frac from picophyt natural mort to bSi
   REAL(wp) ::  zz_rate_micro_Rm       ! micro-phyto natural mort
   REAL(wp) ::  zz_frac_waste_DNM_NH   ! waste frac from micro-phyto natural mort to NH
   REAL(wp) ::  zz_frac_waste_DNM_DON  ! waste frac from micro-phyto natural mort to DON
   REAL(wp) ::  zz_frac_waste_DNM_PON  ! waste frac from micro-phyto natural mort to PON
   REAL(wp) ::  zz_frac_waste_DNM_Ref  ! waste frac from micro-phyto natural mort to refractory
   REAL(wp) ::  zz_frac_waste_DNM_Bsi  ! waste frac from micro-phyto natural mort to bSi
      

   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: p4zmort.F90 5385 2015-06-09 13:50:42Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE p4z_mort( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_mort  ***
      !!
      !! ** Purpose :   Calls the different subroutine to initialize and compute
      !!                the different phytoplankton mortality terms
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt ! ocean time step
      !!---------------------------------------------------------------------

      CALL p4z_nano            ! nanophytoplankton
      
      CALL p4z_myri            ! M. rubrum

      CALL p4z_diat            ! diatoms

   END SUBROUTINE p4z_mort


   SUBROUTINE p4z_nano
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_nano  ***
      !!
      !! ** Purpose :   Compute the mortality terms for nanophytoplankton
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER  :: ji, jj, jk
      CHARACTER (len=25) :: charout
      REAL(wp), POINTER, DIMENSION(:,:,:) ::  zz_NatMort_pico(:,:,:)
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_nano')
      !
      CALL wrk_alloc(jpi, jpj, jpk, zz_NatMort_pico)
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
            
               zz_NatMort_pico(ji,jj,jk) = zz_rate_pico_Rm * tgfunc(ji,jj,jk) * trb(ji,jj,jk,jpphy) * rfact2 
               tra(ji,jj,jk,jpphy) = tra(ji,jj,jk,jpphy) - zz_NatMort_pico(ji,jj,jk)
               tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) + zz_frac_waste_FNM_NH * zz_NatMort_pico(ji,jj,jk)
               tra(ji,jj,jk,jpdon) = tra(ji,jj,jk,jpdon) + zz_frac_waste_FNM_DON * zz_NatMort_pico(ji,jj,jk)
               tra(ji,jj,jk,jppon) = tra(ji,jj,jk,jppon) + zz_frac_waste_FNM_PON * zz_NatMort_pico(ji,jj,jk)
               tra(ji,jj,jk,jpdsi) = tra(ji,jj,jk,jpdsi) + zz_frac_waste_FNM_Bsi * zz_NatMort_pico(ji,jj,jk) * zz_rate_Si_ratio_nano

            END DO
         END DO
      END DO

      CALL iom_put( "MORTPHY", zz_NatMort_pico/rfact2)
      CALL wrk_dealloc(jpi, jpj, jpk, zz_NatMort_pico)
       IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('nano')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
       ENDIF

      IF( nn_timing == 1 )  CALL timing_stop('p4z_nano')
      !
   END SUBROUTINE p4z_nano

   SUBROUTINE p4z_myri
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_myri  ***
      !!
      !! ** Purpose :   Compute the mort terms for M. rubrum
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER  :: ji, jj, jk
      CHARACTER (len=25) :: charout
      ! SOG vars:

      REAL(wp), POINTER, DIMENSION(:,:,:) ::  zz_NatMort_myri(:,:,:)
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_myri')
      !
      CALL wrk_alloc(jpi, jpj, jpk, zz_NatMort_myri)
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
            
               zz_NatMort_myri(ji,jj,jk) = zz_rate_myri_Rm * tgfunc(ji,jj,jk) * trb(ji,jj,jk,jpmes) * rfact2
               tra(ji,jj,jk,jpmes) = tra(ji,jj,jk,jpmes) - zz_NatMort_myri(ji,jj,jk)
               tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) + zz_frac_waste_NNM_NH * zz_NatMort_myri(ji,jj,jk)
               tra(ji,jj,jk,jpdon) = tra(ji,jj,jk,jpdon) + zz_frac_waste_NNM_DON * zz_NatMort_myri(ji,jj,jk)
               tra(ji,jj,jk,jppon) = tra(ji,jj,jk,jppon) + zz_frac_waste_NNM_PON * zz_NatMort_myri(ji,jj,jk)
               tra(ji,jj,jk,jpdsi) = tra(ji,jj,jk,jpdsi) + zz_frac_waste_NNM_Bsi * zz_NatMort_myri(ji,jj,jk) * zz_rate_Si_ratio_myri

            END DO
         END DO
      END DO
      !
      CALL iom_put( "MORTMRUB", zz_NatMort_myri/rfact2)
      CALL wrk_dealloc(jpi, jpj, jpk, zz_NatMort_myri)
       IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('myri')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
       ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop('p4z_myri')
      !
   END SUBROUTINE p4z_myri

   SUBROUTINE p4z_diat
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_diat  ***
      !!
      !! ** Purpose :   Compute the mort terms for diatoms
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER  ::  ji, jj, jk
      CHARACTER (len=25) :: charout
      ! SOG vars:
      REAL(wp), POINTER, DIMENSION(:,:,:) ::  zz_NatMort_micro(:,:,:)
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_diat')
      !
      CALL wrk_alloc(jpi, jpj, jpk, zz_NatMort_micro)

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi

               zz_NatMort_micro(ji,jj,jk) = zz_rate_micro_Rm * tgfunc(ji,jj,jk)*trb(ji,jj,jk,jpdia) * rfact2
               tra(ji,jj,jk,jpdia) = tra(ji,jj,jk,jpdia) - zz_NatMort_micro(ji,jj,jk)
               tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) + zz_frac_waste_DNM_NH * zz_NatMort_micro(ji,jj,jk)
               tra(ji,jj,jk,jpdon) = tra(ji,jj,jk,jpdon) + zz_frac_waste_DNM_DON * zz_NatMort_micro(ji,jj,jk)
               tra(ji,jj,jk,jppon) = tra(ji,jj,jk,jppon) + zz_frac_waste_DNM_PON * zz_NatMort_micro(ji,jj,jk)
               tra(ji,jj,jk,jpdsi) = tra(ji,jj,jk,jpdsi) + zz_frac_waste_DNM_Bsi * zz_NatMort_micro(ji,jj,jk) * zz_rate_Si_ratio_diat
               
            END DO
         END DO
      END DO

      CALL iom_put( "MORTDIAT", zz_NatMort_micro/rfact2)
      !
      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('diat')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
      CALL wrk_dealloc(jpi, jpj, jpk, zz_NatMort_micro)
      IF( nn_timing == 1 )  CALL timing_stop('p4z_diat')
      !
   END SUBROUTINE p4z_diat


   SUBROUTINE p4z_mort_init

      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_mort_init  ***
      !!
      !! ** Purpose :   Initialization of phytoplankton parameters
      !!
      !! ** Method  :   Read the nampismort namelist and check the parameters
      !!      called at the first timestep
      !!
      !! ** input   :   Namelist nampismort
      !!
      !!----------------------------------------------------------------------

      NAMELIST/nampismort/ zz_rate_pico_Rm, zz_frac_waste_FNM_NH, zz_frac_waste_FNM_DON, &
      &  zz_frac_waste_FNM_PON, zz_frac_waste_FNM_Ref, zz_frac_waste_FNM_Bsi, &
      &                    zz_rate_myri_Rm, zz_frac_waste_NNM_NH, zz_frac_waste_NNM_DON, &
      &  zz_frac_waste_NNM_PON, zz_frac_waste_NNM_Ref, zz_frac_waste_NNM_Bsi,  &
      &                    zz_rate_micro_Rm, zz_frac_waste_DNM_NH, zz_frac_waste_DNM_DON, &
      &  zz_frac_waste_DNM_PON, zz_frac_waste_DNM_Ref, zz_frac_waste_DNM_Bsi
      INTEGER :: ios                 ! Local integer output status for namelist read

      REWIND( numnatp_ref )              ! Namelist nampismort in reference namelist : phytoplankton mortality
      READ  ( numnatp_ref, nampismort, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampismort in reference namelist', lwp )

      REWIND( numnatp_cfg )              ! Namelist nampismort in configuration namelist : phytoplankton mortality
      READ  ( numnatp_cfg, nampismort, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampismort in configuration namelist', lwp )
      IF(lwm) WRITE ( numonp, nampismort )

      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for phytoplankton mort, nampismort'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '  picophyt natural mort                               zz_rate_pico_Rm  =', zz_rate_pico_Rm
         WRITE(numout,*) '  waste frac from picophyt natural mort to NH    zz_frac_waste_FNM_NH  =', zz_frac_waste_FNM_NH
         WRITE(numout,*) '  waste frac from picophyt natural mort to DON   zz_frac_waste_FNM_DON =', zz_frac_waste_FNM_DON
         WRITE(numout,*) '  waste frac from picophyt natural mort to PON   zz_frac_waste_FNM_PON =', zz_frac_waste_FNM_PON
         WRITE(numout,*) '  waste frac picophyt natural mort to refractory zz_frac_waste_FNM_Ref =', zz_frac_waste_FNM_Ref
         WRITE(numout,*) '  waste frac from picophyt natural mort to bSi   zz_frac_waste_FNM_Bsi =', zz_frac_waste_FNM_Bsi
         WRITE(numout,*) '  picophyt si/n ratio                            zz_rate_Si_ratio_nano =', zz_rate_Si_ratio_nano
         WRITE(numout,*) '  picophyt natural mort                               zz_rate_myri_Rm  =', zz_rate_myri_Rm
         WRITE(numout,*) '  waste frac from picophyt natural mort to NH    zz_frac_waste_NNM_NH  =', zz_frac_waste_NNM_NH
         WRITE(numout,*) '  waste frac from picophyt natural mort to DON   zz_frac_waste_NNM_DON =', zz_frac_waste_NNM_DON
         WRITE(numout,*) '  waste frac from picophyt natural mort to PON   zz_frac_waste_NNM_PON =', zz_frac_waste_NNM_PON
         WRITE(numout,*) '  waste frac picophyt natural mort to refractory zz_frac_waste_NNM_Ref =', zz_frac_waste_NNM_Ref
         WRITE(numout,*) '  waste frac from picophyt natural mort to bSi   zz_frac_waste_NNM_Bsi =', zz_frac_waste_NNM_Bsi
         WRITE(numout,*) '  picophyt si/n ratio                            zz_rate_Si_ratio_myri =', zz_rate_Si_ratio_myri
         WRITE(numout,*) '  micro-phyto natural mort                           zz_rate_micro_Rm  =', zz_rate_micro_Rm
         WRITE(numout,*) '  waste frac from microphyt natural mort to NH   zz_frac_waste_DNM_NH  =', zz_frac_waste_DNM_NH
         WRITE(numout,*) '  waste frac from microphyt natural mort to DON  zz_frac_waste_DNM_DON =', zz_frac_waste_DNM_DON
         WRITE(numout,*) '  waste frac from microphyt natural mort to PON  zz_frac_waste_DNM_PON =', zz_frac_waste_DNM_PON
         WRITE(numout,*) '  waste frac microphyt natural mort to refractory zz_frac_waste_DNM_Ref=', zz_frac_waste_DNM_Ref
         WRITE(numout,*) '  waste frac from microphyt natural mort to bSi  zz_frac_waste_DNM_Bsi =', zz_frac_waste_DNM_Bsi
         WRITE(numout,*) '  micro-phyto si/n ratio                        zz_rate_Si_ratio_diat =', zz_rate_Si_ratio_diat
      ENDIF

   END SUBROUTINE p4z_mort_init

#else
   !!======================================================================
   !!  Dummy module :                                   No bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_mort                    ! Empty routine
   END SUBROUTINE p4z_mort
#endif 
   !!======================================================================
END MODULE p4zmort
