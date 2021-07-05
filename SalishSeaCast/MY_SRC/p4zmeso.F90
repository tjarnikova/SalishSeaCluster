MODULE p4zmeso
   !!======================================================================
   !!                         ***  MODULE p4zmeso  ***
   !! TOP :   SMELT Compute the sources/sinks for M. rubrum due to 
   !!          heterotrophy
   !!======================================================================
   !! History :   2015 (E. Olson) adapted from 1-d SOG and PISCES
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   'key_pisces'                               SMELT / PISCES bio-model
   !!----------------------------------------------------------------------
   !!   p4z_meso       : Compute the sources/sinks for M. rubrum grazing
   !!   p4z_meso_init  : Initialization of M. rubrum grazing the parameters 
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  Source Minus Sink variables
   USE p4zint          !  interpolation and computation of various fields
   USE prtctl_trc      !  print control for debugging
   USE iom             !  I/O manager
   USE p4zprod, ONLY: zz_rate_Si_ratio_nano

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_meso              ! called in p4zbio.F90
   PUBLIC   p4z_meso_init         ! called in trcsms_pisces.F90

   !! * Shared module variables
   REAL(wp) ::   zz_rate_mesorub_R              !:
   REAL(wp) ::   zz_rate_mesorub_PicoPredSlope  !:
   REAL(wp) ::   zz_rate_mesorub_PicoHalfSat    !:
   REAL(wp) ::   zz_rate_mesorub_eff            !:
   REAL(wp) ::   zz_frac_waste_FEN_NH           !:
   REAL(wp) ::   zz_frac_waste_FEN_DON          !waste fraction from mesorub grazing picophyto to DON
   REAL(wp) ::   zz_frac_waste_FEN_PON          !waste fraction from mesorub grazing picophyto to PON
   REAL(wp) ::   zz_frac_waste_FEN_BSi          !waste fraction from mesorub grazing picophyto to Bsi
      
   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: p4zmeso.F90 3160 2011-11-20 14:27:18Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE p4z_meso( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_meso  ***
      !!
      !! ** Purpose :   Compute the sources/sinks for M. rubrum grazing
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt, knt ! ocean time step
      INTEGER  :: ji, jj, jk
      CHARACTER (len=25) :: charout
      REAL(wp) :: zz_Mesorub_mort_pico, zz_Ppico, zz_Pnano, zz_was_NH, zz_Mesorub_eat, &
          zz_was_DON, zz_was_PON, zz_was_BSi     
      REAL(wp), POINTER, DIMENSION(:,:,:) :: het_growth_MRub
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_meso')
      
      CALL wrk_alloc( jpi, jpj, jpk, het_growth_MRub )

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
                zz_Ppico = trb(ji,jj,jk,jpphy)
                zz_Pnano = trb(ji,jj,jk,jpmes)
                zz_Mesorub_mort_pico = zz_rate_mesorub_R * (zz_Ppico - zz_rate_mesorub_PicoPredSlope) &
                    / (zz_rate_mesorub_PicoHalfSat + zz_Ppico - zz_rate_mesorub_PicoPredSlope &
                    + epsilon(zz_rate_mesorub_PicoHalfSat)) &
                    * zz_Pnano * tgfunc(ji,jj,jk)
                zz_Mesorub_mort_pico = max(zz_Mesorub_mort_pico, 0.d0) * rfact2

                zz_was_NH = zz_frac_waste_FEN_NH * zz_Mesorub_mort_pico * &
                     (1-zz_rate_mesorub_eff)
                zz_was_DON = zz_frac_waste_FEN_DON * zz_Mesorub_mort_pico * &
                     (1-zz_rate_mesorub_eff)
                zz_was_PON = zz_frac_waste_FEN_PON * zz_Mesorub_mort_pico * &
                     (1-zz_rate_mesorub_eff)
                zz_was_BSi = zz_frac_waste_FEN_BSi * zz_Mesorub_mort_pico * &
                     zz_rate_Si_ratio_nano ! all Si in grazed nanos directly to BSi
                zz_Mesorub_eat = zz_Mesorub_mort_pico * zz_rate_mesorub_eff
                
               !   Update the arrays TRA which contain the biological sources and sinks
               tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) + zz_was_NH
               tra(ji,jj,jk,jpdon) = tra(ji,jj,jk,jpdon) + zz_was_DON
               tra(ji,jj,jk,jppon) = tra(ji,jj,jk,jppon) + zz_was_PON
               tra(ji,jj,jk,jpmes) = tra(ji,jj,jk,jpmes) + zz_Mesorub_eat
               tra(ji,jj,jk,jpphy) = tra(ji,jj,jk,jpphy) - zz_Mesorub_mort_pico
               tra(ji,jj,jk,jpdsi) = tra(ji,jj,jk,jpdsi) + zz_was_BSi
               het_growth_MRub(ji,jj,jk) = zz_Mesorub_eat/rfact2
            END DO
         END DO
      END DO

      CALL iom_put( "HetMRub", het_growth_MRub)
      CALL wrk_dealloc( jpi, jpj, jpk, het_growth_MRub )
      
      IF( nn_timing == 1 )  CALL timing_stop('p4z_meso')
      !
   END SUBROUTINE p4z_meso

   SUBROUTINE p4z_meso_init

      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_meso_init  ***
      !!
      !! ** Purpose :   Initialization of M. rubrum parameters
      !!
      !! ** Method  :   Read the nampismes namelist and check the parameters
      !!      called at the first timestep (nittrc000)
      !!
      !! ** input   :   Namelist nampismes
      !!
      !!----------------------------------------------------------------------

      NAMELIST/nampismes/ zz_rate_mesorub_R, zz_rate_mesorub_PicoPredSlope, &
         & zz_rate_mesorub_PicoHalfSat, zz_rate_mesorub_eff, zz_frac_waste_FEN_NH, &
         & zz_frac_waste_FEN_DON, zz_frac_waste_FEN_PON, zz_frac_waste_FEN_BSi 
      INTEGER :: ios           ! Local integer output status for namelist read

      REWIND( numnatp_ref )    ! Namelist nampismes in namelist_ref :  M. rubrum grazing
      READ  ( numnatp_ref, nampismes, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampismes in reference namelist', lwp )

      REWIND( numnatp_cfg )    ! Namelist nampismes in namelist_cfg : M. rubrum grazing
      READ  ( numnatp_cfg, nampismes, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampismes in configuration namelist', lwp )
      IF(lwm) WRITE ( numonp, nampismes )


      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' ' 
         WRITE(numout,*) ' Namelist parameters for M. rubrum grazing, nampismes'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '                zz_rate_mesorub_R =', zz_rate_mesorub_R
         WRITE(numout,*) '    zz_rate_mesorub_PicoPredSlope =', zz_rate_mesorub_PicoPredSlope
         WRITE(numout,*) '      zz_rate_mesorub_PicoHalfSat =', zz_rate_mesorub_PicoHalfSat
         WRITE(numout,*) '              zz_rate_mesorub_eff =', zz_rate_mesorub_eff
         WRITE(numout,*) '             zz_frac_waste_FEN_NH =', zz_frac_waste_FEN_NH
         WRITE(numout,*) '            zz_frac_waste_FEN_DON =', zz_frac_waste_FEN_DON
         WRITE(numout,*) '            zz_frac_waste_FEN_PON =', zz_frac_waste_FEN_PON
         WRITE(numout,*) '            zz_frac_waste_FEN_BSi =', zz_frac_waste_FEN_BSi
      ENDIF


   END SUBROUTINE p4z_meso_init


#else
   !!======================================================================
   !!  Dummy module :                                   No bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_meso                    ! Empty routine
   END SUBROUTINE p4z_meso
#endif 
   !!======================================================================
END MODULE p4zmeso
