MODULE p4zprod
   !!======================================================================
   !!                         ***  MODULE p4zprod  ***
   !! TOP : SMELT  Growth Rates - autotrophic growth 
   !!======================================================================
   !! History :  2015 (E. Olson) adapted from 1-d SOG and PISCES
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   'key_pisces'                              SMELT / PISCES bio-model
   !!----------------------------------------------------------------------
   !!   p4z_prod       :   Compute the growth Rate of the two phytoplanktons groups
   !!   p4z_prod_init  :   Initialization of the parameters for growth
   !!   p4z_prod_alloc :   Allocate variables for growth
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  Source Minus Sink variables
   USE p4zopt          !  optical model
   USE p4zlim          !  Co-limitations of differents nutrients
   USE prtctl_trc      !  print control for debugging
   USE iom             !  I/O manager
#if defined key_skog
   USE p4zcar 
   USE p4zflx  
#endif

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_prod         ! called in p4zbio.F90
   PUBLIC   p4z_prod_init    ! called in trcsms_pisces.F90
   PUBLIC   p4z_prod_alloc

   !! * Shared module variables
         ! SOG params:
      REAL(wp), PUBLIC  ::   zz_rate_R_diat              !: 1/s
      REAL(wp), PUBLIC  ::   zz_rate_R_myri              !: 1/s
      REAL(wp), PUBLIC  ::   zz_rate_R_nano              !: 1/s
      REAL(wp), PUBLIC  ::   zz_rate_maxtemp_diat        !: deg C
      REAL(wp), PUBLIC  ::   zz_rate_maxtemp_myri        !: deg C
      REAL(wp), PUBLIC  ::   zz_rate_maxtemp_nano        !: deg C
      REAL(wp), PUBLIC  ::   zz_rate_temprange_diat      !: deg C
      REAL(wp), PUBLIC  ::   zz_rate_temprange_myri      !: deg C
      REAL(wp), PUBLIC  ::   zz_rate_temprange_nano      !: deg C
      REAL(wp), PUBLIC  ::   zz_rate_Iopt_diat           !: W/m^2
      REAL(wp), PUBLIC  ::   zz_rate_Iopt_myri           !: W/m^2
      REAL(wp), PUBLIC  ::   zz_rate_Iopt_nano           !: W/m^2
      REAL(wp), PUBLIC  ::   zz_rate_gamma_diat          !:
      REAL(wp), PUBLIC  ::   zz_rate_gamma_myri          !:
      REAL(wp), PUBLIC  ::   zz_rate_gamma_nano          !:
      REAL(wp), PUBLIC  ::   zz_rate_K_Si_diat           !:
      REAL(wp), PUBLIC  ::   zz_rate_K_Si_myri           !:
      REAL(wp), PUBLIC  ::   zz_rate_K_Si_nano           !:
      REAL(wp), PUBLIC  ::   zz_rate_kapa_diat           !:
      REAL(wp), PUBLIC  ::   zz_rate_kapa_myri           !:
      REAL(wp), PUBLIC  ::   zz_rate_kapa_nano           !:
      REAL(wp), PUBLIC  ::   zz_rate_k_diat              !:
      REAL(wp), PUBLIC  ::   zz_rate_k_myri              !:
      REAL(wp), PUBLIC  ::   zz_rate_k_nano              !:
      REAL(wp), PUBLIC  ::   zz_rate_Si_ratio_diat       !:
      REAL(wp), PUBLIC  ::   zz_rate_Si_ratio_myri       !:
      REAL(wp), PUBLIC  ::   zz_rate_Si_ratio_nano       !:
      REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   zz_micro_Nlimit   !: diatom N limitation for p4zsink 
      REAL(wp), PUBLIC  ::   zz_chlr

   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: p4zprod.F90 3160 2011-11-20 14:27:18Z cetlod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_prod( kt , knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_prod  ***
      !!
      !! ** Purpose :   Compute PP depending on
      !!              light, temperature and nutrient availability
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT(in) :: kt, knt
      !
      INTEGER  ::   ji, jj, jk
      
      REAL(wp) ::   zproreg, zproreg2, zproregm
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zprorca, zprorcad, zprorcam, zpronew, zpronewd, zpronewm

      REAL(wp), POINTER, DIMENSION(:,:,:) ::   zz_Si, zz_NH, zz_NO, zz_temp, &
         zz_plank_growth_diat, zz_plank_growth_myri, zz_plank_growth_nano, &
         zz_P_diat, zz_P_myri, zz_P_nano, &
         zz_uptake_NO_diat, zz_uptake_NO_myri, zz_uptake_NO_nano, &
         zz_uptake_NH_diat, zz_uptake_NH_myri, zz_uptake_NH_nano, &
         zz_uptake_PC_diat, zz_uptake_PC_myri, zz_uptake_PC_nano, &
         zz_I_par_diat, zz_I_par_myri, zz_I_par_nano, zz_dummy

      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_prod')
      !
      !  Allocate temporary workspace
      CALL wrk_alloc( jpi, jpj, jpk, zprorca, zprorcad, zprorcam, zpronew, zpronewd, zpronewm )
      CALL wrk_alloc( jpi, jpj, jpk, zz_Si, zz_NH, zz_NO, zz_temp, zz_plank_growth_diat, zz_plank_growth_myri, zz_plank_growth_nano)
      CALL wrk_alloc( jpi, jpj, jpk, zz_P_diat, zz_P_myri, zz_P_nano, zz_uptake_NO_diat, zz_uptake_NO_myri )
      CALL wrk_alloc( jpi, jpj, jpk, zz_uptake_NO_nano, zz_I_par_diat, zz_I_par_myri, zz_I_par_nano, zz_dummy )
      CALL wrk_alloc( jpi, jpj, jpk, zz_uptake_NH_diat, zz_uptake_PC_diat, zz_uptake_NH_myri)
      CALL wrk_alloc( jpi, jpj, jpk, zz_uptake_PC_myri, zz_uptake_NH_nano, zz_uptake_PC_nano)     
      !
      zprorca (:,:,:) = 0._wp
      zprorcad(:,:,:) = 0._wp
      zprorcam(:,:,:) = 0._wp
      zpronew (:,:,:) = 0._wp
      zpronewd(:,:,:) = 0._wp
      zpronewm(:,:,:) = 0._wp
      zz_plank_growth_diat(:,:,:) = 0._wp
      zz_plank_growth_myri(:,:,:) = 0._wp
      zz_plank_growth_nano(:,:,:) = 0._wp
      zz_uptake_NO_diat(:,:,:) = 0._wp
      zz_uptake_NO_myri(:,:,:) = 0._wp
      zz_uptake_NO_nano(:,:,:) = 0._wp

      ! set SMELT vars:
      zz_P_diat(:,:,:)  = trb(:,:,:,jpdia)
      zz_P_myri(:,:,:)  = trb(:,:,:,jpmes)
      zz_P_nano(:,:,:) = trb(:,:,:,jpphy)
      zz_I_par_diat(:,:,:) = ediat(:,:,:)
      zz_I_par_myri(:,:,:) = ediat(:,:,:)
      zz_I_par_nano(:,:,:) = ediat(:,:,:)
      zz_Si(:,:,:) = trb(:,:,:,jpsil)
      zz_NH(:,:,:) = trb(:,:,:,jpnh4)
      zz_NO(:,:,:) = trb(:,:,:,jpno3)
      zz_temp(:,:,:) = tsn(:,:,:,jp_tem)

      !diatoms:
      CALL p4z_growthSOG(zz_NO, zz_NH, zz_Si, zz_P_diat, zz_I_par_diat, zz_temp, zz_rate_R_diat, zz_rate_maxtemp_diat, &
        zz_rate_temprange_diat, zz_rate_Iopt_diat, zz_rate_gamma_diat, zz_rate_K_Si_diat, &
        zz_rate_kapa_diat, zz_rate_k_diat, zz_rate_Si_ratio_diat, zz_plank_growth_diat, zz_uptake_NO_diat, zz_uptake_NH_diat, zz_uptake_PC_diat, zz_micro_Nlimit)
      !M. rubrum:
      CALL p4z_growthSOG(zz_NO, zz_NH, zz_Si, zz_P_myri, zz_I_par_myri, zz_temp, zz_rate_R_myri, zz_rate_maxtemp_myri, &
        zz_rate_temprange_myri, zz_rate_Iopt_myri, zz_rate_gamma_myri, zz_rate_K_Si_myri, &
        zz_rate_kapa_myri, zz_rate_k_myri, zz_rate_Si_ratio_myri, zz_plank_growth_myri, zz_uptake_NO_myri, zz_uptake_NH_myri, zz_uptake_PC_myri, zz_dummy)
      !nano/pico:
      CALL p4z_growthSOG(zz_NO, zz_NH, zz_Si, zz_P_nano, zz_I_par_nano, zz_temp, zz_rate_R_nano, zz_rate_maxtemp_nano, &
        zz_rate_temprange_nano, zz_rate_Iopt_nano, zz_rate_gamma_nano, zz_rate_K_Si_nano, &
        zz_rate_kapa_nano, zz_rate_k_nano, zz_rate_Si_ratio_nano, zz_plank_growth_nano, zz_uptake_NO_nano, zz_uptake_NH_nano, zz_uptake_PC_nano, zz_dummy)

#if defined key_skog      
      CALL p4z_car( kt, knt, zz_uptake_NO_diat, zz_uptake_NO_nano, zz_uptake_NO_myri, &
        zz_uptake_NH_diat, zz_uptake_NH_nano, zz_uptake_NH_myri, & 
        zz_uptake_PC_diat, zz_uptake_PC_nano, zz_uptake_PC_myri)
  
      CALL p4z_ta( kt, knt, zz_uptake_NO_diat, zz_uptake_NO_nano, zz_uptake_NO_myri, &
        zz_uptake_NH_diat, zz_uptake_NH_nano, zz_uptake_NH_myri )    
      
      CALL p4z_oxy( kt, knt, zz_uptake_NO_diat, zz_uptake_NO_nano, zz_uptake_NO_myri, &
        zz_uptake_NH_diat, zz_uptake_NH_nano, zz_uptake_NH_myri)

      CALL p4z_flx( kt, knt)
#endif

      ! Computation of the various production terms
!CDIR NOVERRCHK
      DO jk = 1, jpkm1
!CDIR NOVERRCHK
         DO jj = 1, jpj
!CDIR NOVERRCHK
            DO ji = 1, jpi
               !  production terms for nanophyto. zprbio=zz_plank_growth_nano
               zprorca(ji,jj,jk) = zz_plank_growth_nano(ji,jj,jk)  * zz_P_nano(ji,jj,jk) * rfact2
               zpronew(ji,jj,jk) = zz_uptake_NO_nano(ji,jj,jk) * rfact2

               zprorcam(ji,jj,jk) = zz_plank_growth_myri(ji,jj,jk) * zz_P_myri(ji,jj,jk) * rfact2
               zpronewm(ji,jj,jk) = zz_uptake_NO_myri(ji,jj,jk) *rfact2
                 
               !  production terms for diatomees
               zprorcad(ji,jj,jk) = zz_plank_growth_diat(ji,jj,jk) * zz_P_diat(ji,jj,jk) * rfact2
               zpronewd(ji,jj,jk) = zz_uptake_NO_diat(ji,jj,jk) * rfact2
            END DO
         END DO
      END DO

      CALL iom_put( "PAR",ediat) ! PAR (as used in bio model)
      CALL iom_put( "PPDIAT",zz_plank_growth_diat * zz_P_diat) ! Diatom primary productivity (umol N/s)
      CALL iom_put( "PPPHY",zz_plank_growth_nano * zz_P_nano) ! Small phyto primary productivity (umol N/s)
      CALL iom_put( "PPMRUB",zz_plank_growth_myri * zz_P_myri) ! M. Rubrum primary productivity (umol N/s)
      CALL iom_put( "PPDIATNO3",zz_uptake_NO_diat) ! Diatom primary productivity (umol N/s)
      CALL iom_put( "PPPHYNO3",zz_uptake_NO_nano) ! Small phyto primary productivity (umol N/s)
      CALL iom_put( "PPMRUBNO3",zz_uptake_NO_myri) ! M. Rubrum primary productivity (umol N/s)

      !   Update the arrays TRA which contain the biological sources and sinks
      DO jk = 1, jpkm1
         DO jj = 1, jpj
           DO ji =1 ,jpi
              zproreg  = zprorca(ji,jj,jk) - zpronew(ji,jj,jk)
              zproreg2 = zprorcad(ji,jj,jk) - zpronewd(ji,jj,jk)
              zproregm = zprorcam(ji,jj,jk) - zpronewm(ji,jj,jk)
              tra(ji,jj,jk,jpno3) = tra(ji,jj,jk,jpno3) - zpronew(ji,jj,jk) - zpronewd(ji,jj,jk) - zpronewm(ji,jj,jk)
              tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) - zproreg - zproreg2 - zproregm
              tra(ji,jj,jk,jpphy) = tra(ji,jj,jk,jpphy) + zprorca(ji,jj,jk)
              tra(ji,jj,jk,jpmes) = tra(ji,jj,jk,jpmes) + zprorcam(ji,jj,jk)
              tra(ji,jj,jk,jpdia) = tra(ji,jj,jk,jpdia) + zprorcad(ji,jj,jk)
              tra(ji,jj,jk,jpsil) = tra(ji,jj,jk,jpsil) - (zprorcad(ji,jj,jk)*zz_rate_Si_ratio_diat + zprorca(ji,jj,jk)*zz_rate_Si_ratio_nano + zprorcam(ji,jj,jk)*zz_rate_Si_ratio_myri)
          END DO
        END DO
     END DO

      !
      CALL wrk_dealloc( jpi, jpj, jpk, zprorca, zprorcad, zprorcam, zpronew, zpronewd, zpronewm )
      CALL wrk_dealloc( jpi, jpj, jpk, zz_Si, zz_NH, zz_NO, zz_temp, zz_plank_growth_diat, zz_plank_growth_myri, zz_plank_growth_nano)
      CALL wrk_dealloc( jpi, jpj, jpk,  zz_P_diat, zz_P_myri, zz_P_nano, zz_uptake_NO_diat, zz_uptake_NO_myri )
      CALL wrk_dealloc( jpi, jpj, jpk,  zz_uptake_NO_nano, zz_I_par_diat, zz_I_par_myri, zz_I_par_nano, zz_dummy )
      CALL wrk_dealloc(jpi, jpj, jpk, zz_uptake_NH_diat, zz_uptake_PC_diat, zz_uptake_NH_myri, zz_uptake_PC_myri, zz_uptake_NH_nano, zz_uptake_PC_nano)

     IF( nn_timing == 1 )  CALL timing_stop('p4z_prod')
     !
   END SUBROUTINE p4z_prod

   SUBROUTINE p4z_growthSOG(zz_NO, zz_NH, zz_Si, zz_P, zz_I_par, zz_temp, zz_rate_R, zz_rate_maxtemp, &
        zz_rate_temprange, zz_rate_Iopt, zz_rate_gamma, zz_rate_K_Si, &
        zz_rate_kapa, zz_rate_k, zz_rate_Si_ratio, zz_plank_growth, zz_uptake_NO, &
        zz_uptake_NH, zz_uptake_PC, zz_plank_Nlimit)
     !! calculate growth rates and nutrient utilization
        !! light and nutrient limitation
     INTEGER  ::   ji, jj, jk
     REAL(wp), INTENT(in)   :: zz_NO(:,:,:), zz_NH(:,:,:), zz_Si(:,:,:), zz_P(:,:,:), &
        zz_I_par(:,:,:), zz_temp(:,:,:)        ! 3d input arrays, DIMENSION(jpi,jpj,jpk)
     REAL(wp), INTENT(in)  :: zz_rate_R, zz_rate_maxtemp, zz_rate_temprange, zz_rate_Iopt, & 
        zz_rate_gamma, zz_rate_K_Si, zz_rate_kapa, zz_rate_k, &
        zz_rate_Si_ratio  ! parameters
     REAL(wp), INTENT(out), DIMENSION(jpi,jpj,jpk)   :: zz_plank_growth, zz_uptake_NO, &
        zz_uptake_NH, zz_uptake_PC, zz_plank_Nlimit ! 3d output arrays
   
     REAL(wp), POINTER, DIMENSION(:,:,:) ::   zz_Rmax(:,:,:), zz_plank_growth_light(:,:,:), zz_Uc(:,:,:), &
        zz_Sc(:,:,:), zz_Oup_cell(:,:,:), zz_Hup_cell(:,:,:)
   
         ! Allocate temporary workspace vars
     CALL wrk_alloc( jpi, jpj, jpk, zz_Rmax, zz_plank_growth_light, zz_Uc, zz_Sc, zz_Oup_cell, zz_Hup_cell)
        
   ! Q10 effect SOG; replaced temp_Q10 with tgfunc and
      ! KtoC(temp(j)) with tsn(ji,jj,jk,jp_tem)
      ! diatoms:
      zz_Rmax(:,:,:)= zz_rate_R * tgfunc(:,:,:) &
            * min(max(zz_rate_maxtemp - zz_temp, 0.0_wp), &
                  zz_rate_temprange) / &
            (zz_rate_temprange + epsilon(zz_rate_temprange))

      
   ! Computation of SOG limitation terms: light, nutrients
      DO jk = 1, jpkm1
         DO jj = 1, jpj
           DO ji =1, jpi
                  ! LIGHT
                  zz_plank_growth_light(ji,jj,jk) = &
                       ! Steeles scheme like but with extended high light range
                       ! (Steeles scheme has built in light inhibition)
                       ! much broader pattern for a bucket of organisms
                       (1.0_wp - exp(-zz_I_par(ji,jj,jk) / (0.33_wp * zz_rate_Iopt)) ) * &
                       (exp(-zz_I_par(ji,jj,jk) / (30._wp * zz_rate_Iopt))) * 1.06_wp
                  zz_Uc(ji,jj,jk) = (1.0_wp - zz_rate_gamma) * zz_plank_growth_light(ji,jj,jk)

                  ! Si
                  zz_Sc(ji,jj,jk) = zz_Si(ji,jj,jk) / (zz_rate_K_Si + zz_Si(ji,jj,jk)+epsilon(zz_rate_K_Si))

                  ! Nitrate and Ammonium
                  IF (zz_NO(ji,jj,jk) > epsilon(zz_NO(ji,jj,jk))) THEN
                     zz_Oup_cell(ji,jj,jk) = zz_NO(ji,jj,jk) * zz_rate_kapa / &
                          (zz_rate_k + zz_NO(ji,jj,jk) * zz_rate_kapa + &
                          zz_NH(ji,jj,jk)+epsilon(zz_NH))
                  ELSE
                     zz_Oup_cell(ji,jj,jk) = 0._wp
                  ENDIF
                  IF (zz_NH(ji,jj,jk) > epsilon(zz_NH(ji,jj,jk))) THEN
                     zz_Hup_cell(ji,jj,jk) = zz_NH(ji,jj,jk) / &
                          (zz_rate_k + zz_NO(ji,jj,jk) * zz_rate_kapa + &
                          zz_NH(ji,jj,jk))
                  ELSE
                     zz_Hup_cell(ji,jj,jk) = 0._wp
                  ENDIF

                  IF (zz_Oup_cell(ji,jj,jk) < 0._wp) THEN
                     WRITE(numout,*) "Oup_cell(ji,jj,jk) < 0. in NPZD.f90"
                     WRITE(numout,*) zz_Oup_cell(ji,jj,jk)
                     CALL EXIT(1)
                  ENDIF
                  IF (zz_Hup_cell(ji,jj,jk) < 0._wp) THEN
                     WRITE(numout,*) "Hup_cell(ji,jj,jk) < 0. in NPZD.f90"
                     WRITE(numout,*) zz_Hup_cell(ji,jj,jk)
                     CALL EXIT(1)
                  ENDIF

                  ! exponent of 1/5 follows Alain
                  zz_plank_Nlimit(ji,jj,jk) = min((zz_Oup_cell(ji,jj,jk) + &
                     zz_Hup_cell(ji,jj,jk)), zz_Sc(ji,jj,jk))**0.2_wp

                  ! Choose light limitation or nutrient limitation
                  IF (zz_Uc(ji,jj,jk) < 0._wp ) THEN
                     zz_plank_growth(ji,jj,jk) = 0._wp
                  ELSE

                     IF (min(zz_Uc(ji,jj,jk),zz_Sc(ji,jj,jk)) >= &
                          zz_Oup_cell(ji,jj,jk) + zz_Hup_cell(ji,jj,jk)) THEN
                        !N LIMITING
                        zz_plank_growth(ji,jj,jk) = zz_Rmax(ji,jj,jk) * &
                             (zz_Oup_cell(ji,jj,jk) + zz_Hup_cell(ji,jj,jk))

                        IF (zz_plank_growth(ji,jj,jk) < 0._wp) THEN
                           zz_plank_growth(ji,jj,jk) = 0._wp
                        ENDIF
                        zz_uptake_NO(ji,jj,jk) = zz_Rmax(ji,jj,jk) * zz_Oup_cell(ji,jj,jk) * &
                             zz_P(ji,jj,jk)
                        zz_uptake_NH(ji,jj,jk) = zz_Rmax(ji,jj,jk) * zz_Hup_cell(ji,jj,jk) * &
                             zz_P(ji,jj,jk)

                     ELSE
                        IF (zz_Uc(ji,jj,jk) < zz_Sc(ji,jj,jk)) THEN
                        !LIGHT LIMITING
                           zz_plank_growth(ji,jj,jk) = zz_Rmax(ji,jj,jk) * zz_Uc(ji,jj,jk)
                        ELSE
                        ! Si limitation
                           zz_plank_growth(ji,jj,jk) = zz_Rmax(ji,jj,jk) * zz_Sc(ji,jj,jk)
                        ENDIF
                        ! split the nitrogen uptake between NH and NO
                        IF (zz_plank_growth(ji,jj,jk) <= &
                             zz_Rmax(ji,jj,jk) * zz_Hup_cell(ji,jj,jk)) THEN
                           zz_uptake_NH(ji,jj,jk) = zz_plank_growth(ji,jj,jk) * &
                                zz_P(ji,jj,jk) 
                           zz_uptake_NO(ji,jj,jk) = 0.0_wp 
                        ELSE
                           zz_uptake_NH(ji,jj,jk) = zz_Rmax(ji,jj,jk) * &
                                zz_Hup_cell(ji,jj,jk) * zz_P(ji,jj,jk)
                           zz_uptake_NO(ji,jj,jk) = (zz_plank_growth(ji,jj,jk) - &
                                zz_Rmax(ji,jj,jk) * zz_Hup_cell(ji,jj,jk)) &
                                * zz_P(ji,jj,jk)

                        ENDIF
                     ENDIF

                     ! PC Differential Carbon Uptake
                     ! Chlr Redux Factor = 0.2 (Ianson and Allen 2002)
                     zz_uptake_PC(ji,jj,jk) = (zz_Rmax(ji,jj,jk) * zz_Uc(ji,jj,jk) - &
                          zz_plank_growth(ji,jj,jk)) &
                          * 0.2_wp * zz_P(ji,jj,jk)
                  ENDIF
           END DO
        END DO
     END DO

     
      ! Deallocate temporary workspace SOG vars 
     CALL wrk_dealloc( jpi, jpj, jpk, zz_Rmax, zz_plank_growth_light, zz_Uc, zz_Sc, zz_Oup_cell, zz_Hup_cell)
   END SUBROUTINE P4z_growthSOG

   SUBROUTINE p4z_prod_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_prod_init  ***
      !!
      !! ** Purpose :   Initialization of PP parameters
      !!
      !! ** Method  :   Read the nampisprod namelist and check the parameters
      !!      called at the first timestep (nittrc000)
      !!
      !! ** input   :   Namelist nampisprod
      !!----------------------------------------------------------------------
      !
      NAMELIST/nampisprod/ zz_rate_R_diat, zz_rate_R_myri, zz_rate_R_nano, zz_rate_maxtemp_diat, &
      &        zz_rate_maxtemp_myri, zz_rate_maxtemp_nano, zz_rate_temprange_diat, zz_rate_temprange_myri, &
      &        zz_rate_temprange_nano, zz_rate_Iopt_diat, zz_rate_Iopt_myri, zz_rate_Iopt_nano, &
      &        zz_rate_gamma_diat, zz_rate_gamma_myri, zz_rate_gamma_nano, zz_rate_K_Si_diat, zz_rate_K_Si_myri, &
      &        zz_rate_K_Si_nano, &
      &        zz_rate_kapa_diat, zz_rate_kapa_myri, zz_rate_kapa_nano, zz_rate_k_diat, zz_rate_k_myri, &
      &        zz_rate_k_nano, &
      &        zz_rate_Si_ratio_diat, zz_rate_Si_ratio_myri, zz_rate_Si_ratio_nano, zz_chlr
      INTEGER :: ios                 ! Local integer output status for namelist read
      !!----------------------------------------------------------------------

      REWIND( numnatp_ref )              ! Namelist nampisprod in reference namelist : SMELT phytoplankton production
      READ  ( numnatp_ref, nampisprod, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampisprod in reference namelist', lwp )

      REWIND( numnatp_cfg )              ! Namelist nampisprod in configuration namelist : SMELT phytoplankton production
      READ  ( numnatp_cfg, nampisprod, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampisprod in configuration namelist', lwp )
      IF(lwm) WRITE ( numonp, nampisprod )

      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for phytoplankton growth, nampisprod'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) ' zz_rate_R_diat =', zz_rate_R_diat
         WRITE(numout,*) ' zz_rate_R_myri =',  zz_rate_R_myri
         WRITE(numout,*) ' zz_rate_R_nano =', zz_rate_R_nano
         WRITE(numout,*) ' zz_rate_maxtemp_diat =', zz_rate_maxtemp_diat
         WRITE(numout,*) ' zz_rate_maxtemp_myri =', zz_rate_maxtemp_myri
         WRITE(numout,*) ' zz_rate_maxtemp_nano =', zz_rate_maxtemp_nano
         WRITE(numout,*) ' zz_rate_temprange_diat =', zz_rate_temprange_diat
         WRITE(numout,*) ' zz_rate_temprange_myri =', zz_rate_temprange_myri
         WRITE(numout,*) ' zz_rate_temprange_nano =', zz_rate_temprange_nano
         WRITE(numout,*) ' zz_rate_Iopt_diat =', zz_rate_Iopt_diat
         WRITE(numout,*) ' zz_rate_Iopt_myri =', zz_rate_Iopt_myri
         WRITE(numout,*) ' zz_rate_Iopt_nano =', zz_rate_Iopt_nano
         WRITE(numout,*) ' zz_rate_gamma_diat =', zz_rate_gamma_diat
         WRITE(numout,*) ' zz_rate_gamma_myri =', zz_rate_gamma_myri
         WRITE(numout,*) ' zz_rate_gamma_nano =', zz_rate_gamma_nano
         WRITE(numout,*) ' zz_rate_K_Si_diat =', zz_rate_K_Si_diat
         WRITE(numout,*) ' zz_rate_K_Si_myri =', zz_rate_K_Si_myri
         WRITE(numout,*) ' zz_rate_K_Si_nano =', zz_rate_K_Si_nano
         WRITE(numout,*) ' zz_rate_kapa_diat =', zz_rate_kapa_diat
         WRITE(numout,*) ' zz_rate_kapa_myri =', zz_rate_kapa_myri
         WRITE(numout,*) ' zz_rate_kapa_nano =', zz_rate_kapa_nano
         WRITE(numout,*) ' zz_rate_k_diat =', zz_rate_k_diat
         WRITE(numout,*) ' zz_rate_k_myri =', zz_rate_k_myri
         WRITE(numout,*) ' zz_rate_k_nano =', zz_rate_k_nano
         WRITE(numout,*) ' zz_rate_Si_ratio_diat =', zz_rate_Si_ratio_diat
         WRITE(numout,*) ' zz_rate_Si_ratio_myri =', zz_rate_Si_ratio_myri
         WRITE(numout,*) ' zz_rate_Si_ratio_nano =', zz_rate_Si_ratio_nano
         WRITE(numout,*) ' zz_chlr =', zz_chlr
      ENDIF

   END SUBROUTINE p4z_prod_init


   INTEGER FUNCTION p4z_prod_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_prod_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( zz_micro_Nlimit(jpi,jpj,jpk), STAT = p4z_prod_alloc )
      !
      IF( p4z_prod_alloc /= 0 ) CALL ctl_warn('p4z_prod_alloc : failed to allocate arrays.')
      !
   END FUNCTION p4z_prod_alloc



#else
   !!======================================================================
   !!  Dummy module :                                   No bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_prod                    ! Empty routine
   END SUBROUTINE p4z_prod
#endif 
   !!======================================================================
END MODULE p4zprod
