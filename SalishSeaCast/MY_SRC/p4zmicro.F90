MODULE p4zmicro
   !!======================================================================
   !!                         ***  MODULE p4zmicro  ***
   !! TOP :   SMELT Compute the sources/sinks for microzooplankton
   !!======================================================================
   !! History :   2015 (E. Olson) adapted from 1-d SOG and PISCES
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   'key_pisces'                               SMELT / PISCES bio-model
   !!----------------------------------------------------------------------
   !!   p4z_micro       :   Compute the sources/sinks for microzooplankton
   !!   p4z_micro_init  :   Initialize and read the appropriate namelist
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  SMELT Source Minus Sink variables
   USE p4zint          !  interpolation and computation of various fields
   USE iom             !  I/O manager
   USE prtctl_trc      !  print control for debugging
   USE p4zprod, ONLY:  zz_rate_Si_ratio_diat, zz_rate_Si_ratio_myri, zz_rate_Si_ratio_nano

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_micro         ! called in p4zbio.F90
   PUBLIC   p4z_micro_init    ! called in trcsms_pisces.F90

   !! * Shared module variables
   REAL(wp) ::   zz_rate_uzoo_Rm                 ! uM N / s microzo natural mortality rate @ 20 deg C
   REAL(wp) ::   zz_rate_uzoo_excr               ! uM N / s microzo excretion rate @ 20 deg C
   REAL(wp) ::   zz_frac_waste_ZNM_NH            
   REAL(wp) ::   zz_frac_waste_ZEX_NH            
   REAL(wp) ::   zz_rate_uzoo_PredSlope          ! uM N microzo total grazing limit
   REAL(wp) ::   zz_rate_uzoo_HalfSat            ! uM N microzo total grazing half sat const
   REAL(wp) ::   zz_rate_uzoo_MicroPref          ! microzo preference for micro-phyto
   REAL(wp) ::   zz_rate_uzoo_MicroPredslope     ! uM N microzo micro-phyto grazing half sat const
   REAL(wp) ::   zz_rate_uzoo_MicroHalfSat       ! uM N microzo micro-phyto grazing half sat const
   REAL(wp) ::   zz_rate_uzoo_NanoPref           ! microzo preference for nano-phyto
   REAL(wp) ::   zz_rate_uzoo_NanoPredslope      ! uM N microzo nano-phyto grazing limit
   REAL(wp) ::   zz_rate_uzoo_NanoHalfSat        ! uM N microzo nano-phyto grazing half sat const
   REAL(wp) ::   zz_rate_uzoo_PicoPref           ! microzo preference for pico-phyto
   REAL(wp) ::   zz_rate_uzoo_PicoPredslope      ! uM N microzo pico-phyto grazing limit
   REAL(wp) ::   zz_rate_uzoo_PicoHalfSat        ! uM N microzo pico-phyto grazing half sat const
   REAL(wp) ::   zz_rate_uzoo_PON_Pref           ! microzo preference for PON
   REAL(wp) ::   zz_rate_uzoo_PON_Predslope      ! uM N microzo microzoo grazing limit
   REAL(wp) ::   zz_rate_uzoo_PON_HalfSat        ! uM N microzo PON grazing half sat const
   REAL(wp) ::   zz_rate_uzoo_Z_Pref             ! microzo preference for microzoo cannibalism
   REAL(wp) ::   zz_rate_uzoo_Z_Predslope        ! uM N microzo microzoo grazing limit
   REAL(wp) ::   zz_rate_uzoo_Z_HalfSat          ! uM N microzo microzoo grazing half sat const
   REAL(wp) ::   zz_rate_uzoo_R                  ! match SOG 0.600E-04_wp       ! 1/s micro-phyto maximum growth rate # Hitchcock 1980, 1.4 d-1 for T. nordenskelii at 10degrees
   REAL(wp) ::   zz_frac_waste_PEZ_NH            ! waste frac from microzoo grazing PON to NH
   REAL(wp) ::   zz_frac_waste_DEZ_NH            ! waste frac from microzoo grazing microphyto to NH
   REAL(wp) ::   zz_frac_waste_NEZ_NH            ! waste frac from microzoo grazing nanophyto to NH
   REAL(wp) ::   zz_frac_waste_FEZ_NH            ! waste frac from microzoo grazing picophyto to NH
   REAL(wp) ::   zz_frac_waste_ZEZ_NH            ! match SOG 1._wp               ! waste frac from micro-zoo excretion to NH
   REAL(wp) ::   zz_rate_uZoo_eff                ! match SOG
   REAL(wp) ::   zz_frac_waste_ZNM_DON           ! waste frac from micro-zoo natural mortality to DON
   REAL(wp) ::   zz_frac_waste_ZEX_DON           ! waste frac from micro-zoo excretion to DON
   REAL(wp) ::   zz_frac_waste_DEZ_DON           ! waste frac from microzoo grazing microphyto to DON
   REAL(wp) ::   zz_frac_waste_FEZ_DON           ! waste frac from microzoo grazing picophyto to DON
   REAL(wp) ::   zz_frac_waste_NEZ_DON           ! waste frac from microzoo grazing nanophyto to DON
   REAL(wp) ::   zz_frac_waste_PEZ_DON           ! waste frac from microzoo grazing PON to DON
   REAL(wp) ::   zz_frac_waste_ZEZ_DON           ! waste frac from microzoo grazing microzoo to DON
   REAL(wp) ::   zz_frac_waste_ZNM_PON           ! waste frac from micro-zoo natural mortality to PON
   REAL(wp) ::   zz_frac_waste_ZEX_PON           ! waste frac from micro-zoo excretion to PON
   REAL(wp) ::   zz_frac_waste_DEZ_PON           ! waste frac from microzoo grazing microphyto to PON
   REAL(wp) ::   zz_frac_waste_FEZ_PON           ! waste frac from microzoo grazing picophyto to PON
   REAL(wp) ::   zz_frac_waste_NEZ_PON           ! waste frac from microzoo grazing nanophyto to PON
   REAL(wp) ::   zz_frac_waste_PEZ_PON           ! waste frac from microzoo grazing PON to PON
   REAL(wp) ::   zz_frac_waste_ZEZ_PON           ! waste frac from microzoo grazing microzoo to PON
   REAL(wp) ::   zz_frac_waste_ZNM_BSi           ! waste frac from micro-zoo natural mortality to bSi
   REAL(wp) ::   zz_frac_waste_ZEX_Bsi           ! waste frac from micro-zoo excretion to bSi
   REAL(wp) ::   zz_frac_waste_DEZ_Bsi           ! waste frac from microzoo grazing microphyto to Bsi
   REAL(wp) ::   zz_frac_waste_FEZ_BSi           ! waste frac from microzoo grazing picophyto to Bsi
   REAL(wp) ::   zz_frac_waste_NEZ_BSi           ! waste frac from microzoo grazing nanophyto to Bsi
   REAL(wp) ::   zz_frac_waste_PEZ_BSi           ! waste frac from microzoo grazing PON to Bsi
   REAL(wp) ::   zz_frac_waste_ZEZ_BSi           ! waste frac from microzoo grazing microzoo to Bsi
   
   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: p4zmicro.F90 3160 2011-11-20 14:27:18Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE p4z_micro( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_micro  ***
      !!
      !! ** Purpose :   Compute the sources/sinks for microzooplankton
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::  kt  ! ocean time step
      INTEGER, INTENT(in) ::  knt 
      INTEGER  :: ji, jj, jk
      REAL(wp) ::   zz_NatMort_uzoo, zz_Excr_uzoo, zz_food_limitation, zz_denominator, &
        zz_MicroZoo_eat
      REAL(wp), POINTER, DIMENSION(:,:,:) ::  zz_P_diat, zz_P_myri, zz_P_nano, zz_Z, zz_D_PON, zz_uZoo_mort_micro, &
        zz_uZoo_mort_nano, zz_uZoo_mort_pico, zz_uZoo_graz_PON, zz_uZoo_graz_Z
       
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_micro')
      CALL wrk_alloc( jpi, jpj, jpk, zz_P_diat, zz_P_myri, zz_P_nano, zz_Z, zz_D_PON, zz_uZoo_mort_micro )
      CALL wrk_alloc( jpi, jpj, jpk, zz_uZoo_mort_nano, zz_uZoo_mort_pico, zz_uZoo_graz_PON, zz_uZoo_graz_Z)
      
      zz_P_diat(:,:,:) = trb(:,:,:,jpdia)
      zz_P_myri(:,:,:) = trb(:,:,:,jpmes)
      zz_P_nano(:,:,:) = trb(:,:,:,jpphy)
      zz_Z(:,:,:)      = trb(:,:,:,jpzoo)
      zz_D_PON(:,:,:)  = trb(:,:,:,jppon)

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi

                ! microzoo natural mortality:
                zz_NatMort_uzoo = zz_rate_uzoo_Rm * tgfunc(ji,jj,jk) * zz_Z(ji,jj,jk)
                zz_Excr_uzoo = zz_rate_uzoo_excr * tgfunc(ji,jj,jk) * zz_Z(ji,jj,jk)

                tra(ji,jj,jk,jpzoo) = tra(ji,jj,jk,jpzoo) - (zz_NatMort_uzoo + zz_Excr_uzoo) * rfact2
                tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) + (zz_frac_waste_ZNM_NH * zz_NatMort_uzoo &
                    + zz_frac_waste_ZEX_NH * zz_Excr_uzoo) * rfact2
                    
                tra(ji,jj,jk,jpdon) = tra(ji,jj,jk,jpdon) + (zz_frac_waste_ZNM_DON * zz_NatMort_uzoo &
                    + zz_frac_waste_ZEX_DON * zz_Excr_uzoo) * rfact2
                tra(ji,jj,jk,jppon) = tra(ji,jj,jk,jppon) + (zz_frac_waste_ZNM_PON * zz_NatMort_uzoo &
                    + zz_frac_waste_ZEX_PON * zz_Excr_uzoo) * rfact2
                !tra(ji,jj,jk,(Ref)) = tra(ji,jj,jk,(Ref)) + zz_frac_waste_ZNM_Ref * zz_NatMort_uzoo &
                !    + zz_frac_waste_ZEX_Ref * zz_Excr_uzoo
                ! no transfer to BSi

                ! microzoo grazing:
                zz_food_limitation = (zz_P_diat(ji,jj,jk) + zz_P_myri(ji,jj,jk) + zz_P_nano(ji,jj,jk) + &
                    zz_D_PON(ji,jj,jk) + zz_Z(ji,jj,jk) - zz_rate_uzoo_PredSlope) / &
                    (zz_rate_uzoo_HalfSat + zz_P_diat(ji,jj,jk) + zz_P_myri(ji,jj,jk) + zz_P_nano(ji,jj,jk) &
                    + zz_D_PON(ji,jj,jk) + zz_Z(ji,jj,jk) - zz_rate_uzoo_Predslope &
                    + epsilon(zz_rate_uzoo_HalfSat))

                zz_denominator = (zz_rate_uzoo_MicroPref * zz_P_diat(ji,jj,jk) + &
                    zz_rate_uzoo_NanoPref * zz_P_myri(ji,jj,jk) + &
                    zz_rate_uzoo_PicoPref * zz_P_nano(ji,jj,jk) + &
                    zz_rate_uzoo_PON_Pref * zz_D_PON(ji,jj,jk) + &
                    zz_rate_uzoo_Z_Pref * zz_Z(ji,jj,jk) + &
                    epsilon(zz_P_diat(ji,jj,jk)) )

                zz_uZoo_mort_micro(ji,jj,jk) = min(zz_rate_uzoo_MicroPref * zz_food_limitation &
                    * zz_P_diat(ji,jj,jk) / zz_denominator, &
                    (zz_P_diat(ji,jj,jk) - zz_rate_uzoo_MicroPredslope) / &
                    (zz_rate_uzoo_MicroHalfSat + zz_P_diat(ji,jj,jk) &
                    - zz_rate_uzoo_MicroPredslope + &
                    epsilon(zz_rate_uzoo_MicroHalfSat)) )

                zz_uZoo_mort_nano(ji,jj,jk) = min(zz_rate_uzoo_NanoPref * zz_food_limitation &
                    * zz_P_myri(ji,jj,jk) / zz_denominator, &
                    (zz_P_myri(ji,jj,jk) - zz_rate_uzoo_NanoPredslope) / &
                    (zz_rate_uzoo_NanoHalfSat + zz_P_myri(ji,jj,jk) &
                    - zz_rate_uzoo_NanoPredslope + &
                    epsilon(zz_rate_uzoo_NanoHalfSat)) )

                 zz_uZoo_mort_pico(ji,jj,jk) = min(zz_rate_uzoo_PicoPref * zz_food_limitation &
                    * zz_P_nano(ji,jj,jk) / zz_denominator, &
                    (zz_P_nano(ji,jj,jk) - zz_rate_uzoo_PicoPredslope) / &
                    (zz_rate_uzoo_PicoHalfSat + zz_P_nano(ji,jj,jk) &
                    - zz_rate_uzoo_PicoPredslope + &
                    epsilon(zz_rate_uzoo_PicoHalfSat)) )

                zz_uZoo_graz_PON(ji,jj,jk) = min(zz_rate_uzoo_PON_Pref * zz_food_limitation &
                    * zz_D_PON(ji,jj,jk) / zz_denominator, &
                    (zz_D_PON(ji,jj,jk) - zz_rate_uzoo_PON_Predslope) / &
                    (zz_rate_uzoo_PON_HalfSat + zz_D_PON(ji,jj,jk) &
                    - zz_rate_uzoo_PON_Predslope + &
                    epsilon(zz_rate_uzoo_PON_HalfSat)) )

                zz_uZoo_graz_Z(ji,jj,jk) = min(zz_rate_uzoo_Z_Pref * zz_food_limitation &
                    * zz_Z(ji,jj,jk) / zz_denominator, &
                    (zz_Z(ji,jj,jk) - zz_rate_uzoo_Z_Predslope) / &
                    (zz_rate_uzoo_Z_HalfSat + zz_Z(ji,jj,jk) &
                    - zz_rate_uzoo_Z_Predslope + &
                    epsilon(zz_rate_uzoo_Z_HalfSat)) )

                zz_uZoo_mort_micro(ji,jj,jk) = zz_rate_uzoo_R * tgfunc(ji,jj,jk) &
                    * zz_Z(ji,jj,jk) * max(0._wp, zz_uZoo_mort_micro(ji,jj,jk))
                zz_uZoo_mort_nano(ji,jj,jk) = zz_rate_uzoo_R * tgfunc(ji,jj,jk) &
                    * zz_Z(ji,jj,jk) * max(0._wp, zz_uZoo_mort_nano(ji,jj,jk))
                zz_uZoo_mort_pico(ji,jj,jk) = zz_rate_uzoo_R * tgfunc(ji,jj,jk) &
                    * zz_Z(ji,jj,jk) * max(0._wp, zz_uZoo_mort_pico(ji,jj,jk))
                zz_uZoo_graz_PON(ji,jj,jk) = zz_rate_uzoo_R * tgfunc(ji,jj,jk) &
                    * zz_Z(ji,jj,jk) * max(0._wp, zz_uZoo_graz_PON(ji,jj,jk))
                zz_uZoo_graz_Z(ji,jj,jk) = zz_rate_uzoo_R * tgfunc(ji,jj,jk) &
                    * zz_Z(ji,jj,jk) * max(0._wp, zz_uZoo_graz_Z(ji,jj,jk))
                
                zz_MicroZoo_eat = (zz_uZoo_mort_micro(ji,jj,jk) + zz_uZoo_mort_nano(ji,jj,jk) + zz_uZoo_mort_pico(ji,jj,jk) &
                     + zz_uZoo_graz_PON(ji,jj,jk) + zz_uZoo_graz_Z(ji,jj,jk)) * zz_rate_uZoo_eff
                tra(ji,jj,jk,jpzoo) = tra(ji,jj,jk,jpzoo) + (zz_MicroZoo_eat - zz_uZoo_graz_Z(ji,jj,jk)) * rfact2
                tra(ji,jj,jk,jpdia) = tra(ji,jj,jk,jpdia) - zz_uZoo_mort_micro(ji,jj,jk) * rfact2
                tra(ji,jj,jk,jpphy) = tra(ji,jj,jk,jpphy) - zz_uZoo_mort_pico(ji,jj,jk) * rfact2
                tra(ji,jj,jk,jpmes) = tra(ji,jj,jk,jpmes) - zz_uZoo_mort_nano(ji,jj,jk) * rfact2

                tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) + ( &
                     zz_frac_waste_PEZ_NH * zz_UZoo_graz_PON(ji,jj,jk) + &
                     zz_frac_waste_DEZ_NH * zz_uZoo_mort_micro(ji,jj,jk) + &
                     zz_frac_waste_NEZ_NH * zz_uZoo_mort_nano(ji,jj,jk) + &
                     zz_frac_waste_FEZ_NH * zz_uZoo_mort_pico(ji,jj,jk) + &
                     zz_frac_waste_ZEZ_NH * zz_uZoo_graz_Z(ji,jj,jk)) * (1-zz_rate_uZoo_eff) * rfact2
                tra(ji,jj,jk,jpdon) = tra(ji,jj,jk,jpdon) + ( &
                     zz_frac_waste_PEZ_DON * zz_UZoo_graz_PON(ji,jj,jk) + &
                     zz_frac_waste_DEZ_DON * zz_uZoo_mort_micro(ji,jj,jk) + &
                     zz_frac_waste_NEZ_DON * zz_uZoo_mort_nano(ji,jj,jk) + &
                     zz_frac_waste_FEZ_DON * zz_uZoo_mort_pico(ji,jj,jk) + &
                     zz_frac_waste_ZEZ_DON * zz_uZoo_graz_Z(ji,jj,jk)) * (1-zz_rate_uZoo_eff) * rfact2
                tra(ji,jj,jk,jppon) = tra(ji,jj,jk,jppon) + ( &
                     zz_frac_waste_PEZ_PON * zz_UZoo_graz_PON(ji,jj,jk) + &
                     zz_frac_waste_DEZ_PON * zz_uZoo_mort_micro(ji,jj,jk) + &
                     zz_frac_waste_NEZ_PON * zz_uZoo_mort_nano(ji,jj,jk) + &
                     zz_frac_waste_FEZ_PON * zz_uZoo_mort_pico(ji,jj,jk) + &
                     zz_frac_waste_ZEZ_PON * zz_uZoo_graz_Z(ji,jj,jk))  * &
                     (1-zz_rate_uZoo_eff) * rfact2 - zz_uZoo_graz_PON(ji,jj,jk) * rfact2
                tra(ji,jj,jk,jpdsi) = tra(ji,jj,jk,jpdsi) + ( &
                     zz_frac_waste_DEZ_BSi * zz_uZoo_mort_micro(ji,jj,jk) * zz_rate_Si_ratio_diat + &
                     zz_frac_waste_NEZ_BSi * zz_uZoo_mort_nano(ji,jj,jk) * zz_rate_Si_ratio_myri + &
                     zz_frac_waste_FEZ_BSi * zz_uZoo_mort_pico(ji,jj,jk) * zz_rate_Si_ratio_nano &
                     ) * rfact2

            END DO
         END DO
      END DO
      
      CALL iom_put( "GRMICZDIAT", zz_uZoo_mort_micro)
      CALL iom_put( "GRMICZMRUB", zz_uZoo_mort_nano)
      CALL iom_put( "GRMICZPHY", zz_uZoo_mort_pico)
      CALL iom_put( "GRMICZPON", zz_uZoo_graz_PON)
      CALL iom_put( "GRMICZMICZ", zz_uZoo_graz_Z)

      CALL wrk_dealloc( jpi, jpj, jpk, zz_P_diat, zz_P_myri, zz_P_nano, zz_Z, zz_D_PON, zz_uZoo_mort_micro )
      CALL wrk_dealloc( jpi, jpj, jpk, zz_uZoo_mort_nano, zz_uZoo_mort_pico, zz_uZoo_graz_PON, zz_uZoo_graz_Z )
      !
      IF( nn_timing == 1 )  CALL timing_stop('p4z_micro')
      !
   END SUBROUTINE p4z_micro


   SUBROUTINE p4z_micro_init

      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_micro_init  ***
      !!
      !! ** Purpose :   Initialization of microzooplankton parameters
      !!
      !! ** Method  :   Read the nampiszoo namelist and check the parameters
      !!                called at the first timestep (nittrc000)
      !!
      !! ** input   :   Namelist nampiszoo
      !!
      !!----------------------------------------------------------------------

      NAMELIST/nampiszoo/ zz_rate_uzoo_Rm, zz_rate_uzoo_excr, zz_frac_waste_ZNM_NH, zz_frac_waste_ZEX_NH, &
         &   zz_rate_uzoo_PredSlope, zz_rate_uzoo_HalfSat, zz_rate_uzoo_MicroPref, zz_rate_uzoo_MicroPredslope, &
         &   zz_rate_uzoo_MicroHalfSat, zz_rate_uzoo_NanoPref, zz_rate_uzoo_NanoPredslope, zz_rate_uzoo_NanoHalfSat, &
         &   zz_rate_uzoo_PicoPref, zz_rate_uzoo_PicoPredslope, zz_rate_uzoo_PicoHalfSat, zz_rate_uzoo_PON_Pref, &
         &   zz_rate_uzoo_PON_Predslope, zz_rate_uzoo_PON_HalfSat, zz_rate_uzoo_Z_Pref, zz_rate_uzoo_Z_Predslope, &
         &   zz_rate_uzoo_Z_HalfSat, zz_rate_uzoo_R, zz_frac_waste_PEZ_NH, zz_frac_waste_DEZ_NH, zz_frac_waste_NEZ_NH, &
         &   zz_frac_waste_FEZ_NH, zz_frac_waste_ZEZ_NH, zz_rate_uZoo_eff, zz_frac_waste_ZNM_DON, zz_frac_waste_ZEX_DON, &
         &   zz_frac_waste_DEZ_DON, zz_frac_waste_FEZ_DON, zz_frac_waste_NEZ_DON, zz_frac_waste_PEZ_DON, &
         &   zz_frac_waste_ZEZ_DON, zz_frac_waste_ZNM_PON, zz_frac_waste_ZEX_PON, zz_frac_waste_DEZ_PON, &
         &   zz_frac_waste_FEZ_PON, zz_frac_waste_NEZ_PON, zz_frac_waste_PEZ_PON, zz_frac_waste_ZEZ_PON, &
         &   zz_frac_waste_ZNM_BSi, zz_frac_waste_ZEX_Bsi, zz_frac_waste_DEZ_Bsi, zz_frac_waste_FEZ_BSi, &
         &   zz_frac_waste_NEZ_BSi, zz_frac_waste_PEZ_BSi, zz_frac_waste_ZEZ_BSi 
      INTEGER :: ios                 ! Local integer output status for namelist read

      REWIND( numnatp_ref )              ! Namelist nampiszoo in reference namelist : microzooplankton
      READ  ( numnatp_ref, nampiszoo, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampiszoo in reference namelist', lwp )

      REWIND( numnatp_cfg )              ! Namelist nampiszoo in configuration namelist : microzooplankton
      READ  ( numnatp_cfg, nampiszoo, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampiszoo in configuration namelist', lwp )
      IF(lwm) WRITE ( numonp, nampiszoo )

      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for microzo, nampiszoo'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '  uM N / s microzo natural mort @ 20 deg C            zz_rate_uzoo_Rm=', zz_rate_uzoo_Rm
         WRITE(numout,*) '  uM N / s microzo excretion rate @ 20 deg C        zz_rate_uzoo_excr=', zz_rate_uzoo_excr
         WRITE(numout,*) '                                                 zz_frac_waste_ZNM_NH=', zz_frac_waste_ZNM_NH
         WRITE(numout,*) '                                                 zz_frac_waste_ZEX_NH=', zz_frac_waste_ZEX_NH
         WRITE(numout,*) '  uM N microzo total grazing limit             zz_rate_uzoo_PredSlope=', zz_rate_uzoo_PredSlope
         WRITE(numout,*) '  uM N microzo total grazing half sat const      zz_rate_uzoo_HalfSat=', zz_rate_uzoo_HalfSat
         WRITE(numout,*) '  microzo preference for micro-phyto           zz_rate_uzoo_MicroPref=', zz_rate_uzoo_MicroPref
         WRITE(numout,*) '  uM N miczo micphyt graz half sat const  zz_rate_uzoo_MicroPredslope=', zz_rate_uzoo_MicroPredslope
         WRITE(numout,*) '  uM N miczo micphyt grazing half sat const zz_rate_uzoo_MicroHalfSat=', zz_rate_uzoo_MicroHalfSat
         WRITE(numout,*) '  microzo preference for nanophyto              zz_rate_uzoo_NanoPref=', zz_rate_uzoo_NanoPref
         WRITE(numout,*) '  uM N microzo nanophyt grazing limit      zz_rate_uzoo_NanoPredslope=', zz_rate_uzoo_NanoPredslope
         WRITE(numout,*) '  uM N microzo nanophyt graz half sat const  zz_rate_uzoo_NanoHalfSat=', zz_rate_uzoo_NanoHalfSat
         WRITE(numout,*) '  microzo preference for picophyto              zz_rate_uzoo_PicoPref=', zz_rate_uzoo_PicoPref
         WRITE(numout,*) '  uM N microzo picophyt grazing limit      zz_rate_uzoo_PicoPredslope=', zz_rate_uzoo_PicoPredslope
         WRITE(numout,*) '  uM N microzo picophyt graz half sat const  zz_rate_uzoo_PicoHalfSat=', zz_rate_uzoo_PicoHalfSat
         WRITE(numout,*) '  microzo preference for PON                    zz_rate_uzoo_PON_Pref=', zz_rate_uzoo_PON_Pref
         WRITE(numout,*) '  uM N microzo microzo grazing limit       zz_rate_uzoo_PON_Predslope=', zz_rate_uzoo_PON_Predslope
         WRITE(numout,*) '  uM N microzo PON grazing half sat const    zz_rate_uzoo_PON_HalfSat=', zz_rate_uzoo_PON_HalfSat
         WRITE(numout,*) '  microzo preference for microzo cannibalism      zz_rate_uzoo_Z_Pref=', zz_rate_uzoo_Z_Pref
         WRITE(numout,*) '  uM N microzo microzo grazing limit         zz_rate_uzoo_Z_Predslope=', zz_rate_uzoo_Z_Predslope
         WRITE(numout,*) '  uM N microzo microzo grazing half sat const  zz_rate_uzoo_Z_HalfSat=', zz_rate_uzoo_Z_HalfSat
         WRITE(numout,*) '  match SOG 0.600E-04  1/s micro-phyto max growth rate zz_rate_uzoo_R=', zz_rate_uzoo_R
         WRITE(numout,*) '  waste frac from microzo graz PON to NH         zz_frac_waste_PEZ_NH=', zz_frac_waste_PEZ_NH
         WRITE(numout,*) '  waste frac from microzo graz microphyt to NH   zz_frac_waste_DEZ_NH=', zz_frac_waste_DEZ_NH
         WRITE(numout,*) '  waste frac from microzo graz nanophyt to NH    zz_frac_waste_NEZ_NH=', zz_frac_waste_NEZ_NH
         WRITE(numout,*) '  waste frac from microzo graz picophyt to NH    zz_frac_waste_FEZ_NH=', zz_frac_waste_FEZ_NH
         WRITE(numout,*) '  waste frac from microzo excretion to NH        zz_frac_waste_ZEZ_NH=', zz_frac_waste_ZEZ_NH
         WRITE(numout,*) '                                                     zz_rate_uZoo_eff=', zz_rate_uZoo_eff
         WRITE(numout,*) '  waste frac from microzo natural mort to DON   zz_frac_waste_ZNM_DON=', zz_frac_waste_ZNM_DON
         WRITE(numout,*) '  waste frac from microzo excretion to DON      zz_frac_waste_ZEX_DON=', zz_frac_waste_ZEX_DON
         WRITE(numout,*) '  waste frac from microzo graz microphyt to DON zz_frac_waste_DEZ_DON=', zz_frac_waste_DEZ_DON
         WRITE(numout,*) '  waste frac from microzo graz picophyt to DON  zz_frac_waste_FEZ_DON=', zz_frac_waste_FEZ_DON
         WRITE(numout,*) '  waste frac from microzo graz nanophyt to DON  zz_frac_waste_NEZ_DON=', zz_frac_waste_NEZ_DON
         WRITE(numout,*) '  waste frac from microzo graz PON to DON       zz_frac_waste_PEZ_DON=', zz_frac_waste_PEZ_DON
         WRITE(numout,*) '  waste frac from microzo graz microzoo to DON  zz_frac_waste_ZEZ_DON=', zz_frac_waste_ZEZ_DON
         WRITE(numout,*) '  waste frac from microzo natural mort to PON   zz_frac_waste_ZNM_PON=', zz_frac_waste_ZNM_PON
         WRITE(numout,*) '  waste frac from microzo excretion to PON      zz_frac_waste_ZEX_PON=', zz_frac_waste_ZEX_PON
         WRITE(numout,*) '  waste frac from microzo graz microphyt to PON zz_frac_waste_DEZ_PON=', zz_frac_waste_DEZ_PON
         WRITE(numout,*) '  waste frac from microzo graz picophyt to PON  zz_frac_waste_FEZ_PON=', zz_frac_waste_FEZ_PON
         WRITE(numout,*) '  waste frac from microzo graz nanophyt to PON  zz_frac_waste_NEZ_PON=', zz_frac_waste_NEZ_PON
         WRITE(numout,*) '  waste frac from microzo graz PON to PON       zz_frac_waste_PEZ_PON=', zz_frac_waste_PEZ_PON
         WRITE(numout,*) '  waste frac from microzo graz microzoo to PON  zz_frac_waste_ZEZ_PON=', zz_frac_waste_ZEZ_PON
         WRITE(numout,*) '  waste frac from microzo natural mort to bSi   zz_frac_waste_ZNM_BSi=', zz_frac_waste_ZNM_BSi
         WRITE(numout,*) '  waste frac from microzo excretion to bSi      zz_frac_waste_ZEX_Bsi=', zz_frac_waste_ZEX_Bsi
         WRITE(numout,*) '  waste frac from microzo graz microphyt to Bsi zz_frac_waste_DEZ_Bsi=', zz_frac_waste_DEZ_Bsi
         WRITE(numout,*) '  waste frac from microzo graz picophyt to Bsi  zz_frac_waste_FEZ_BSi=', zz_frac_waste_FEZ_BSi
         WRITE(numout,*) '  waste frac from microzo graz nanophyt to Bsi  zz_frac_waste_NEZ_BSi=', zz_frac_waste_NEZ_BSi
         WRITE(numout,*) '  waste frac from microzo graz PON to Bsi       zz_frac_waste_PEZ_BSi=', zz_frac_waste_PEZ_BSi
         WRITE(numout,*) '  waste frac from microzo graz microzoo to Bsi  zz_frac_waste_ZEZ_BSi=', zz_frac_waste_ZEZ_BSi
         WRITE(numout,*) '  micro-phyto silicon/nitrogen ratio           zz_rate_micro_Si_ratio=', zz_rate_Si_ratio_diat
         WRITE(numout,*) '  nano-phyto silicon/nitrogen ratio             zz_rate_nano_Si_ratio=', zz_rate_Si_ratio_myri
         WRITE(numout,*) '  pico-phyto silicon/nitrogen ratio             zz_rate_pico_Si_ratio=', zz_rate_Si_ratio_nano
      ENDIF

   END SUBROUTINE p4z_micro_init

#else
   !!======================================================================
   !!  Dummy module :                                   No bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_micro                    ! Empty routine
   END SUBROUTINE p4z_micro
#endif 
   !!======================================================================
END MODULE p4zmicro
