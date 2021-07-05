MODULE p4zcar
   !!======================================================================
   !!                         ***  MODULE p4zrem  ***
   !! TOP :   SMELT Compute remineralization/dissolution of organic compounds
   !!=========================================================================
   !! History :   2014 (E. Olson) adapted from 1-d SOG and PISCES
   !!----------------------------------------------------------------------
#if defined key_skog
   !!----------------------------------------------------------------------
   !!   'key_top'       and                                      TOP models
   !!   'key_pisces'                           SMELT / PISCES bio-model
   !!----------------------------------------------------------------------
   !!   p4z_rem       :  Compute remineralization/dissolution of organic compounds
   !!   p4z_rem_init  :  Initialisation of parameters for remineralisation
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  SMELT Source Minus Sink variables
   USE p4zopt          !  optical model
   !USE p4zprod ! ONLY: zz_remin_NH         !  Growth rate of the 2 phyto groups
   USE p4zint          !  interpolation and computation of various fields
   USE prtctl_trc      !  print control for debugging
   USE iom             !  I/O manager
   USE p4zrem, ONLY : zz_remin_NH, zz_remin_D_DON, zz_remin_D_PON
   USE par_oce, ONLY : jp_tem
   IMPLICIT NONE
   PRIVATE
   
   PUBLIC   p4z_car_init    ! called in trcini_pisces.F90
   PUBLIC   p4z_car         ! called in p4zprod
   PUBLIC   p4z_ta          ! TJSJ, p4zprod
   PUBLIC   p4z_oxy         ! TJSJ d
   REAL(wp) ::   zz_redfield_c_n, zz_redfield_p_n, zz_redfield_o_n !zz_remin_nh  ! redfield ratio of 106:16     	  

   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: p4zrem.F90 3160 2011-11-20 14:27:18Z cetlod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_car( kt, knt, zz_uptake_NO_diat, zz_uptake_NO_nano, zz_uptake_NO_myri, &
        zz_uptake_NH_diat, zz_uptake_NH_nano, zz_uptake_NH_myri, &
        zz_uptake_PC_diat, zz_uptake_PC_nano, zz_uptake_PC_myri)
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_car  ***
      !!
      !! ** Purpose :   Compute change in DIC
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER  ::   ji, jj, jk
      INTEGER, INTENT(in) ::   kt, knt ! ocean time step
      REAL(wp), INTENT(in)   :: zz_uptake_NO_diat(:,:,:), zz_uptake_NO_nano(:,:,:), zz_uptake_NO_myri(:,:,:), &
      zz_uptake_NH_diat(:,:,:), zz_uptake_NH_nano(:,:,:), zz_uptake_NH_myri(:,:,:), &
      zz_uptake_PC_diat(:,:,:), zz_uptake_PC_nano(:,:,:), zz_uptake_PC_myri(:,:,:)

      REAL(wp) :: zz_loc_remin_NH4
      IF( nn_timing == 1 )  CALL timing_start('p4z_car')

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
 

                zz_loc_remin_NH4 = (zz_remin_D_DON*trb(ji,jj,jk,jpdon) + zz_remin_D_PON*trb(ji,jj,jk,jppon)) *tgfunc(ji,jj,jk)
               
               tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) + ( zz_loc_remin_NH4 - &
               (zz_uptake_NO_diat(ji,jj,jk)+zz_uptake_NH_diat(ji,jj,jk)+zz_uptake_PC_diat(ji,jj,jk) + &
               zz_uptake_NO_nano(ji,jj,jk)+zz_uptake_NH_nano(ji,jj,jk)+zz_uptake_PC_nano(ji,jj,jk) + &
               zz_uptake_NO_myri(ji,jj,jk)+zz_uptake_NH_myri(ji,jj,jk)+zz_uptake_PC_myri(ji,jj,jk)) &
                ) * zz_redfield_c_n * rfact2  
            END DO
         END DO
      END DO
      IF( nn_timing == 1 )  CALL timing_stop('p4z_car')
      

   END SUBROUTINE p4z_car


   SUBROUTINE p4z_oxy( kt, knt, zz_uptake_NO_diat, zz_uptake_NO_nano, zz_uptake_NO_myri, &
                        zz_uptake_NH_diat, zz_uptake_NH_nano, zz_uptake_NH_myri)
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_oxy  ***
      !!
      !! ** Purpose :   Compute change in O2
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER  ::   ji, jj, jk
      INTEGER, INTENT(in) ::   kt, knt ! ocean time step
      REAL(wp), INTENT(in)   :: zz_uptake_NO_diat(:,:,:), zz_uptake_NO_nano(:,:,:), zz_uptake_NO_myri(:,:,:), &
                                zz_uptake_NH_diat(:,:,:), zz_uptake_NH_nano(:,:,:), zz_uptake_NH_myri(:,:,:)

      REAL(wp) :: zz_loc_remin_NH4
      REAL(wp) :: zz_loc_oxy_NH4
      REAL(wp) :: zz_loc_temp_Q10
      IF( nn_timing == 1 )  CALL timing_start('p4z_oxy')

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi


                zz_loc_remin_NH4 = (zz_remin_D_DON*trb(ji,jj,jk,jpdon) + zz_remin_D_PON*trb(ji,jj,jk,jppon)) *tgfunc(ji,jj,jk)
                zz_loc_oxy_NH4 = zz_remin_NH * trb(ji,jj,jk,jpnh4) ** 2 * tgfunc(ji,jj,jk)

! update o2 arrays
               tra(ji,jj,jk,jpo2) = tra(ji,jj,jk, jpo2) + ( & 
               (zz_uptake_NO_diat(ji,jj,jk)+zz_uptake_NH_diat(ji,jj,jk) + &
               zz_uptake_NO_nano(ji,jj,jk)+zz_uptake_NH_nano(ji,jj,jk) + &
               zz_uptake_NO_myri(ji,jj,jk)+zz_uptake_NH_myri(ji,jj,jk)) &
                - zz_loc_remin_NH4 - zz_loc_oxy_NH4 ) * zz_redfield_o_n * rfact2

            END DO
         END DO
      END DO
      IF( nn_timing == 1 )  CALL timing_stop('p4z_oxy')


   END SUBROUTINE p4z_oxy


  SUBROUTINE p4z_car_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_car_init  ***
      !!
      !! ** Purpose :   Initialization of carbonate parameters
      !!
      !! ** Method  :   Read the nampiscarbon namelist and check the parameters
      !!      called at the first timestep
      !!
      !! ** input   :   Namelist nampiscarbon
      !!
      !!----------------------------------------------------------------------
      NAMELIST/nampiscarbon/ zz_redfield_c_n, zz_redfield_p_n, zz_redfield_o_n
     INTEGER :: ios                 ! Local integer output status for namelist read

      REWIND( numnatp_ref )      ! Namelist nampiscarbon in reference namelist : (but I dont think there is one TJSJ)
      READ  ( numnatp_ref, nampiscarbon, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampiscarbon in reference namelist', lwp )

      REWIND( numnatp_cfg )      ! Namelist nampiscarbon in configuration namelist : carbonate chemistry
      READ  ( numnatp_cfg, nampiscarbon, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampiscarbon in configuration namelist', lwp )
      IF(lwm) WRITE ( numonp, nampiscarbon )

      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for carbonate chemistry, nampiscarbon'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) ' Redfield ratio of carbon to nitrogen    zz_redfield_c_n   =', zz_redfield_c_n
         WRITE(numout,*) ' Redfield ratio of phosphorus to nitrogen    zz_redfield_p_n  =', zz_redfield_p_n
         WRITE(numout,*) ' Redfield ratio of phosphorus to nitrogen    zz_redfield_o_n  =', zz_redfield_o_n

      ENDIF

   END SUBROUTINE p4z_car_init

   SUBROUTINE p4z_ta( kt, knt, zz_uptake_NO_diat, zz_uptake_NO_nano, zz_uptake_NO_myri, &
                        zz_uptake_NH_diat, zz_uptake_NH_nano, zz_uptake_NH_myri )

      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_car  ***
      !!
      !! ** Purpose :   Compute change in DIC
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      
      INTEGER  ::   ji, jj, jk
      
      INTEGER, INTENT(in) ::   kt, knt ! ocean time step
      REAL(wp), INTENT(in)   :: zz_uptake_NO_diat(:,:,:), zz_uptake_NO_nano(:,:,:), zz_uptake_NO_myri(:,:,:), &
                                zz_uptake_NH_diat(:,:,:), zz_uptake_NH_nano(:,:,:), zz_uptake_NH_myri(:,:,:)
      REAL(wp) ::   zz_loc_remin_NH4, zz_loc_oxy_NH4

     IF( nn_timing == 1 )  CALL timing_start('p4z_ta')

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi

                 zz_loc_remin_NH4 = (zz_remin_D_DON*trb(ji,jj,jk,jpdon) + zz_remin_D_PON*trb(ji,jj,jk,jppon)) *tgfunc(ji,jj,jk)
                 zz_loc_oxy_NH4 = zz_remin_NH * trb(ji,jj,jk,jpnh4) ** 2 * tgfunc(ji,jj,jk)

                tra(ji,jj,jk,jpta) =tra(ji,jj,jk,jpta)+ ((zz_redfield_p_n + 1.0) * &
                (zz_uptake_NO_diat(ji,jj,jk)+ &
                zz_uptake_NO_myri(ji,jj,jk)+&
                zz_uptake_NO_nano(ji,jj,jk))+ &  
                (zz_redfield_p_n - 1.0) * &
                (zz_uptake_NH_diat(ji,jj,jk)+ &
                 zz_uptake_NH_myri(ji,jj,jk)+&
                 zz_uptake_NH_nano(ji,jj,jk) &
                 - zz_loc_remin_NH4)  - 2 * zz_loc_oxy_NH4) * rfact2

                 
            END DO
         END DO
      END DO
     
      IF( nn_timing == 1 )  CALL timing_stop('p4z_ta')
      

   END SUBROUTINE p4z_ta


#else
   !!======================================================================
   !!  Dummy module :                                   No bio-model
   !!======================================================================

CONTAINS
   SUBROUTINE p4z_car                    ! Empty routine
   END SUBROUTINE p4z_car

   SUBROUTINE p4z_car_init                    ! Empty routine
   END SUBROUTINE p4z_car_init
 

   SUBROUTINE p4z_ta
   END SUBROUTINE p4z_ta

   SUBROUTINE p4z_oxy
   END SUBROUTINE p4z_oxy
#endif 

   !!======================================================================
END MODULE p4zcar
