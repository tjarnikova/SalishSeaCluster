MODULE p4zsink
   !!======================================================================
   !!                         ***  MODULE p4zsink  ***
   !! TOP :  SMELT  vertical fluxes due to gravitational sinking
   !!======================================================================
   !! History :   2014 (E. Olson) adapted from 1-d SOG and PISCES
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   p4z_sink       :  Compute vertical flux of particulate matter due to gravitational sinking
   !!   p4z_sink_init  :  Unitialisation of sinking speed parameters
   !!   p4z_sink_alloc :  Allocate sinking speed variables
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  Source Minus Sink variables
   USE prtctl_trc      !  print control for debugging
   USE iom             !  I/O manager
   USE lib_mpp
   USE p4zprod, ONLY: zz_micro_Nlimit
   USE p4zriv, ONLY: wsink

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_sink         ! called in p4zbio.F90
   PUBLIC   p4z_sink_init    ! called in trcsms_pisces.F90

   REAL(wp), PUBLIC  ::   zz_w_sink_Pmicro_min  ! m/d  microphyto minimum sinking rate # Alain
   REAL(wp), PUBLIC  ::   zz_w_sink_Pmicro_max  ! m/d microphyto maximum sinking rate # Alain
   REAL(wp), PUBLIC  ::   zz_w_sink_D_PON       ! m/s PON detritus sinking rate # Jeffery quoting Dune and Bacon
   REAL(wp), PUBLIC  ::   zz_w_sink_D_bSi       !  m/s  biogenic silicon detritus sinking rate # match NO3 particles
   REAL(wp), PUBLIC  ::   zz_alpha_b_Si         !  fraction of sinking flux reflected back at bottom - Si
   REAL(wp), PUBLIC  ::   zz_alpha_b_D         !  fraction of sinking flux reflected back at bottom - diatoms
   REAL(wp), PUBLIC  ::   zz_alpha_b_N          !  fraction of sinking flux reflected back at bottom - N
   REAL(wp), PUBLIC  ::   zz_alpha_b_T          !  fraction of sinking flux reflected back at bottom - turbidity

   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: p4zsink.F90 3160 2011-11-20 14:27:18Z cetlod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_sink ( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_sink  ***
      !!
      !! ** Purpose :   Compute vertical flux of particulate matter due to 
      !!                gravitational sinking
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) :: kt, knt
      INTEGER  ::   ji, jj, jk
      !!---------------------------------------------------------------------
      !
      ! sinking defined by nutrient status
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zz_Pmicro_w_sink, zz_unvec, zz_Pmicro, zz_D_PON, &
                                                 zz_D_bSi, zz_Pmicro_new, zz_D_PON_new, &
                                                 zz_D_bSi_new, zz_TRA_new
    
      IF( nn_timing == 1 )  CALL timing_start('p4z_sink')



      CALL wrk_alloc(jpi, jpj, jpk, zz_Pmicro_w_sink, zz_unvec, zz_Pmicro, zz_D_PON )
      CALL wrk_alloc(jpi, jpj, jpk, zz_D_bSi, zz_Pmicro_new, zz_D_PON_new, zz_D_bSi_new)
      CALL wrk_alloc(jpi, jpj, jpk, zz_TRA_new)
      
      zz_Pmicro(:,:,:)=trb(:,:,:,jpdia)
      zz_D_PON(:,:,:) =trb(:,:,:,jppon)
      zz_D_bSi(:,:,:) =trb(:,:,:,jpdsi)
      zz_unvec(:,:,:) = 1.0_wp
      
      DO jk = 1,jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
              ! Calculate the sinking term for the quantities that sink
              ! to  calculate sinking at the interfaces used the values above
              zz_Pmicro_w_sink(ji,jj,jk) = zz_w_sink_Pmicro_min * zz_micro_Nlimit(ji,jj,jk) &
                  + zz_w_sink_Pmicro_max * (1.0_wp - zz_micro_Nlimit(ji,jj,jk))
            END DO
         END DO
      END DO
      
      ! limit the values of the sinking speeds to avoid numerical instabilities  
      !
      ! OA Below, this is garbage. the ideal would be to find a time-splitting 
      ! OA algorithm that does not increase the computing cost by too much
      ! OA In ROMS, I have included a time-splitting procedure. But it is 
      ! OA too expensive as the loop is computed globally. Thus, a small e3t
      ! OA at one place determines the number of subtimesteps globally
      ! OA AWFULLY EXPENSIVE !! Not able to find a better approach. Damned !!
      !DO jk = 1,jpkm1
      !   DO jj = 1, jpj
      !      DO ji = 1, jpi
      !         zwsmax = 0.5 * fse3t(ji,jj,jk) / xstep
      !         zz_Pmicro_w_sink(ji,jj,jk) = MIN( zz_Pmicro_w_sink(ji,jj,jk), zwsmax )
      !      END DO
      !   END DO
      !END DO
      !zz_w_sink_D_PON = MIN( zz_w_sink_D_PON, zwsmax )
      !zz_w_sink_D_bSi = MIN( zz_w_sink_D_bSi(ji,jj,jk), zwsmax )
     call zz_sinking_advection(zz_Pmicro, zz_Pmicro_w_sink, zz_Pmicro_new,zz_alpha_b_D)
     call zz_sinking_advection(zz_D_PON, zz_w_sink_D_PON*zz_unvec, zz_D_PON_new,zz_alpha_b_N)
     call zz_sinking_advection(zz_D_bSi, zz_w_sink_D_bSi*zz_unvec, zz_D_bSi_new,zz_alpha_b_Si)
     call zz_sinking_advection(trb(:,:,:,jptra), wsink*zz_unvec, zz_TRA_new,zz_alpha_b_T)
     WHERE (tmask==1)
         trb(:,:,:,jpdia)=zz_Pmicro_new
         trb(:,:,:,jppon)=zz_D_PON_new
         trb(:,:,:,jpdsi)=zz_D_bSi_new
         trb(:,:,:,jptra)=zz_TRA_new
     END WHERE 
     CALL wrk_dealloc(jpi, jpj, jpk, zz_Pmicro_w_sink, zz_unvec, zz_Pmicro, zz_D_PON )
     CALL wrk_dealloc(jpi, jpj, jpk, zz_D_bSi, zz_Pmicro_new, zz_D_PON_new, zz_D_bSi_new)
     CALL wrk_dealloc(jpi, jpj, jpk, zz_TRA_new)
     IF( nn_timing == 1 )  CALL timing_stop('p4z_sink')
      !
   END SUBROUTINE p4z_sink
   
     subroutine zz_sinking_advection(zz_qty, zz_w_sink, zz_qty_new, zz_alpha_b) !zz_dz, zz_dt, 
    ! Calculate the sinking term of the semi-implicit PDEs for the biology
    ! quantities that sink (some classes of plankton, and some of detritus).
    ! use upwind advection

    ! Arguments:
    !real(kind=wp), dimension(1:), intent(in) :: &
    !     zz_dz  ! Grid arrays
    !real(kind=wp) :: &
    !     zz_dt  ! Time step [s]
    INTEGER  ::   ji, jj, jk
    real(kind=wp), intent(in) :: &
         zz_alpha_b  ! bb reflection coefficient
    real(kind=wp), dimension(jpi,jpj,jpk), intent(in) :: &
         zz_qty  ! Profile array of sinking quantity
    real(kind=wp), dimension(jpi,jpj,jpk), intent(in) :: &
         zz_w_sink  ! Sinking velocity [m/s] on interfaces
    real(kind=wp), dimension(jpi,jpj,jpk), intent(out) :: &
         zz_qty_new  ! RHS term vector for semi-implicit diffusion/advection PDE

    ! local variables
    integer :: zz_i ! grid index 
    real(kind=wp), dimension(jpi,jpj,0:jpk) :: zz_flux ! positive is upward

    zz_flux(:,:,0) = 0._wp ! no flux in or out from above
    
      DO jk = 1,jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               if (zz_w_sink(ji,jj,jk).lt.0) then
                  write (*,*) "You have got biology sinking upward!"
                  call exit(1)
               else
                  zz_flux(ji,jj,jk) = - zz_w_sink(ji,jj,jk) * zz_qty(ji,jj,jk) * &
                          & ((1-zz_alpha_b)*tmask(ji,jj,jk)+zz_alpha_b*tmask(ji,jj,jk+1))
               endif
            END DO
         END DO
      END DO
      zz_flux(:,:,jpk)=0._wp
      !WRITE(numout,*) "sinkvel: ", zz_w_sink(3,3,:)
      !WRITE(numout,*) "zz_flux: ", zz_flux(3,3,:)
      !WRITE(numout,*) "zz_qty: ", zz_qty(3,3,:)
      DO jk = 1,jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
                zz_qty_new(ji,jj,jk) = zz_qty(ji,jj,jk) - rfact2 * (zz_flux(ji,jj,jk-1) - zz_flux(ji,jj,jk))/fse3t(ji,jj,jk)
            END DO
         END DO
      END DO
      !WRITE(numout,*) "flux: ", -rfact2 * (zz_flux(2,2,0:19) - zz_flux(2,2,1:20))/fse3t(2,2,1:20)
  end subroutine zz_sinking_advection


   SUBROUTINE p4z_sink_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_sink_init  ***
      !!----------------------------------------------------------------------
      INTEGER  ::   ios                 ! Local integer output status for namelist read
      !
      NAMELIST/nampissink/ zz_w_sink_Pmicro_min, zz_w_sink_Pmicro_max, zz_w_sink_D_PON, zz_w_sink_D_bSi, &
                           & zz_alpha_b_Si, zz_alpha_b_D, zz_alpha_b_N, zz_alpha_b_T


      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_sink_init')
      !

      REWIND( numnatp_ref )              ! Namelist nampissink in reference namelist : SOG sinking
      READ  ( numnatp_ref, nampissink, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampissink in reference namelist', lwp )

      REWIND( numnatp_cfg )              ! Namelist nampissink in configuration namelist : SOG sinking
      READ  ( numnatp_cfg, nampissink, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampissink in configuration namelist', lwp )
      IF(lwm) WRITE ( numonp, nampissink )

      zz_w_sink_Pmicro_min=zz_w_sink_Pmicro_min/rday
      zz_w_sink_Pmicro_max=zz_w_sink_Pmicro_max/rday

      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for sinking of biological material, nampissink'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '  m/d  microphyto minimum sinking rate zz_w_sink_Pmicro_min*86400 =', zz_w_sink_Pmicro_min*rday
         WRITE(numout,*) '  m/d microphyto maximum sinking rate  zz_w_sink_Pmicro_max*86400 =', zz_w_sink_Pmicro_max*rday
         WRITE(numout,*) '  m/s PON detritus sinking rate        zz_w_sink_D_PON =', zz_w_sink_D_PON
         WRITE(numout,*) '  m/s  bio si detritus sinking rate    zz_w_sink_D_bSi =', zz_w_sink_D_bSi
         WRITE(numout,*) '  BB reflection parameter Si           zz_alpha_b_Si =', zz_alpha_b_Si
         WRITE(numout,*) '  BB reflection parameter diatoms      zz_alpha_b_D =', zz_alpha_b_D
         WRITE(numout,*) '  BB reflection parameter N            zz_alpha_b_N =', zz_alpha_b_N
         WRITE(numout,*) '  BB reflection parameter turbidity    zz_alpha_b_T =', zz_alpha_b_T
      ENDIF


      
      !
      IF( nn_timing == 1 )  CALL timing_stop('p4z_sink_init')
      !
  END SUBROUTINE p4z_sink_init
   
#else
   !!======================================================================
   !!  Dummy module :                                   No bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_sink                    ! Empty routine
   END SUBROUTINE p4z_sink
#endif 
   !!======================================================================
END MODULE p4zsink
