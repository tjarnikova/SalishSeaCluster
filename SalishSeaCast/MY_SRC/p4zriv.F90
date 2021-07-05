MODULE p4zriv
   !!======================================================================
   !!                         ***  MODULE p4riv  ***
   !! TOP :   SMELT Compute input of nutrients/tracers from rivers
   !!======================================================================
   !! History :   2016 !  (E. Olson) added module
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   'key_pisces'                              SMELT / PISCES bio-model
   !!----------------------------------------------------------------------
   !!   p4z_sed        :  Compute loss of organic matter in the sediments
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  Source Minus Sink variables
   !USE p4zsink         !  vertical flux of particulate matter due to sinking
   USE p4zopt          !  optical model
   USE p4zlim          !  Co-limitations of differents nutrients
   USE p4zsbc          !  External source of nutrients 
   USE p4zint          !  interpolation and computation of various fields
   USE iom             !  I/O manager
   USE prtctl_trc      !  print control for debugging

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_riv  
   PUBLIC   p4z_riv_init ! call from trcini_pisces.F90
   REAL(wp) :: Wd
   REAL(wp) :: alpha_C
   REAL(wp) :: dpers
   REAL(wp), PUBLIC :: wsink

   REAL(wp) :: r1_rday                  !: inverse of rday
   

   !!* Substitution
#  include "top_substitute.h90"

CONTAINS

   SUBROUTINE p4z_riv( kt)
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_riv  ***
      !!
      !! ** Purpose : Compute input of nutrients/tracer through rivers
      !!              May need to balance this with loss terms in future?
      !!              Currently loss is through zoo closure & open boundaries
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt ! ocean time step
      INTEGER  ::   ji, jj, jk, jn
      REAL(wp) ::   znfact, zfact, zdep
      !
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_riv')   
         ! Add the external input of nutrients from river
      ! ----------------------------------------------------------
      IF( ln_river ) THEN
         zfact = 0.5_wp
         DO jn = jp_pcs0, jp_pcs1
           DO jj = 1, jpj
             DO ji = 1, jpi
               IF( rnf(ji,jj) /= 0._wp ) THEN
                 zdep=zfact / h_rnf(ji,jj)
                 DO jk = 1, nk_rnf(ji,jj)
                   tra(ji,jj,jk,jn) = tra(ji,jj,jk,jn) &
                    &                + (rnf_pis_b(ji,jj,jn) + rnf_pis(ji,jj,jn) ) *zdep
                 ENDDO
               ENDIF
             ENDDO
           ENDDO
         ENDDO
      ENDIF
      IF( ln_turb ) THEN
         tra(:,:,:,jpriv) = -1_wp*alpha_C*trb(:,:,:,jpriv)*trb(:,:,:,jpriv)*dpers*rfact2
         DO jk = 1, jpk
            WHERE (sf_turb(1)%fnow(:,:,1) >= 0._wp)
               tra(:,:,jk,jpriv) = 0.0_wp 
               trn(:,:,jk,jpriv) = sf_turb(1)%fnow(:,:,1)
            END WHERE
         END DO
      ENDIF
   END SUBROUTINE p4z_riv

   SUBROUTINE p4z_riv_init
      !!-----------------------------------------------------------------------
      !!                   *** ROUTINE p4z_riv_init ***
      !! ** Purpose :  Initilization of turbidity decay timescales
      !!-----------------------------------------------------------------------
      !
      INTEGER :: ios       ! Local integer output status for namelist read
      INTEGER :: jn
      NAMELIST/nampisriv/ Wd, alpha_C 
      !!-----------------------------------------------------------------------

      IF( nn_timing ==1) CALL timing_start('p4z_riv_init')

      REWIND( numnatp_ref)            ! Namelist nampisriv in reference namelist
      READ ( numnatp_ref, nampisriv, IOSTAT = ios, ERR = 901)
901   IF( ios/= 0 ) CALL ctl_nam (ios, 'nampisriv in reference namelist', lwp )

      REWIND( numnatp_cfg)            ! Namelist nampisriv in configuration namelist
      READ ( numnatp_cfg, nampisriv, IOSTAT = ios, ERR = 902 )
902   IF( ios/= 0 ) CALL ctl_nam ( ios , 'nampisriv in configuration namelist', lwp )
      IF(lwm) WRITE ( numonp, nampisriv )

      dpers=1.0_wp/rday
      wsink = Wd * dpers

      IF(lwp) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) ' namelist : nampisriv '
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~ '
         WRITE(numout,*) '     Wd m/d, alpha_C 1/(d NTU)'
         WRITE(numout,*) '      ', Wd, alpha_C
       ENDIF

       IF( nn_timing == 1 ) CALL timing_stop('p4z_riv_init')
       !
    END SUBROUTINE p4z_riv_init
 
#endif

END MODULE p4zriv
