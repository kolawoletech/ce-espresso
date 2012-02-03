!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE stop_lr( )
  !----------------------------------------------------------------------------
  !
  ! ... Synchronize processes before stopping.
  !
  ! Modified by O. Baris Malcioglu (2009)
  USE kinds, ONLY : DP
  USE mp,    ONLY : mp_end, mp_barrier
  !
  USE parallel_include
  USE lr_variables,         ONLY : n_ipol, LR_polarization, beta_store
  USE lr_variables,         ONLY :  gamma_store, zeta_store, norm0, rho_1_tot
  USE lr_variables,         ONLY : lr_verbosity, itermax, bgz_suffix
  USE io_global,            ONLY : ionode
  USE io_files,             ONLY : tmp_dir, prefix
  USE io_global,            ONLY : stdout
  USE environment,          ONLY : environment_end
  ! For gaussian cube file
  USE ions_base,  ONLY : nat, ityp, atm, ntyp => nsp, tau
  USE cell_base,  ONLY : celldm, at, bg
  !
  IMPLICIT NONE
  !
  CHARACTER(len=6), EXTERNAL :: int_to_char
  !
  CHARACTER(len=256) :: filename
  !
  INTEGER :: ip,i,j
  !
  !

  IF (lr_verbosity > 5) THEN
    WRITE(stdout,'("<stop_lr>")')
  ENDIF
  ! I write the beta gamma and z coefficents to output directory for
  ! easier post processing. These can also be read from the output log file
#ifdef __MPI
  IF (ionode) THEN
#endif
  !
  DO ip=1,n_ipol
   IF (n_ipol==3) filename = trim(prefix) // trim(bgz_suffix) // trim(int_to_char(ip))
   IF (n_ipol==1) filename = trim(prefix) // trim(bgz_suffix) // trim(int_to_char(LR_polarization))
   filename = trim(tmp_dir) // trim(filename)
  !
  !
  OPEN (158, file = filename, form = 'formatted', status = 'replace')
  !
  WRITE(158,*) itermax
  !
  WRITE(158,*) norm0(ip)
  !
  DO i=1,itermax
     !
     WRITE(158,*) beta_store(ip,i)
     WRITE(158,*) gamma_store(ip,i)
     !This is absolutely necessary for cross platform compatibilty
     DO j=1,n_ipol
      WRITE(158,*) zeta_store (ip,j,i)
     ENDDO
     !
  ENDDO
  !
  CLOSE(158)
  !
  ENDDO
#ifdef __MPI
  ENDIF
#endif
  !
  !   Deallocate lr variables
  !
  CALL lr_dealloc()
  !
  CALL environment_end('TDDFPT')
  !
  CALL mp_barrier()
  !
  CALL mp_end()
  !
#if defined (__T3E)
  !
  ! ... set streambuffers off
  !
  CALL set_d_stream( 0 )
  !
#endif
  !
  STOP
  !
  !
END SUBROUTINE stop_lr
