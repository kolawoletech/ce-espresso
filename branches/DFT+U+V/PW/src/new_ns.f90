!
! Copyright (C) 2001-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE new_ns(ns)
  !-----------------------------------------------------------------------
  !
  ! This routine computes the new value for ns (the occupation numbers of
  ! ortogonalized atomic wfcs).
  ! These quantities are defined as follows: ns_{I,s,m1,m2} = \sum_{k,v}
  ! f_{kv} <\fi^{at}_{I,m1}|\psi_{k,v,s}><\psi_{k,v,s}|\fi^{at}_{I,m2}>
  !
  USE io_global,            ONLY : stdout
  USE kinds,                ONLY : DP
  USE parameters,           ONLY : lmaxx
  USE ions_base,            ONLY : nat, ityp
  USE basis,                ONLY : natomwfc
  USE klist,                ONLY : nks, ngk
  USE ldaU,                 ONLY : oatwfc, swfcatom, ilm
  USE symm_base,            ONLY : d1, d2, d3
  USE lsda_mod,             ONLY : lsda, current_spin, nspin, isk
  USE symm_base,            ONLY : nsym, irt
  USE wvfct,                ONLY : nbnd, npw, npwx, igk, wg
  USE control_flags,        ONLY : gamma_only
  USE wavefunctions_module, ONLY : evc
  USE io_files,             ONLY : iunigk, nwordwfc, iunwfc, nwordatwfc, iunsat
  USE buffers,              ONLY : get_buffer
  USE mp_global,            ONLY : intra_pool_comm, inter_pool_comm
  USE mp,                   ONLY : mp_sum
  USE becmod,               ONLY : bec_type, calbec, &
                                   allocate_bec_type, deallocate_bec_type

  IMPLICIT NONE
  !
  REAL(DP), INTENT(OUT) :: ns((lmaxx+1)**2,(lmaxx+1)**2,nspin,nat)
  !
  TYPE (bec_type) :: proj     ! proj(natomwfc,nbnd)
  INTEGER :: ik, ibnd, is, i, na, nb, nt, isym, m1, m2, m0, m00, ldim
  ! counter on k points
  !    "    "  bands
  !    "    "  spins
  ! in the natomwfc ordering
  REAL(DP) , ALLOCATABLE :: nr (:,:,:,:)
  REAL(DP) :: psum
  REAL(DP), EXTERNAL :: ddot
  COMPLEX(DP) :: zdotc
  COMPLEX(DP) , ALLOCATABLE :: proj(:,:)

  CALL start_clock('new_ns')
  ldim = (lmaxx+1)**2
  ALLOCATE( nr(ldim,ldim,nspin,nat) )
  ldim = 2 * Hubbard_lmax + 1
  ALLOCATE( nr(ldim,ldim,nspin,nat) )  
  CALL allocate_bec_type ( natomwfc, nbnd, proj ) 
  !
  ! D_Sl for l=1, l=2 and l=3 are already initialized, for l=0 D_S0 is 1
  !
  ! Offset of atomic wavefunctions initialized in setup and stored in oatwfc
  !
  nr (:,:,:,:) = 0.d0
  ns (:,:,:,:) = 0.d0
  !
  !    we start a loop on k points
  !
  IF (nks > 1) REWIND (iunigk)
  DO ik = 1, nks
     current_spin = 1
     IF (lsda) current_spin = isk(ik)
     npw = ngk (ik)
     IF (nks > 1) THEN
        READ (iunigk) igk
        CALL get_buffer  (evc, nwordwfc, iunwfc, ik)
     END IF
     CALL davcio (swfcatom, nwordatwfc, iunsat, ik, - 1)
     !
     ! make the projection
     !
     CALL calbec ( npw, swfcatom, evc, proj )
     !
     ! compute the occupation numbers (the quantities n(m1,m2)) of the
     ! atomic orbitals
     !
     DO na = 1, nat  
        DO il = 0, lmaxx
          if (oatwfc(na,il) == -1) cycle
          DO m1 = 1, 2*il+1  
            DO m2 = m1, 2*il+1
              ilm1 = ilm(il,m1)
              ilm2 = ilm(il,m2)
              IF ( gamma_only ) THEN
                DO ibnd = 1, nbnd  
                  nr(ilm1,ilm2,current_spin,na) = nr(ilm1,ilm2,current_spin,na) + &
                    proj%r(oatwfc(na,il)+m2,ibnd) * &
                    proj%r(oatwfc(na,il)+m1,ibnd) * wg(ibnd,ik)
                ENDDO
              ELSE
                DO ibnd = 1, nbnd  
                  nr(ilm1,ilm2,current_spin,na) = nr(ilm1,ilm2,current_spin,na) + &
                    DBLE( proj%k(oatwfc(na,il)+m2,ibnd)) * &
                    CONJG(proj%k(oatwfc(na,il)+m1,ibnd)) * wg(ibnd,ik
                ENDDO
              ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
     ENDDO

! on k-points
  ENDDO
  CALL deallocate_bec_type (proj) 
#ifdef __PARA
  CALL mp_sum( nr, inter_pool_comm )
#endif
  IF (nspin.EQ.1) nr = 0.5d0 * nr
  !
  ! impose hermiticity of n_{m1,m2}
  !
  DO na = 1, nat  
    DO is = 1, nspin  
      DO il = 0, lmaxx  
        DO m1 = 1, 2*il+1  
          DO m2 = m1+1, 2*il+1
            ilm1 = ilm(il,m1)
            ilm2 = ilm(il,m2)
            nr(ilm2, ilm1,is, na) = nr(ilm1, ilm2, is, na)  
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  ! symmetrize the quantities nr -> ns
  DO na = 1, nat  
    DO is = 1, nspin  
      DO il = 0, lmaxx  
        if (oatwfc(na,il) == -1) cycle
        DO m1 = 1, 2*il+1  
          DO m2 = 1, 2*il+1
            ilm1 = ilm(il,m1)
            ilm2 = ilm(il,m2)

            DO isym = 1, nsym  
              nb = irt (isym, na)  
              DO m0 = 1, 2*il+1
                DO m00 = 1, 2*il+1
                  IF (il.EQ.0) THEN
                     ns(ilm1,ilm2,is,na) = ns(ilm1,ilm2,is,na) +  &
                                   nr(il*il+m0,il*il+m00,is,nb) / nsym
                  ELSE IF (il.EQ.1) THEN
                     ns(ilm1,ilm2,is,na) = ns(ilm1,ilm2,is,na) +  &
                                   d1(m0 ,m1,isym) * nr(ilm(il,m0),ilm(il,m00),is,nb) * &
                                   d1(m00,m2,isym) / nsym
                  ELSE IF (il.EQ.2) THEN
                     ns(ilm1,ilm2,is,na) = ns(ilm1,ilm2,is,na) +  &
                                   d2(m0 ,m1,isym) * nr(ilm(il,m0),ilm(il,m00),is,nb) * &
                                   d2(m00,m2,isym) / nsym
                  ELSE IF (il.EQ.3) THEN
                     ns(ilm1,ilm2,is,na) = ns(ilm1,ilm2,is,na) +  &
                                   d3(m0 ,m1,isym) * nr(ilm(il,m0),ilm(il,m00),is,nb) * &
                                   d3(m00,m2,isym) / nsym
                  ELSE
                     CALL errore ('new_ns', 'angular momentum not implemented', il)
                  END IF
                ENDDO
              ENDDO
            ENDDO

          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  ! Now we make the matrix ns(m1,m2) strictly hermitean
  DO na = 1, nat  
    DO is = 1, nspin  
      DO il = 0, lmaxx  
        DO m1 = 1, 2*il+1  
          DO m2 = m1, 2*il+1
            ilm1 = ilm(il,m1)
            ilm2 = ilm(il,m2)
            psum = ABS ( ns(ilm1,ilm2,is,na) - ns(ilm2,ilm1,is,na) )  
            IF (psum.GT.1.d-10) THEN  
              WRITE( stdout, * ) na, is, il, m1, m2  
              WRITE( stdout, * ) ns (ilm1, ilm2, is, na)  
              WRITE( stdout, * ) ns (ilm2, ilm1, is, na)  
              CALL errore ('new_ns', 'non hermitean matrix', 1)  
            ELSE  
              ns(ilm1,ilm2,is,na) = 0.5d0 * (ns(ilm1,ilm2,is,na) + ns(ilm2,ilm1,is,na) )
              ns(ilm2,ilm1,is,na) = ns(ilm1,ilm2,is,na)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  DEALLOCATE ( nr )
  CALL stop_clock('new_ns')

  RETURN

END SUBROUTINE new_ns